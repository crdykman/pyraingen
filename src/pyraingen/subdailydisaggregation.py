import numpy as np
import random
from numba import njit

from .computesubdailykde import computeSubDailyKDE
from .getdaysperyear import getDaysPerYear
from .getwetstate import getWetState
from .getseasonfromday import getSeasonFromDay

from .global_ import missingDay
from .global_ import ndaysYear, ndaysYearLeap
from .global_ import idxfebTwentyNine
from .global_ import recordsPerDay
from .global_ import yesterday, today, tomorrow

# nogil=True matters: regionalisedsubdailysim dispatches this through
# joblib's threading backend, and a numba nopython function holds the GIL
# unless told otherwise, so without it the "parallel" simulations run one at a
# time. Safe here because param, nGoodDays and the three fragment arrays are
# only read, every write targets a function-local array, and numba keeps its
# random state per thread.
@njit(nogil=True)
def subDailyDisaggregation(targetDailyRain, param,
            nGoodDays, fragments, fragmentsState, fragmentsDailyDepth):
    """Sub-daily disaggregation based on
    Westra et al 2012 "Continuous rainfall simulation 1: A regionalised subdaily
    disaggregation approach.

    Parameters
    ----------
    targetDailyRain :  array
        Array of size nDaysKnown filled with daily depths.
    param : dict where keys are words and values are float
        Dictionary of the run parameters.
    nGoodDays : array
        Number of good days per year per day per season.
        A "good day" is one that is not of state bad and has
        a depth greater than the dryWetCutoff. 
    fragments : list
        A list of arrays (one for each season) containing a count  
        of the number of good fragments available for sampling.
    fragmentsState : list 
        A list of arrays (one for each season) containing the daily wetState sequences
        of the good fragments available for sampling.
    fragmentsDailyDepth : list 
        A list of arrays (one for each season) containing the daily depth sequences
        of the good fragments available for sampling.\n

    Returns
    ----------
    array
        Array of subdaily disaggregation data (recordsPerDay, nDays).
    """
    ## Step b) Loop Over The Days - This is what takes time
    # Allocate RAM
    nDays = targetDailyRain.size
    subDailySim = np.ones((recordsPerDay, nDays)) * missingDay

    # This is a variable to hold yesterdays, todays and tomorrows depths.
    workingDailyDepth = np.zeros(3)

    # Counter into the linear rainfall array second dimension
    idxDayLinear = 0

    # We know in advance that we want a maximum of param['maxNearNeighb']
    # fragments that we are going to randomly select from for any given day.
    # Therefore, we can save some RAM by only allocating for that and sorting
    # as we go.  Hence, the following potential array of fragments is of size:
    #   1) param.maxNearNeighb by;
    #   2) 3
    # with the columns in second dimension being:
    #   1) absolute difference of depth between potential day and target day
    #   2) index of the fragment into the day of the year, i.e. along the first
    #       dimension of nearbyPool.fragments
    #   3) index of the fragment into the number of good days, i.e. along the
    #       second dimension of nearbyPool.fragments.
    #
    # I initialise it high as we are looking for small differences
    potFrags = np.ones((int(param['maxNearNeighb']), 3)) * 1e9
    # Define some parameters to index into this array
    potFragsIdxAbsDiff = 0
    potFragsIdxDayYear = 1
    potFragsIdxDayGood = 2

    # These two are used in the insert function and are only computed once per
    # run:
    idxShiftOldEnd = potFrags[:,0].size -1   #double check
    idxShiftNewEnd = potFrags[:,0].size      #double check
    # With this array we will sort as we go and only keep fragments that we are
    # interested in.
    #
    # This is the array that I'll be holding the indices that I'll use to
    # slice through the first dimension of nearbyPool.fragments when searching
    # for the nearest neighbours.
    idxDayWin = np.zeros(((2*int(param['halfWinLen']))+1, 1))

    # Pre-compute the probability vector as it is only a function of
    # param['maxNearNeighb'].  This is Equation 2 in the paper.
    probSampleVecFull = computeSubDailyKDE(param['maxNearNeighb'])

    # Counter to track the number of days that we did not find a fragment.
    warningCounter = 0
    warningSum = 0
    check = 0

    # This double loop is analogous to Step 3 in the paper.
    for loopYear in range(int(param['simYearStart']), int(param['simYearEnd'])+1):
        nDaysCurrYear = getDaysPerYear(loopYear)
        for loopDay in range(ndaysYearLeap):
            # Extract some conveniance variables, first to hold some working rain
            # falls:
            #For GenSeqOpt < 3
            if loopDay < param['DayStart'] and loopYear == param['simYearStart']:
                idxDayLinear += 1
                continue
            if loopDay > param['DayEnd'] and loopYear == param['simYearEnd']:
                idxDayLinear += 1
                continue
            if (loopDay == 0 or loopDay == param['DayStart']) and loopYear == param['simYearStart']:
                # Similarly to above assume the first day was dry before we
                # started.
                workingDailyDepth[yesterday] = 0.0
                workingDailyDepth[today] = targetDailyRain[idxDayLinear]
                workingDailyDepth[tomorrow] = targetDailyRain[idxDayLinear+1]
            elif (loopDay == (ndaysYearLeap-1) or loopDay == param['DayEnd']) and loopYear ==param['simYearEnd']:
                # And assume the last day is also dry.
                #
                # The test is against ndaysYearLeap-1, not nDaysCurrYear-1: the
                # loop always runs over ndaysYearLeap and skips 29 February in
                # non-leap years, so the final day sits at index 365 either way.
                # Using nDaysCurrYear-1 made this branch fire a day early *and*
                # again at 365 in non-leap years, disaggregating the last two
                # days as if each were the end of the record.
                #
                # The window is shifted here as in every other branch, so that
                # today holds this day's depth rather than the previous day's.
                workingDailyDepth[yesterday] = workingDailyDepth[today]
                workingDailyDepth[today]     = workingDailyDepth[tomorrow]
                workingDailyDepth[tomorrow] = 0.0
            elif (loopDay == idxfebTwentyNine and
                nDaysCurrYear != ndaysYearLeap): 
                # If we are on the 60th day of the year but this is
                # not a leap year increment loopDay by one as that
                # is the index we are looping through our 2D array
                # year with.  Or, just continue to the next
                # iteration of the day loop
                #idxDayLinear += 1
                continue
            else:
                # shuffle the data through for this loop
                #idxDayLinear += 1
                workingDailyDepth[yesterday] = workingDailyDepth[today]
                workingDailyDepth[today]     = workingDailyDepth[tomorrow]
                workingDailyDepth[tomorrow]  = targetDailyRain[idxDayLinear + 1]
            # check if today is wet or dry?
            if workingDailyDepth[today] < param['dryWetCutoff']:
                # It is dry so no further work.
                subDailySim[:, idxDayLinear] = 0.0
                idxDayLinear += 1
                continue

            # Otherwise we have a wet day, so for the remainder of Step 3 we have
            # to:
            #   1) clasifiy today's wet state, eg WetWet, DryWet, etc, see above
            #   2) find which season we are in to confine the search.
            #   3) look for days within the moving window +-param.halfWinLen long,
            #       e.g. to the window itself is length=(2*param.halfWinLen+1).
            #
            # 1) Classify
            todayWetState = (
                 getWetState(workingDailyDepth, param['dryWetCutoff'])
                 )
            todaysRain = workingDailyDepth[today]

            # 2) Get the Season
            todaysSeason = int(getSeasonFromDay(loopDay+1, nDaysCurrYear))

            # 3) Window Search
            # Compute the indices of the moving window into the year days and
            # still handle leap years.  First generate the sequence as we only
            # need to, potentially, change it if this is a not a leap year
            idxDayWin = (
                np.arange((loopDay - param['halfWinLen']),(loopDay + param['halfWinLen']+1),1)
            )
            if nDaysCurrYear == ndaysYear:
                # Now we have to do a bit more work as this is NOT a leap
                # year...
                if loopDay < idxfebTwentyNine:
                    # This handles when our day of interest is coming up to Feb
                    # 29, so add one to the indices that are at Feb 29 or
                    # greater.
                    idxDayWin[idxDayWin >= idxfebTwentyNine] = (
                        idxDayWin[idxDayWin >= idxfebTwentyNine] + 1
                    )
                else:
                    # This is the alternate case where our day of interest is
                    # past Feb 29 and going higher, therefore it is the
                    # opposite of the alternate case.
                    idxDayWin[idxDayWin <= idxfebTwentyNine] = (
                        idxDayWin[idxDayWin <= idxfebTwentyNine] - 1
                    )
            
            # "unwrap" idxDay win to handle the wrap around at 31/12<->1/1
            idxDayWin[idxDayWin < 0] = idxDayWin[idxDayWin < 0] + nDaysCurrYear
            idxDayWin[idxDayWin > (nDaysCurrYear-1)] = (
                 idxDayWin[idxDayWin > (nDaysCurrYear-1)] - nDaysCurrYear
            )

            # Re-initialise (Why?)
            potFrags = np.ones((int(param['maxNearNeighb']), 3)) * 1e9
            fragFoundCounter = 0

            #BREAK POINT

            # Now with our indices computed loop and search 
            for loopSearchDay in range(idxDayWin.size):
                 currDayIdx = int(idxDayWin[loopSearchDay])
                 # get the maximum number of fragments to search on this day -
                 # that is across dimension two of nearbyPool['fragments'].
                 potFragCountDay = nGoodDays[todaysSeason][currDayIdx]
                 for loopFrag in range(int(potFragCountDay)):
                     potFragState = (
                      fragmentsState[todaysSeason][currDayIdx,loopFrag]
                     )
                     potFragDailyDepth = (
                        fragmentsDailyDepth[todaysSeason][currDayIdx,loopFrag]
                     )
                     if potFragState == todayWetState:
                         # it is of the same state so now check the depths
                         absDiffDepth = abs(todaysRain - potFragDailyDepth)

                         # Attempt to insert if it is smaller than something in
                         # our current potential fragments.  But only if the
                         # absolute differnce is less than param.absDiffTol
                         if absDiffDepth <= param['absDiffTol'] * todaysRain:
                             for loopInsert in range(potFrags[:,0].size):
                                 if absDiffDepth < potFrags[loopInsert,potFragsIdxAbsDiff]:
                                     # Increment our found counter
                                     fragFoundCounter += 1

                                     # We have a new good candidate, shift the
                                     # others down and insert
                                     if loopInsert == potFrags[:,0].size-1:
                                         potFrags[-1, potFragsIdxAbsDiff] = absDiffDepth
                                         potFrags[-1, potFragsIdxDayYear] = currDayIdx
                                         potFrags[-1, potFragsIdxDayGood] = loopFrag
                                     else:
                                         idxShiftOldStart = loopInsert
                                         idxShiftNewStart = loopInsert + 1 

                                         # Shift Down #Check index referencing
                                         potFrags[idxShiftNewStart:idxShiftNewEnd, potFragsIdxAbsDiff] = (
                                             potFrags[idxShiftOldStart:idxShiftOldEnd, potFragsIdxAbsDiff]
                                         )
                                         potFrags[idxShiftNewStart:idxShiftNewEnd, potFragsIdxDayYear] = (
                                             potFrags[idxShiftOldStart:idxShiftOldEnd, potFragsIdxDayYear]
                                         )
                                         potFrags[idxShiftNewStart:idxShiftNewEnd, potFragsIdxDayGood] = (
                                             potFrags[idxShiftOldStart:idxShiftOldEnd, potFragsIdxDayGood]
                                         )

                                         # Insert the new
                                         potFrags[loopInsert, potFragsIdxAbsDiff] = absDiffDepth
                                         potFrags[loopInsert, potFragsIdxDayYear] = currDayIdx
                                         potFrags[loopInsert, potFragsIdxDayGood] = loopFrag
                                    
                                     # We found a target, inserted it so break out
                                     break
            
            if fragFoundCounter == 0:
                # We did not find any potential fragments, so continue without
                # inserting anything.
                warningCounter = warningCounter + 1
                warningSum += todaysRain
                
                # If no potential fragments are found, find the frag with the closest fragmentDailyDepth
                # if the difference between the closest fragmentDailyDepth and todays rain is within 100%
                # use it. If not assume rainfall is spread evenly across day.
                minAbsDiffDepth = 1e9
                for loopSearchDay in range(idxDayWin.size):
                    currDayIdx = int(idxDayWin[loopSearchDay])
                 # get the maximum number of fragments to search on this day -
                 # that is across dimension two of nearbyPool['fragments'].
                    potFragCountDay = nGoodDays[todaysSeason][currDayIdx]
                    for loopFrag in range(int(potFragCountDay)):
                        potFragState = (
                        fragmentsState[todaysSeason][currDayIdx,loopFrag]
                        )
                        potFragDailyDepth = (
                            fragmentsDailyDepth[todaysSeason][currDayIdx,loopFrag]
                        )
                        if potFragState == todayWetState:
                            # it is of the same state so now check the depths
                            absDiffDepth = abs(todaysRain - potFragDailyDepth)
                            if absDiffDepth < minAbsDiffDepth:
                                minAbsDiffDepth = absDiffDepth
                                tmpIdxDay = currDayIdx 
                                tmpIdxFrg = loopFrag

                if minAbsDiffDepth < todaysRain:
                    tmpFragDepth = (
                        fragmentsDailyDepth[todaysSeason][tmpIdxDay,tmpIdxFrg]
                    )
                    selectedFrag = (
                        fragments[todaysSeason][tmpIdxDay,tmpIdxFrg]
                    )
                    # Scale:
                    selectedFrag = selectedFrag / tmpFragDepth * todaysRain
                    # Insert
                    subDailySim[:, idxDayLinear] = selectedFrag
                else:
                    subDailySim[:, idxDayLinear] = (todaysRain/recordsPerDay)
                
                idxDayLinear += 1
                continue
                

            # Woo HOO! We now have, hopefully, param.maxNearNeighb similar
            # fragments to sample from.  So, grab a uniform random number and
            # then select our fragment:
            #
            # First check if we found at least param.maxNearNeighb, if so no
            # futher action is required.
            currProbSampleVec = probSampleVecFull
            if fragFoundCounter <= param['maxNearNeighb']:
                currProbSampleVec[0:fragFoundCounter] = (
                    computeSubDailyKDE(fragFoundCounter)
                )
                # And force the maximum to avoid dropping through
                currProbSampleVec[fragFoundCounter-1] = 1.1
            
            randSelector = random.random()
            for loopProb in range(int(param['maxNearNeighb'])):
                if randSelector <= currProbSampleVec[loopProb]:
                    # We will use the fragment identified at this loopProb
                    break

            
            # Finally, we have a fragment so scale it and insert it
            tmpIdxDay = int(potFrags[loopProb, 1])
            tmpIdxFrg = int(potFrags[loopProb, 2])
            tmpFragDepth = (
                fragmentsDailyDepth[todaysSeason][tmpIdxDay,tmpIdxFrg]
            )
            selectedFrag = (
                fragments[todaysSeason][tmpIdxDay,tmpIdxFrg]
            )
            # With the conveniance variables above, scale:
            selectedFrag = selectedFrag / tmpFragDepth * todaysRain
            # Insert
            subDailySim[:, idxDayLinear] = selectedFrag
            idxDayLinear += 1    

    if warningCounter > 0:
        days_perc = warningCounter / nDays * 100
        days_perc = int(round(days_perc,0))
        avrf = warningSum / warningCounter
        avrf = int(round(avrf,0))

        print(
            'There were ' +
            str(warningCounter) +
            ' (' + str(days_perc) + '%)' +
            ' days where no potential fragment was found.' +
            ' The average rainfall depth on those days was ' +
            str(avrf) + ' mm'
        )
    return subDailySim