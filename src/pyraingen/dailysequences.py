# Packages & Libraries
import numpy as np
import netCDF4 as nc
from datetime import date
from numba.typed import List
#import nvtx

# Defined Functions
from .jdtodatevec import jdToDateVec
from .getdaysperyear import getDaysPerYear
from .getwetstate import getWetState

## Global constants
from .global_ import missingDay
from .global_ import stateBad
from .global_ import ndaysYearLeap
from .global_ import idxfebTwentyNine
from .global_ import yesterday, today, tomorrow

#@nvtx.annotate()
def dailySequences(nSeasons, 
                    nYearsPool,  
                    stnDetails, 
                    nearStationIdx, 
                    param, 
                    param_path):
    """Computes the daily sequences.
    Loop through the number of years in the seasonal
    pool and load the daily sequences.

    Parameters
    ----------
    nSeaons : int
        Number of seasons.
    nYearsPool : int
        Number of years in each seasonsal pool.
    stnDetails : dict where keys are words and values are float
        The record of stations from index.txt.  This is
        passed back as it may have changed shape to add according to
        genSeqOption3.
    nearStationIdx : array
        A set of indices into station data for the
        seasonal nearest like stations. This array is nominally two
        dimensional with the first dimension being the seasons and the
        second the number of stations returned.  This is because there
        are possible different correlations between stations over
        different seasons. 
    param : dict where keys are words and values are float
        Dictionary of the run parameters. 
    param_path : dict where keys are words and values are str
        Dictionary of the necessary paths.\n

    Returns
    ----------
    dailyDepth : list
        A list of arrays (one for each season) containing the daily depth sequences. 
    dailyWetState : list
        A list of arrays (one for each season) containing the daily wetState sequences.
    """                
    dailyDepth = List()
    dailyWetState = List()
    workingRainDepth = np.zeros((3,1))
    for loopSeason in range(nSeasons):
        # Allocate RAM for this daily series:
        dailyDepth.append(np.ones(
            (ndaysYearLeap,
            int(nYearsPool[loopSeason])))
            *missingDay)
        dailyWetState.append(np.ones(
            (ndaysYearLeap,
            int(nYearsPool[loopSeason])))
            *stateBad)
        
        #Counter through the year dimension of our arrays
        idxYear = 0
        
        for loopStation in range(nearStationIdx[0,:].size):
            # Grab a conveniance variable:
            currStnIndex = int(nearStationIdx[loopSeason, loopStation])

            if currStnIndex == 0:
            # There are no more stations for this season
                break
            else:
                # Get the start and end years from the NetCDF then compute the
                # number of years in the sequence.  We also want to reshape our
                # linear vector of rainfall data from the NetCDF into the
                # year-on-year array herein.
                fnameNC = ('{}/plv{:06}.nc'.format(param_path['pathSubDaily'], 
                    int(stnDetails['stnIndex'][currStnIndex]))
                )
                ds=nc.Dataset(fnameNC)
                daySeries = ds['day'][:].data
                dayVecStart = jdToDateVec(daySeries[0])
                dayVecEnd = jdToDateVec(daySeries[-1])
                yearStart = int(dayVecStart[0])
                yearEnd = int(dayVecEnd[0])
                tmpSubDaily = ds['rainfall'][:].data/10
                # The algorithm below works on the assumption that the
                # tmpSubDaily array is populated with full years.  So pad out
                # the data array to make full years with missingDay values.
                nDaysKnown = (date.toordinal(date(yearEnd,12,31))
                    - date.toordinal(date(yearStart,1,1))+1)
                tmpDaily = np.ones((1,nDaysKnown)) * missingDay
                # Get an index into the tmpDaily array where our known data
                # should start.
                dataIdxStart = (date.toordinal(date(
                    int(dayVecStart[0]), int(dayVecStart[1]), int(dayVecStart[2])))
                    - date.toordinal(date(yearStart,1,1)))
                dataIdxEnd = (date.toordinal(date(
                    int(dayVecEnd[0]),int(dayVecEnd[1]),int(dayVecEnd[2])))
                    - date.toordinal(date(yearStart,1,1))+1)

                # The sum here works because the method described in the paper
                # is based on increments per day.  We do a bit of a fudge when
                # dealing with daily read data aggregated to 9am.
                tmpDaily[0,dataIdxStart:dataIdxEnd] = tmpSubDaily.sum(axis=2)

                # This is the index into the linear array tmpDaily.
                idxDayLinear = 0
                # Loop over each year in this station
                for loopYear in range(yearStart,yearEnd+1):
                    # Loop over the days and compute the wetness index.  I
                    # don't think that this can be parallelised as it works on
                    # offset indexing.
                    nDaysCurrYear =  getDaysPerYear(loopYear)
                    #idxYear += 1

                    #for loopDay in range(nDaysCurrYear):
                    for loopDay in range(ndaysYearLeap):         
                        # We are now looping through the days of a particular
                        # year.
                        if loopDay == 0 and loopYear == yearStart:
                            # No previous day, assume dry as its a pretty good bet
                            # in Australia, except in the tropics. :-[, got to
                            # assume something.
                            #
                            # This is the first time through, so populate our
                            # working array
                            workingRainDepth[yesterday] = 0
                            workingRainDepth[today]     = tmpDaily[0, idxDayLinear]
                            workingRainDepth[tomorrow]  = tmpDaily[0, idxDayLinear + 1]
                        elif (loopDay == idxfebTwentyNine and 
                        nDaysCurrYear != ndaysYearLeap):
                            # If we are on the 29th day of the year but this is
                            # not a leap year increment loopDay by one as that
                            # is the index we are looping through our 2D array
                            # year with.  Or, just continue to the next
                            # iteration of the day loop.  But don't index our
                            # linear counter.
                            continue
                        elif loopDay == (ndaysYearLeap-1) and loopYear == yearEnd:
                            # This is the end of our sequence, assume dry as
                            # per the start
                            workingRainDepth[tomorrow] = 0
                        else:
                            # This is a normal day out in the middle of the
                            # sequence, so get ready for the next loop:
                            idxDayLinear = idxDayLinear + 1
                            workingRainDepth[yesterday] = workingRainDepth[today]
                            workingRainDepth[today]     = workingRainDepth[tomorrow]
                            workingRainDepth[tomorrow] = tmpDaily[0, idxDayLinear + 1]
                            
                    
                        # Store the depth for this day
                        dailyDepth[loopSeason][loopDay, idxYear] = (
                                tmpDaily[0, idxDayLinear]
                                )

                        # Decide our wetness state for this day
                        if np.all(workingRainDepth >= 0.0 - np.spacing(abs(workingRainDepth))):
                            # But only if we have "good" data.  The default is
                            # to leave it as zero
                            dailyWetState[loopSeason][loopDay,idxYear] = (
                                getWetState(workingRainDepth, param['dryWetCutoff'])
                            )
                    idxYear += 1
    return dailyDepth, dailyWetState                