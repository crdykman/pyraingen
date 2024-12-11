# Packages & Libraries
import numpy as np
import netCDF4 as nc
from datetime import date
import matplotlib.pyplot as plt
from numba.typed import List
#import nvtx

# Defined Functions
from .jdtodatevec import jdToDateVec
from .getdaysperyear import getDaysPerYear

## Global constants
from .global_ import missingDay
from .global_ import stateBad
from .global_ import ndaysYearLeap
from .global_ import idxfebTwentyNine
from .global_ import recordsPerDay

#@nvtx.annotate()
def getFragments(nSeasons, nGoodDays, dailyWetState, dailyDepth, stnDetails, nearStationIdx, param, param_path):
    """Loops over the stations andload and store only the fragments whose
    wetstate != 0 (i.e. some possibly good data).

    Parameters
    ----------
    nSeaons : int
        Number of seasons.
    nGoodDays : array
        Number of good days per year per day per season.
        A "good day" is one that is not of state bad and has
        a depth greater than the dryWetCutoff.
    dailyDepth : list
        A list of arrays (one for each season) containing the daily depth sequences. 
    dailyWetState : list
        A list of arrays (one for each season) containing the daily wetState sequences.
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
    fragments : list
        A list of arrays (one for each season) containing a count  
        of the number of good fragments available for sampling.
    fragmentsState : list
        A list of arrays (one for each season) containing the daily wetState sequences
        of the good fragments available for sampling.
    fragmentsDailyDepth : list
        A list of arrays (one for each season) containing the daily depth sequences
        of the good fragments available for sampling.
    """
    
    # Allocate RAM
    fragments = List()
    fragmentsState = List()
    fragmentsDailyDepth = List()

    for loopSeason in range(nSeasons):
        # Append Season Array to list
        maxNGoodDays = int(np.max(nGoodDays[loopSeason]))
        fragments.append(
            np.zeros((ndaysYearLeap, maxNGoodDays, recordsPerDay)))
        fragmentsState.append(
            np.zeros((ndaysYearLeap, maxNGoodDays)))
        fragmentsDailyDepth.append(
            np.zeros((ndaysYearLeap, maxNGoodDays)))
        fragmentCounter = np.zeros((ndaysYearLeap, 1))
        #Counter through the year dimension of our arrays
        idxYear = 0
        for loopStation in range(nearStationIdx[0,:].size):
            # Grab a convenience variable:
            currStnIndex = int(nearStationIdx[loopSeason, loopStation])
            if currStnIndex == 0:
                # There are no more stations for this season
                break
            else:
                fnameNC = ('{}/plv{:06}.nc'.format(param_path['pathSubDaily'], 
                    int(stnDetails['stnIndex'][currStnIndex])))
                ds=nc.Dataset(fnameNC)
                daySeries = ds['day'][:].data
                dayVecStart = jdToDateVec(daySeries[0])
                dayVecEnd = jdToDateVec(daySeries[-1])
                yearStart = int(dayVecStart[0])
                yearEnd = int(dayVecEnd[0])
                # As per above, we need to pad out to Jan 1 and Dec 31.
                nDaysKnown = (date.toordinal(date(yearEnd,12,31))
                    - date.toordinal(date(yearStart,1,1))+1)
                # Get an index into the SubDaily array where our known data
                # should start.
                dataIdxStart = (date.toordinal(date(
                    int(dayVecStart[0]),int(dayVecStart[1]),int(dayVecStart[2])))
                    - date.toordinal(date(yearStart,1,1)))
                dataIdxEnd = (date.toordinal(date(
                    int(dayVecEnd[0]),int(dayVecEnd[1]),int(dayVecEnd[2])))
                    - date.toordinal(date(yearStart,1,1))+1)

                tmpSubDaily = np.ones((nDaysKnown, recordsPerDay,)) * missingDay #dimensions flipped
                tmpSubDaily[dataIdxStart:dataIdxEnd, :] = ds['rainfall'][:].data/10 #
                
                # This is the index into the days dimension of tmpSubDaily
                idxDayLinear = 0
                # Loop over each year in this station
                for loopYear in range(yearStart,yearEnd+1):
                    nDaysCurrYear = int(getDaysPerYear(loopYear))
                    #for loopDay in range(nDaysCurrYear):
                    for loopDay in range(ndaysYearLeap):
                        # still need to handle a leap year...
                        if (loopDay == idxfebTwentyNine and 
                        nDaysCurrYear != ndaysYearLeap):    
                            continue
                        # But otherwise work through:
                        if (dailyWetState[loopSeason][loopDay,idxYear] != stateBad and 
                        dailyDepth[loopSeason][loopDay,idxYear] > param['dryWetCutoff']):
                            
                            fragments[loopSeason][loopDay,int(fragmentCounter[loopDay]),:] = \
                                tmpSubDaily[idxDayLinear,:]
                            
                            fragmentsState[loopSeason][loopDay,int(fragmentCounter[loopDay])]= \
                                dailyWetState[loopSeason][loopDay, idxYear]
                            
                            fragmentsDailyDepth[loopSeason][loopDay,int(fragmentCounter[loopDay])]= \
                                dailyDepth[loopSeason][loopDay,idxYear]

                            # Sanity Check:
                            if (abs(np.sum((
                                fragments[loopSeason][loopDay,int(fragmentCounter[loopDay]),:]))
                                - dailyDepth[loopSeason][loopDay,idxYear]) > 1):
                                
                                plt.plot((
                                fragments[loopSeason][loopDay,int(fragmentCounter[loopDay]),:]
                                ))
                                raise ValueError('Sum fail')
                            fragmentCounter[loopDay] += 1
                        idxDayLinear += 1
                    idxYear += 1
    
    return fragments, fragmentsState, fragmentsDailyDepth