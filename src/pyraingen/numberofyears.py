import numpy as np
import netCDF4 as nc
#import nvtx

# Defined Functions
from .jdtodatevec import jdToDateVec

#@nvtx.annotate()
def numberOfYears(nSeasons, stnDetails, nearStationIdx, param_path):
    """Computes the number of years in each seasonal pool.

    Parameters
    ----------
    nSeaons : int
        Number of seasons.
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
    param_path : dict where keys are words and values are str
        Dictionary of the necessary paths.\n

    Returns
    ----------
    array
        Number of years in each seasonsal pool.

    Notes
    -----
    This counts the calendar years *spanned* by each station's record
    (``endYear - startYear + 1``), not the years of usable data. Gaps within a
    record and the completeness fraction ``stnIPer`` are both ignored, so a
    station with a long but sparse record contributes its full span. The pool
    is therefore an upper bound, and ``param['nYearsPool']`` is reached with
    fewer stations than the available data alone would justify.
    """
    nYearsPool = np.zeros((nSeasons,1))
    for loopSeason in range(nSeasons):    
        # initialise our variables for this season:
        for loopStation in range(nearStationIdx[0,:].size):
        # Grab a conveniance variable:
            currStnIndex = int(nearStationIdx[loopSeason, loopStation])
            if currStnIndex < 0:
                # There are no more stations for this season
                break
            else: 
                # Get the start and end years from the NetCDF then compute the
                # number of years in the sequence.
                fnameNC = ('{}/plv{:06}.nc'.format(param_path['pathSubDaily'], 
                    int(stnDetails['stnIndex'][currStnIndex])))
                ds=nc.Dataset(fnameNC)
                daySeries = ds['day'][:].data
                ds.close()
                dayVecStart = jdToDateVec(daySeries[0])
                dayVecEnd = jdToDateVec(daySeries[-1])
                nYearsPool[loopSeason] += (dayVecEnd[0] 
                    - dayVecStart[0] + 1) 
    return nYearsPool