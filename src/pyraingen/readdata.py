# Packages & Libraries
import netCDF4 as nc
#import nvtx

# Defined Functions
from .jdtodatevec import jdToDateVec

#@nvtx.annotate()
def readData(fnameInput):
    """Reads in daily rainfall data from netCDF file.
    
    Parameters
    ----------
    fnameInput : str
        Location and file name of input netCDF.\n

    Returns
    ----------
    ds : netCDF dataset
        Daily rainfall dataset including daily rainfall depths and julian dates.
    daySeries : array
        Array of julian date series. 
    dayVecStart : array
        Vector of year, month, day for start of simulation.
    dayVecEnd :  array
        Vector of year, month, day for end of simulation.
    simYearStart : int
        Start year for simulation.
    simYearEnd : int
        End year for simulation.
    """
    ds=nc.Dataset(fnameInput)
    daySeries = ds['day'][:].data
    # NB: ds is returned open on purpose -- the caller reads 'rainfall' from it
    # and is responsible for closing it.
    dayVecStart = jdToDateVec(daySeries[0])
    dayVecEnd = jdToDateVec(daySeries[-1])
    simYearStart = int(dayVecStart[0])
    simYearEnd = int(dayVecEnd[0])

    return ds, daySeries, dayVecStart, dayVecEnd, simYearStart, simYearEnd