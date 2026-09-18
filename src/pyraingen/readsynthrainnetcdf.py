import netCDF4 as nc

def readSynthRainNetCDF(fname):
    """Reads formated sythetic rain series from NetCDF
    This function reads a sequence of synthetic rainfall series from a
    NetCDF formatted, compressed file.

    Parameters
    ----------
    fname : str
        The name and location of the file to read from.\n

    Returns
    ----------
    yearStart : int 
        The year that the record starts in.
    yearEnd : int 
        The year that the record ends in.
    data : array 
        The data in the shape of (240, nDays, nRecordsToCondition)
    """
    #   Dr Peter Brady <peter.brady@wmawater.com.au>
    #   2015-07-15
    #
    #   Caleb Dykman
    #   2022-05-23

    #Defined Functions
    from .jdtodatevec import jdToDateVec

    # Read the headers
    ds = nc.Dataset(fname)
    
    daySeries = ds['day'][:].data
    dayVecStart = jdToDateVec(daySeries[0])
    dayVecEnd = jdToDateVec(daySeries[-1])
    yearStart = int(dayVecStart[0])
    yearEnd = int(dayVecEnd[0])

    #Preallocate RAM
    data = ds['rainfall'][:].data
    ds.close()

    return yearStart, yearEnd, data