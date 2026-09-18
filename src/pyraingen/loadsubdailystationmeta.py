import numpy as np
import netCDF4 as nc

def loadSubDailyStationMeta(fname):
    """Loads sub-daily station meta data
    Specifically this function loads the reference data for the stations
    used in the logistic regression to find the "nearest" similar stations.

    Parameters
    ----------
    fname : str 
        The name and location of the NetCDF to open the data from.\n

    Returns
    ----------
    dict where keys are words and values are arrays
        A dictionary of all stations holding the data fields:\n
        -) Index\n
        -) Longitude\n
        -) Latitude\n
        -) Start month\n
        -) Start year\n
        -) End month\n
        -) End year\n
        -) Data length\n
        -) Iper\n
        -) Elevation\n
        -) Distance to coast\n
        -) Average annual maximum daily temperature
    """
    #   Dr Peter Brady <peter.brady@wmawater.com.au>
    #   2016-08-09

    #   Caleb Dykman
    #   2021-09-16
    
    ds = nc.Dataset(fname)

    stnIndex      = ds['index'][:].data
    stnLon        = ds['lon'][:].data
    stnLat        = ds['lat'][:].data
    stnStartMonth = ds['startmonth'][:].data
    stnStartYear  = ds['startyear'][:].data
    stnEndMonth   = ds['endmonth'][:].data
    stnEndYear    = ds['endyear'][:].data
    stnDataL      = ds['datal'][:].data
    stnIPer       = ds['iper'][:].data
    stnElevation  = ds['elevation'][:].data
    stnDistCoast  = ds['disttocoast'][:].data
    stnTemp       = ds['avg_annual_max_daily_temp'][:].data
    ds.close()

    # Convert into a structure (matlab) / dictionary (python)
    
    stnData = {}

    stnData['stnIndex'] =      stnIndex
    stnData['stnLon'] =        stnLon
    stnData['stnLat'] =        stnLat
    stnData['stnStartMonth'] = stnStartMonth
    stnData['stnStartYear'] =  stnStartYear
    stnData['stnEndMonth'] =   stnEndMonth
    stnData['stnEndYear'] =    stnEndYear
    stnData['stnDataL'] =      stnDataL
    stnData['stnIPer'] =       stnIPer
    stnData['stnElevation'] =  stnElevation
    stnData['stnDistCoast'] =  stnDistCoast
    stnData['stnTemp'] =       stnTemp

    return stnData