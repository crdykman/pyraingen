from importlib import resources

from .loadsubdailystationmeta import loadSubDailyStationMeta
from .targetstations import targetStations

def nearbystations(pathIndex, pathCoeff, targetIndex, 
                    lat=None, lon=None, elev=None, 
                    distcoast=None, anrf=None, temp=None, 
                    minYears=10, nYearsPool=500):
    """Finds nearby subdaily stations
    This seeks to implement the nearby station routines of
    Westra et al 2011, section 3.1 and 3.3.
    
    Parameters
    ----------
    pathIndex : str
        Path to station details data. Use get_index() function.
    pathCoeff : str
        Path to the logistic regression coefficients for selection of nearby stations.
        Use get_coeffs() function.
    targetIndex : int
        Index of the site to perform sub-daily disaggregation. Can be an existing
        station in the station details dataset or a chosen 6 digit index linked to 
        the gso3 data inputs.
    minYears : float
        Minimum number of years data to consider. 
        Default is 10.
    nYearsPool : float
        Total number of years to pool from adjacent stations. 
        Default is 500.
    lat : float
        Latitude of target site, if not in station
        details dataset.
        Default is None
    lon : float
        Longitude of target site, if not in station
        details dataset.
        Defaults is None.
    elev : float
        Elevation of target site, if not in station
        details dataset.
        Default is None.
    distcoast : float
        Distance to coast of target site, if not in station
        details dataset.
        Default is None. 
    anrf : float
        Average annual rainfall of target site, if not in station
        details dataset.
        Default is None.
    temp : float
        Average annual maximum daily temperature of target site, if not in station
        details dataset.
        Default is None.\n

    Returns
    ----------
    subdailystns : list
        A list of the set of subdaily station Ids for all the
        seasonal nearest like stations.
    """
    param = {}
    #Minimum number of years data to consider
    param['minYears']        = minYears #keep
    #total number of years to pool from adjacent stations
    param['nYearsPool']      = nYearsPool #keep
    param['targetIndex']     = targetIndex #keep

    # This is the path to the input file name,  It can be either a single
    # sub-daily data set, a single daily data or a block of daily simulations
    # depending on the value of param.genSeqOption.
    # genSeqOption = 0
    # param['fnameInput']      = '../data/sub_daily_netcdf/plv066037.nc'
    # genSeqOption = 1, 2
    # param['fnameInput']      = '../data/daily_netcdf/rev_dr066037.nc'
    # genSeqOption = 3/4

    param_path = {}
    # Get Data
    if pathIndex == None:
        with resources.path("pyraingen.data", "index.nc") as f:
            pathIndex = str(f)
    if pathCoeff == None:
        with resources.path("pyraingen.data", "coefficients.dat") as f:
            pathCoeff = str(f)
    # Path to index data
    param_path['pathIndex'] = pathIndex
    # Path to Coefficients
    param_path['pathCoeff'] = pathCoeff

    # Additional parameters for param.genSeqOption == 3:
    # These are read in the subroutine target_station if required

    genSeqOption3 = {}
    genSeqOption3['lat']             = lat
    genSeqOption3['lon']             = lon
    genSeqOption3['elevation']       = elev
    genSeqOption3['distToCoast']     = distcoast
    genSeqOption3['annualRainDepth'] = anrf
    genSeqOption3['temperature']     = temp

    ## Preparation
    stnDetails = loadSubDailyStationMeta(param_path['pathIndex'])
    stnDetails, nearStationIdx = targetStations(param, param_path,
                                                genSeqOption3, 
                                                stnDetails)

    s1 = ([stnDetails['stnIndex'][id] for id in nearStationIdx[0, :].astype(int) if id > 0])
    s2 = ([stnDetails['stnIndex'][id] for id in nearStationIdx[1, :].astype(int) if id > 0])
    s3 = ([stnDetails['stnIndex'][id] for id in nearStationIdx[2, :].astype(int) if id > 0])
    s4 = ([stnDetails['stnIndex'][id] for id in nearStationIdx[3, :].astype(int) if id > 0])
    
    subdailystns = list(set(s1 + s2 + s3 + s4))

    return subdailystns