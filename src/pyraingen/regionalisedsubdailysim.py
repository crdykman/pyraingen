## Regionalised Sub-Daily Code
# This code is part of the ARR Project 4 port from the research quality
# code written by UNSW into a production ready code, that is its the
# development component of R&D.  Really, this is more of a ground up
# rewrite than a port as we can take significant advantage of modern
# compilers dynamic memory allocation, data storage and data structures
# with native compression.
#
# Specifically, this is a front end to the regionalised sub-daily
# disaggregation code.  As this is a prototype a number of values that
# would normally be entered by the user at run time are hard coded during
# development.  This limitation would be removed prior to release.
#
# In this code when I refer to "the paper" I am refering to Westra et al
# 2012 "Continuous rainfall simulation 1: A regionalised subdaily
# disaggregation approach"
#
# Dr Peter Brady <peter.brady@wmawater.com.au>
# 2016-08-09

# Caleb Dykman
# 2021-09-17

# Reset Interpreter
#from IPython import get_ipython
#get_ipython().magic('reset -sf')

# Close all figure winows?
#cls = lambda: print("\033c\033[3J", end='')
#cls()

# Packages & Libraries
import numpy as np
import random
from datetime import date
from numba.core import types
from numba.typed import Dict

# Defined Functions
from .loadsubdailystationmeta import loadSubDailyStationMeta
from .targetstations import targetStations
from .getndailysims import getNDailySims
from .datevectojd import datevecToJD
from .producesubdailynetcdf import produceSubDailyNetCDF

## Global constants
from .global_ import nSeasons
from .global_ import ndaysYearLeap

def regionalisedsubdailysim(fnameInput, pathSubDaily, pathIndex, 
                            pathCoeff, targetIndex, pathReference=' ', 
                            fnameSubDaily='subdaily.nc',
                            minYears=10, nYearsPool=500, dryWetCutoff=0.30,
                            halfWinLen=15, maxNearNeighb=10, nSims=10,
                            genSeqOption=3, nYearsRef=50,
                            absDiffTol=0.1, gso3_lat=None, gso3_lon=None,
                            gso3_elev=None, gso3_distcoast=None, 
                            gso3_anrf=None, gso3_temp=None):
    """Front end to the regionalised sub-daily
    disaggregation code.

    Options (genSeqOption) included:\n
    0. Only sub daily data is available. 
        - Daily data is formed from the subdaily record and is used for disaggregation.\n
    1. Daily and sub-daily rainfall data at target station is available \n(daily record comes from other source and therefore may be of different length to sub-daily record).
        - using sub-daily record at target station to disaggregate the daily record at the target station.\n
    2. Only daily rainfall at target station is available
        - using nearby locations sub-daily information disaggregate daily rainfall at target location.\n
    3. No daily or sub-daily rainfall information is available at target station.
        - first derive the daily rainfall series at the target location using daily rainfall simulator.
        - and subsequently perform daily rainfall disaggreation using sub-daily information from nearby locations.
        - under this option many realisations of daily rainfall are available.\n
    4. Multiple simulated sequences of daily rainfall at the target station are available. 
        - sub daily rainfall sequences are generated using sub-daily information at the target station.\n
    
    Step 1: Pre-Compute the Pool - Mix of Step 1 and 2 in Westra et al
    2012 "Continuous rainfall simulation 1: A regionalised subdaily
    disaggregation approach". This depends on the user's input of 
    genSeqOption. Different sequences are generated for each season, because if 
    genSeqOption = 3 then there may be different stations per
    season.  However, if the target data is known (genSeqOption <= 2)
    then it is the same data set for each season.  Finally, if if there is 
    simulated daily data then it is used exclusively for each season with one daily
    simulation per sub-daily simulation. If there are more sub-daily
    simulations than there are daily then daily simulations are recycled.

    First is to compute the nearby station data. The key to note here is that 
    the method of fragments is only concerned with one parameter that can be computed up
    front:
    
      1) the day of the year the rain fell, e.g. November 23rd, and not which
    year and the where is already handled with the "nearby" stations".
    
    The second logical must be computed sliding through the target daily
    series:

      2) a logical if the wet states for the previous and next day are the
    same as from the target daily sequence.
    
    The general process is then to loop over the seasons and then each station
    within that season. From there we compute the sequences.
    
    Step a) Compute the number of years in the seasonal pool.

    Step b) Compute the daily sequences in the seasonal pool.

    Step c) Work out the maximum number of good days per year per day per season.

    Step d) Load and store only the fragments whose wetstate != 0 
    (i.e. some possibly good data)

    Step 2 a) If Required Load Daily Reference Data.

    Step 2 b) Dissagregation Loop

    Save Data

    Parameters
    ----------
    fnameInput : str
        This is the path to the input file name.  It can be either a single
        sub-daily data set, a single daily data or a block of daily simulations
        depending on the value of genSeqOption.
    pathSubDaily : str
        Path to the sub-daily pluviograph records
    pathIndex : str
        Path to station details data. Use get_index() function.
    pathCoeff : str
        Path to the logistic regression coefficients for selection of nearby stations.
        Use get_coeffs() function.
    targetIndex : int
        Index of the site to perform sub-daily disaggregation. Can be an existing
        station in the station details dataset or a chosen 6 digit index linked to 
        the gso3 data inputs.
    pathReference : str
        Path to the observed single location daily data, for genSeqOption<=2.
        Default is ' '.
    fnameSubDaily : str
        File name to dump the resultant disaggregated data into.
        Default is 'subdaily.nc'.
    minYears : float
        Minimum number of years data to consider. 
        Default is 10.
    nYearsPool : float
        Total number of years to pool from adjacent stations. 
        Default is 500.
    dryWetCutoff : float 
        Threshold for determing whether day is wet or dry. 
        Default is 0.30.
    halfWinLen : float
        Half the length of the search window in days for finding fragments from nearby days. 
        Default is 15.
    maxNearNeighb : float
        The maximum number of nearest neighbour fragments to sample from. 
        Default is 10.
    nSims : float
        Number of simulations to be performed.
        Default is 10.
    genSeqOption : int 
        genSeqOption = 0: no daily read is required, convert a known
        subdaily sequence.
        genSeqOption = 1: there is observed daily data at this
        location, so read it in and compute the requisite series from that
        data set.
        genSeqOption = 2: as per param.genSeqOption = 1.
        genSeqOption = 3: simulated daily data is available at this
        station with multiple simualtion, e.g. from the Daily section of
        this project, repeat the process for each daily simulation as a
        basis.
        genSeqOption = 4: as per param.genSeqOption = 3.
        Default is 3.
    nYearsRef : float
        number of years of daily/sub-daily data available. 
        Default is 50.
    absDiffTol : float
        The difference in the target/source rainfall, decimal percent.
        Default is 0.1.
    gso3_lat : float
        Latitude of simulation site for genSeqOption = 3, if not in station
        details dataset.
        Default is None
    gso3_lon : float
        Longitude of simulation site for genSeqOption = 3, if not in station
        details dataset.
        Defaults is None.
    gso3_elev : float
        Elevation of simulation site for genSeqOption = 3, if not in station
        details dataset.
        Default is None.
    gso3_distcoast : float
        Distance to coast of simulation site for genSeqOption = 3, if not in station
        details dataset.
        Default is None. 
    gso3_anrf : float
        Average annual rainfall of simulation site for genSeqOption = 3, if not in station
        details dataset.
        Default is None.
    gso3_temp : float
        Average annual maximum daily temperature of simulation site for genSeqOption = 3, if not in station
        details dataset.
        Default is None.\n

    Returns
    ----------
    netCDF
        Saves netCDF of disaggregated sub-daily simulations
        to specified file path.
    """
    print('Initialising...')
    # Set the random seed:
    # Python uses the Mersenne Twister as the core generator.
    # current system time is used
    random.seed()

    ## Hard Coded Constants Normally Read from data_inf.dat
    # From the Fortran source the variable ipos is replaced with genSeqOption:
    # Options (genSeqOption) included:
    # 0. Only sub daily data is available. Daily data is formed from the
    #   subdaily record and is used for disaggregation
    # 1. Daily and sub-daily rainfall data at target station is available
    #   (daily record comes from other source and therefore may be of different
    #   length to sub-daily record
    #       - using sub-daily record at target station to disaggregate the
    #           daily record at the target station
    # 2. Only daily rainfall at target station is available
    #   - using nearby locations sub-daily information disaggregate daily
    #       rainfall at target location
    # 3. No daily or sub-daily rainfall information is available at target
    #   station
    #       - first derive the daily rainfall series at the target location
    #           using daily rainfall simulator
    #       - and subsequently perform daily rainfall disaggreation using
    #           sub-daily information from nearby locations
    #       - under this option many realisations of daily rainfall are
    #           available
    # 4. Multiple simulated sequences of daily rainfall at the target station
    #   are available - sub daily rainfall sequences are generated using
    #   sub-daily information at the target station
    

    # param = {}
    param = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.float64,
    )
    #Minimum number of years data to consider
    param['minYears']        = minYears
    #total number of years to pool from adjacent stations
    param['nYearsPool']      = nYearsPool 
    param['dryWetCutoff']    = dryWetCutoff
    param['halfWinLen']      = halfWinLen
    param['maxNearNeighb']   = maxNearNeighb
    param['nSims']           = nSims
    param['targetIndex']     = targetIndex
    #See below:
    param['genSeqOption']    = genSeqOption 
    #number of years of daily/sub-daily data available
    param['nYearsRef']       = nYearsRef
    #the difference in the target/source rainfall, decimal percent
    param['absDiffTol']      = absDiffTol 

    # This is the path to the input file name,  It can be either a single
    # sub-daily data set, a single daily data or a block of daily simulations
    # depending on the value of param.genSeqOption.
    # genSeqOption = 0
    # param['fnameInput']      = '../data/sub_daily_netcdf/plv066037.nc'
    # genSeqOption = 1, 2
    # param['fnameInput']      = '../data/daily_netcdf/rev_dr066037.nc'
    # genSeqOption = 3/4

    param_path = Dict.empty(
        key_type=types.unicode_type,
        value_type=types.unicode_type,
    )
    # param_path = {}
    param_path['fnameInput']      = fnameInput

    #File name to dump the resultant disaggregated data into.
    param_path['fnameSubDaily']   = fnameSubDaily

    # Paths of interest:
    # Path to the sub-daily pluviograph records:
    param_path['pathSubDaily'] = pathSubDaily
    # Path to the observed single location daily data, for genSeqOption<=2
    param_path['pathReference'] = pathReference
    # Path to index data
    param_path['pathIndex'] = pathIndex
    # Path to Coefficients
    param_path['pathCoeff'] = pathCoeff

    # Additional parameters for param.genSeqOption == 3:
    # These are read in the subroutine target_station if required

    genSeqOption3 = {}
    genSeqOption3['lat']             = gso3_lat
    genSeqOption3['lon']             = gso3_lon
    genSeqOption3['elevation']       = gso3_elev
    genSeqOption3['distToCoast']     = gso3_distcoast
    genSeqOption3['annualRainDepth'] = gso3_anrf
    genSeqOption3['temperature']     = gso3_temp

    ## Preparation
    stnDetails = loadSubDailyStationMeta(param_path['pathIndex'])
    if param['genSeqOption'] == 2 or param['genSeqOption'] == 3:
        stnDetails, nearStationIdx = targetStations(param, param_path,
            genSeqOption3, 
            stnDetails
            )
    else:
        # If we get to here we have some form of data and am not calling
        # targetStations to find the nearest stations.  Therefore, we set up
        # some dummy arrays for later use that hold only the target station of
        # interest.
        if param['targetIndex'] not in stnDetails['stnIndex']:
            raise ValueError('Target station not found in reference stations')
        else:
            idxTarget = stnDetails['stnIndex'].index(param['targetIndex'])
            nearStationIdx = np.ones((nSeasons,1)) * idxTarget
            del idxTarget

    ## Dissagregation
    # In the Fortran we'd jump into the dissagregation subroutine but its huge
    # so I'm going to try to split it up a bit.  Actually, I'm going to ignore
    # Sanjeev's code from here on in an rewrite based on the discussion in the
    # paper with reference to the code.

    ## Step 1: Pre-Compute the Pool - Mix of Step 1 and 2 in the Paper
    # This depends on the user's input of param.genSeqOption, specifically:
    #
    #   -) param.genSeqOption = 0: no daily read is required, convert a known
    #       subdaily sequence.
    #   -) param.genSeqOption = 1: there is observed daily data at this
    #       location, so read it in and compute the requisite series from that
    #       data set.
    #   -) param.genSeqOption = 2: as per param.genSeqOption = 1.
    #   -) param.genSeqOption = 3: simulated daily data is available at this
    #       station with multiple simualtion, e.g. from the Daily section of
    #       this project, repeat the process for each daily simulation as a
    #       basis.
    #   -) param.genSeqOption = 4: as per param.genSeqOption = 3.
    #
    #  ******************************
    #  I really should enumerate the above to remove constants in the code.
    #  ******************************
    #
    # We generate different sequences for each season because if we are working
    # with param.genSeqOption = 3 then we may have different stations per
    # season.  However, if we have known target data (param.genSeqOption <= 2)
    # we use the same data set for each season.  Finally, if we have simulated
    # daily data we use it exclusively for each season with one daily
    # simulation per sub-daily simulation.  If we have asked for more sub-daily
    # simulations than we have daily then we recyle the daily simulations.
    #
    # Irrespective, the nearStationsIdx will have the required stations by now
    # so run through it.
    #
    # First is to compute the nearby station data as we may need to use that
    # for our target location.  The key to note here is that the method of
    # fragments is only concerned with one parameter that can be computed up
    # front:
    #
    #   1) the day of the year the rain fell, e.g. November 23rd, and not which
    #       year and we have already handled where with the "nearby" stations".
    #
    # The second logical must be computed as we slide through the target daily
    # series:
    #   2) a logical if the wet states for the previous and next day are the
    #       same as from the target daily sequence.
    #
    # Therefore, we can save a lot of RAM by first computing the daily series
    # and then looking at only the rainfall days to store more information
    # where rain actually fell.

    # The general process is to loop over the seasons and then each station
    # within that season.  From there we compute the sequences.
    #
    # Step a) loop over the stations and compute the number of years in the
    # pool:

    print('Step 1(a) looping over stations and computing number of years')
    from .numberofyears import numberOfYears
    nYearsPool = numberOfYears(nSeasons, stnDetails, nearStationIdx, param_path)

    # Step b) compute the daily sequences.
    # Now as we know the number of years in the seasonal pool we can allocate
    # RAM and loop through the NetCDFs loading the daily sequences.
    #
    # workingRain depth holds three days of daily rainfall with the indices:
    #   (0 indexed for Python)
    #   0 == previous day
    #   1 == current day
    #   2 == next day
    # See the symbolic constants

    print('Step 1(b) Compute daily sequences')
    from .dailysequences import dailySequences
    dailyDepth, dailyWetState = dailySequences(nSeasons, nYearsPool, stnDetails, 
                                            nearStationIdx, param, param_path)

    # Step c) work out the maximum number of good days per year per day per season.
    # For this computation a "good day" is one that is not of state bad and has
    # a depth greater than param['dryWetCutoff'].

    print('Step 1(c) Compute maximum number of good days per year day per season')
    from .maxgooddays import maxGoodDays
    nGoodDays = maxGoodDays(nSeasons, nYearsPool, dailyWetState, dailyDepth, param)

    # Step d) with our daily sequences load and store only the fragments whose
    # wetstate != 0 (i.e. some possibly good data).  This does involve looping
    # over the stations again.

    print('Step 1(d) load and store only possibly good fragments')
    from .getfragments import getFragments
    fragments, fragmentsState, fragmentsDailyDepth = getFragments(nSeasons, 
                                                        nGoodDays, dailyWetState, 
                                                        dailyDepth,stnDetails, 
                                                        nearStationIdx, param, 
                                                        param_path)
    del dailyDepth, dailyWetState

    # Step 2 a) If Required Load Daily Reference Data
    # This section of the code runs differently to the UNSW version because I
    # computed the pool trying to reduce the RAM footprint.  Therefore, for
    # param.genSeqOption < 3, I have to reread and compute the daily sequence
    # that the original program already had computed.
    # So for param.genSeqOption >= 3 we have multiple daily simulations so
    # read each one individually.  The OR loopSim == 1 is to handle where
    # we have a known target, with daily, subdaily or a mix, and we only
    # load it once, compute the statistics and use therein.
    print('Step 2 a) Load Daily Reference Data')
    from .readdata import readData
    from .paddata import padData
    (ds,
    daySeries, 
    dayVecStart,
    dayVecEnd,
    simYearStart, 
    simYearEnd) = readData(param_path['fnameInput'])

    if param['genSeqOption'] == 0:
        # this is reading a known sub-daily that has no daily, so read and
        # aggregate.
        tmpSubDailyRain = ds['rainfall'][0,:,:].data
        # pad out the data array to make full years with missingDay values.
        targetDailyRain, dataIdxStart, dataIdxEnd = \
                    padData(simYearStart, simYearEnd, dayVecStart, dayVecEnd)
        targetDailyRain[dataIdxStart:dataIdxEnd] = tmpSubDailyRain.sum(axis=1)
        
        DayStart = dataIdxStart
        DayEnd = (date.toordinal(date(
            int(dayVecEnd[0]),int(dayVecEnd[1]),int(dayVecEnd[2])))
            - date.toordinal(date(simYearEnd,1,1)))

    elif param['genSeqOption'] == 1 or param['genSeqOption'] == 2:
        # param.genSeqOption == 1: there are known, daily and sub-daily,
        # I've already computed the sub-daily fragments above so read the
        # known daily sequence.
        # param.genSeqOption == 2: There are known daily, so actually this
        # reads the same file as param.genSeqOption == 1.
        tmpDaily = ds['rainfall'][0,:].data
        # pad out the data array to make full years with missingDay values.
        targetDailyRain, dataIdxStart, dataIdxEnd = \
                    padData(simYearStart, simYearEnd, dayVecStart, dayVecEnd)
        targetDailyRain[dataIdxStart:dataIdxEnd] = tmpDaily
        
        DayStart = dataIdxStart
        DayEnd = (date.toordinal(date(
            int(dayVecEnd[0]),int(dayVecEnd[1]),int(dayVecEnd[2])))
            - date.toordinal(date(simYearEnd,1,1)))

    else:
        # this is reading multiple daily simulations, read them per
        # simulation.  It actually covers param.genSeqOption == 3 and
        # param.genSeqOption == 4, which has some form of pre-simulated
        # daily rainfall data.
        
        # Check the number of simulations.  This is to handle recycling if
        # the user has asked for 100 simulations, say, but the source only
        # has 50.
        nDailySims = getNDailySims(param_path['fnameInput'])
        targetDailyRain = ds['rainfall'][:].data 
        DayStart = 0
        DayEnd = ndaysYearLeap-1

    param['simYearStart'] = simYearStart
    param['simYearEnd'] = simYearEnd
    param['DayStart'] = DayStart
    param['DayEnd'] = DayEnd

    ## Step 2 b) Dissagregation Loop
    print('Step 2 b) Performing subdaily disaggregation')
    # from .subdailysimloop0 import subDailySimLoop0
    # from .subdailysimloop1 import subDailySimLoop1
    from joblib import Parallel, delayed
    from .subdailydisaggregation import subDailyDisaggregation


    if param['genSeqOption'] >=0 and param['genSeqOption'] < 3:
        if __name__ == "__main__":
            results = Parallel(n_jobs=-1, prefer="threads")(delayed(subDailyDisaggregation)(
                                                    targetDailyRain, 
                                                    param,nGoodDays, 
                                                    fragments, 
                                                    fragmentsState, 
                                                    fragmentsDailyDepth) for i in range(int(param['nSims']))) 
    elif param['genSeqOption'] >=3 and param['genSeqOption'] <= 4:
        if param['nSims'] > nDailySims:
            # This is to handle recycling if the user has asked for 100 simulations, 
            # say, but the source only has 50.
            diff = param['nSims'] - nDailySims
            nSims = [*range(int(nDailySims))] + random.choices(range(int(nDailySims)), k=diff)
        else:
            nSims = range(int(param['nSims']))
        
        if __name__ == "__main__":
            results = Parallel(n_jobs=-1, prefer="threads")(delayed(subDailyDisaggregation)(
                                                    targetDailyRain[:, int(i)].astype(np.float64), 
                                                    param, nGoodDays, 
                                                    fragments, 
                                                    fragmentsState, 
                                                    fragmentsDailyDepth) for i in nSims)
    else:
        print("Unsupported genSeqOption value. Must be >= 0 and <= 4")

    if param['nSims'] > 1:
        subdailySims = np.dstack(results)
    else:
        subdailySims = results

    print('Saving data')
    jdStart = datevecToJD(date(int(simYearStart), 1, 1))
    simJDSeries  = np.arange(0,np.size(subdailySims, axis=1),1)+ jdStart
    produceSubDailyNetCDF(param_path['fnameSubDaily'],subdailySims ,simJDSeries)
