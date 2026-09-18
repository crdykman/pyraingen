import os

import numpy as np
import pandas as pd
import xarray as xr
from importlib import resources

from .daily.core import run as _run_daily
from .getnearbystations import station
from .producedailynetcdf import produceDailyNetCDF


def regionaliseddailysim(nyears, startyear, nsim,
                            targetidx, targetlat, targetlon,
                            targetelev, targetdcoast, targetanrf,
                            targettemp, data_path,
                            pathStnData=None, pathModelCoeffs=None,
                            output_path_nc='daily.nc',
                            output_stats='stat_.out', output_val='rev_.out',
                            cutoff=0.30, wind=15, nstation=5,
                            nlon=3, lag=1, iamt=1, ival=0, irf=1, rm=1.0,
                            z = [[2,2,2],[90,180,345]],
                            getstations=True, seed=131):
    """Front end to the regionalised daily rainfall generator.
    For a detailed description of the approach refer to:
    Mehrotra et al (2012) "Continuous Rainfall simulation: 2.
    A regionalised daily rainfall generation approach".

    Parameters
    ----------
    nyears : int
        Number of years to be simulatied.
    startyear : int
        Start year of simulation.
    nsim : int
        Number of realisations to be simulated.
    targetidx : int
        Id number of target 'station'. 5 digits.
    targetlat : float
        Latitude of target.
    targetlon : float
        Longitude of target.
    targetelev: float
        Elevation of target in m (AHD).
    targetdcoast : float
        Distance to coast of target in km.
    targetanrf : float
        Annual rainfall at target in mm.
    targettemp : float
        Average annual maximum daily temperature at target in C.
    data_path : str
        Path to historical daily rainfall data.
    pathStnData : str
        Name and location for station record csv. Leave as default unless using own file.
        Default = None.
    pathModelCoeffs : str
        Name and location for daily logistic regression coefficients. Leave as default unless using own file.
        Default = None.
    output_path_nc : str
        Name and location for output netcdf.
        Default is 'daily.nc' in working directory.
    output_stats : str
        NB: No longer used. The Fortran implementation ignored this argument
        and always wrote its diagnostics to 'stat.out'; the pure-Python
        generator writes no diagnostics file. Accepted for backwards
        compatibility.
    output_val : str
        NB: Currently not supported. Leave as default.
    cutoff : float
        Threshold for defining wet and dry days. Days with rainfall above this threshold will be
        defined as wet and below as dry.
        Default is 0.30.
    wind : int
        half the size of window to search around target day.
        Default is 15 (i.e. search 15 days before and 15 days after).
    nstation : int
        Number of nearby stations to find.
        Default is 5.
    nlon : int
        Leave at default.
        Default is 3.
    lag : int
        Leave at default.
        Default is 1.
    iamt : int
        Leave at default.
        Default is 1.
    ival : int
        NB: Currently not supported. Leave as default.
        Default is 0.
    irf : int
        Leave at default.
        Default is 1.
    rm : float
        Rainfall multiplier to scale mean annual rainfall.
        Default is 1.0
    z : list
        Start and end of each long memory window, in days.
        Default is [[2,2,2],[90,180,345]].
    getstations : boolean
        Get nearby stations or use existing "nearby_station_details.out"
        Default is True
    seed : int
        Seed for the random number generator.
        Default is 131, the value the Fortran implementation hard-coded.

        NB: the generator reseeds with ``max(-seed, 1)``, so every value
        ``>= 0`` produces the same stream as the default. Pass a *negative*
        seed to obtain a different realisation set.

    Returns
    ----------
    netCDF
        Saves netCDF of daily simulations to specified file path.
    """

    nearby_path = 'nearby_station_details.out'

    if getstations:
        # Check if file exists
        if os.path.exists(nearby_path):
            os.remove(nearby_path)

        # Get Data
        if pathStnData == None:
            with resources.path("pyraingen.data", "stn_record.csv") as f:
                pathStnData = str(f)
        if pathModelCoeffs == None:
            with resources.path("pyraingen.data", "daily_logreg_coefs.csv") as f:
                pathModelCoeffs = str(f)

        # Get Nearby Stations
        param = {}
        param['nNearStns']       = nstation
        param['pathStnData']     = pathStnData
        param['pathModelCoeffs'] = pathModelCoeffs
        param['pathDailyData']   = data_path

        target = {}
        target['index']           = targetidx
        target['lat']             = targetlat
        target['lon']             = targetlon
        target['elevation']       = targetelev
        target['distToCoast']     = targetdcoast
        target['annualRainDepth'] = targetanrf
        target['temp']            = targettemp

        print("\nFinding nearby stations...\n")
        station(param, target, nAttributes=33, fout=nearby_path)

    # Begin simulation
    dailyData, dayVector = _run_daily(
        data_path=data_path,
        nearby_path=nearby_path,
        nyears=nyears,
        startyear=startyear,
        nsim=nsim,
        nstation=nstation,
        cutoff=cutoff,
        wind=wind,
        nlon=nlon,
        lag=lag,
        iamt=iamt,
        z=z,
        seed=seed,
    )

    if os.path.exists(output_path_nc):
        os.remove(output_path_nc)
    produceDailyNetCDF(output_path_nc, dailyData.astype(np.float64), dayVector)

    # Bias Correct and/or Scale Rainfall
    daily_rain = xr.open_dataset(output_path_nc)
    daily_rain['day'] = pd.to_datetime(daily_rain['day'], unit='D', origin='julian')
    # Mean annual total over all years and realisations. Equivalent to the
    # previous resample(day="A"), but without the frequency alias that pandas
    # renamed in 2.2.
    smanrf = daily_rain['rainfall'].groupby("day.year").sum().mean()
    daily_rain.close()
    del daily_rain
    daily_rain_bc = xr.open_dataset(output_path_nc)
    daily_rain_bc['rainfall'] *= (targetanrf * rm)/smanrf
    daily_rain_bc.close()
    daily_rain_bc.to_netcdf(output_path_nc)
