from pyraingen.regionalisedsubdailysim import regionalisedsubdailysim

fnameInput = "daily.nc"
pathSubDaily = "C:/Users/z3460382/OneDrive - UNSW/PhD/Continuous_Rainfall_Simulation/Caleb/Python/data/sub_daily_netcdf/"
targetIndex = 66037
targetlat = -33.9410 
targetlon = 151.1730 
targetelev = 6.00 
targetdcoast = 0.27 
targetanrf = 1086.71
targettemp = 22.40
fnameSubDaily = "subdaily.nc"

regionalisedsubdailysim(
    fnameInput, pathSubDaily, targetIndex, 
    pathIndex=None, pathCoeff=None, pathReference=' ', 
    fnameSubDaily=fnameSubDaily,
    minYears=10, nYearsPool=500, dryWetCutoff=0.30,
    halfWinLen=15, maxNearNeighb=10, nSims=10,
    genSeqOption=3, nYearsRef=50,
    absDiffTol=0.1, gso3_lat=targetlat, gso3_lon=targetlon,
    gso3_elev=targetelev, gso3_distcoast=targetdcoast, 
    gso3_anrf=targetanrf, gso3_temp=targettemp
)