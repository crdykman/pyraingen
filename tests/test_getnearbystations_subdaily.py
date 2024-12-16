from pyraingen.nearbystations import nearbystations

pathIndex = None
pathCoeff = None
targetIndex = 66037
targetlat = -33.9410 
targetlon = 151.1730 
targetelev = 6.00 
targetdcoast = 0.27 
targetanrf = 1086.71
targettemp = 22.40


subdailystns = nearbystations(pathIndex, pathCoeff, targetIndex, 
                                lat=targetlat, lon=targetlon, elev=targetelev, 
                                distcoast=targetdcoast, anrf=targetanrf, temp=targettemp) 

print(subdailystns)