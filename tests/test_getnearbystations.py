from pyraingen.getnearbystations import station

# Get Nearby Stations
param = {}
param['nNearStns']       = 5
param['pathStnData']     = '../src/pyraingen/data/stn_record.csv'
param['pathModelCoeffs'] = '../src/pyraingen/data/daily_logreg_coefs.csv'
param['pathDailyData']   = "C:/Users/z3460382/OneDrive - UNSW/PhD/Continuous_Rainfall_Simulation/Caleb/Python/reference_implementation/daily/daily_data/"

target = {}
target['index']           = 66037      
target['lat']             = -33.9410      
target['lon']             = 151.1730    
target['elevation']       = 6.00     
target['distToCoast']     = 0.27      
target['annualRainDepth'] = 1086.71      
target['temp']            = 26.40    

station(param, target, nAttributes=33,
    fout=f'nearby_station_details.out'
)