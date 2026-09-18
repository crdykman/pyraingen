import numpy as np
import pandas as pd
import warnings
from importlib import resources

def station(param, target, nAttributes=33, fout='nearby_station_details.out'):
    """Algorithm for finding nearby daily stations. 

    Parameters
    ----------
    param : dict
        parameter dictionary requiring 'nNearStns', 
        'pathStnData', 'pathModelCoeffs', 'pathDailyData'.
    target : dict 
        target data dictionary requiring 'index', 
        'lat', 'lon', 'elevation', 'distToCoast', and
        'annualRainDepth'.
    nAttributes : int
        Number of attributes of similarity. Leave as default.
        default = 33.
    fout : str
        Path to save output. Leave as default.
        Default is 'nearby_station_details.out'.

    Returns
    ----------
    str
        Prints nearby daily stations.
        Saves text file of nearby daily stations
        to specified file path.
    """    
    ## Read the Station Data
    # Get Data
    if param['pathStnData'] == None:
        with resources.path("pyraingen.data", "stn_record.csv") as f:
            param['pathStnData'] = str(f)
    stnData = pd.read_csv(param['pathStnData']).to_xarray()

    # Allocate RAM
    canUseStation = np.ones((len(stnData.index)), dtype=bool)

    # Check if our target station is within the data set.  If it is then its
    # correlation will be perfect and we want to exclude it.
    if target['index']  in stnData['INDEX']:
        Targetidx = stnData['INDEX'].where(
            stnData['INDEX'] == target['index'], drop=True).squeeze().index
        canUseStation[Targetidx] = 0

        # As per the original F77 source, also exclude any stations such that:
        #   abs(deltaLon) < 0.001 AND abs(deltaLat) < 0.001
        Targetidx = [i for i in range(len(stnData['INDEX'])) if (
            abs(stnData['LON'][i] - target['lon']) < 0.001 and
            abs(stnData['LAT'][i] - target['lat']) < 0.001
        )]
        canUseStation[Targetidx] = 0
    
    ## Reduce Our Station List
    # That is, work with only those that we can use.
    stnToUse = stnData.sel(index = canUseStation)

    ## Read the Coefficients Data
    if param['pathModelCoeffs'] == None:
        with resources.path("pyraingen.data", "daily_logreg_coefs.csv") as f:
            param['pathModelCoeffs'] = str(f)
    modelCoeffs = np.transpose(pd.read_csv(param['pathModelCoeffs']).values)

    ## Compute the Predictors For Every Station At Once
    # This was a scalar double loop over stations x attributes, indexing into
    # xarray DataArrays one element at a time, which cost roughly 200 s for the
    # 1694 stations in stn_record.csv. The arithmetic below is identical; only
    # the loop is gone.
    stnWeight = np.zeros((param['nNearStns'],1))

    lat   = stnToUse['LAT'].values
    lon   = stnToUse['LON'].values
    coast = stnToUse['DIST_COAST'].values
    elev  = stnToUse['ELEVATION'].values
    tmax  = stnToUse['av_an_tmax'].values

    deltaLat = np.abs(target['lat'] - lat)
    deltaLon = np.abs(target['lon'] - lon)
    if np.any((deltaLat < 0.001) & (deltaLon < 0.001)):
        warnings.warn(
            'one or more candidate stations are co-located with the target'
        )

    deltaDistToCoast = (np.abs(target['distToCoast'] - coast)
                        / ((target['distToCoast'] + coast) / 2))
    deltaElevation = (np.abs(target['elevation'] - elev)
                      / ((target['elevation'] + elev) / 2))
    deltaTemp = np.abs(target['temp'] - tmax)

    # Columns match the coefficient order: intercept, lat, lon, lat*lon,
    # distance to coast, elevation, temperature.
    design = np.column_stack([
        np.ones_like(deltaLat),
        deltaLat,
        deltaLon,
        deltaLat * deltaLon,
        deltaDistToCoast,
        deltaElevation,
        deltaTemp,
    ])
    # einsum rather than the @ operator so this does not depend on a working
    # BLAS; at this size the two are equivalent in speed.
    invPredictor = 1.0 / (1.0 + np.exp(
        -np.einsum('ij,kj->ik', design, modelCoeffs[:nAttributes, :7])
    ))

    # The running maximum in the original loop is just the column maximum.
    invPredMax = invPredictor.max(axis=0)

    ## Now, for this season compute the combined predictor values
    # This is based on the F77 implementation and is a multi-step operation:
    #   -) for each attribute
    #       -) normalise by the maximum
    #       -) sort descending
    #       -) store the rank INDEX, not the value
    #   -) add the
    pValue = (invPredictor / invPredMax).sum(axis=1)[:, None] / nAttributes

    # Now sort and get the rank index
    pValueSort = np.sort(pValue, axis=0)[::-1]
    PvalueSortIdx = np.argsort(pValue, axis=0)[::-1]

    ## Compute the Station Weights
    # But only on the subset;
    weightSum = 0
    for loopStation in range(param['nNearStns']):
        stnWeight[loopStation] = pValueSort[loopStation]
        weightSum += stnWeight[loopStation]

    stnWeight = stnWeight / weightSum

    nearbyStn = {
        'stnIndex': [],
        'weight' : [] ,                      
        'nYears' : [] ,              
        'avAnRain' : [] ,
        'startYear' : [] ,
    } 

    ## As they are now sorted, grab how many we need
    for loopStore in range(param['nNearStns']):
        nearbyStn['stnIndex'].append(stnToUse['INDEX'][PvalueSortIdx[loopStore]].data[0])
        nearbyStn['nYears'].append(stnToUse['NYEAR'][PvalueSortIdx[loopStore]].data[0])
        nearbyStn['avAnRain'].append(stnToUse['AN_RAINFALL'][PvalueSortIdx[loopStore]].data[0])
        
        # The following has a different index as we already use the sort list
        # above when computing the weights.
        nearbyStn['weight'].append(stnWeight[loopStore][0])
    
    ## Read Start year from file
    for i in nearbyStn['stnIndex']:
        with open(param['pathDailyData'] + f'rev_dr{i:06d}.txt') as f:
            nstart = int(f.readline()[43:47])
        nearbyStn['startYear'].append(nstart)

    ## Write nearby station details file
    with open(fout, 'w') as f:
        f.write('    No Index Weight Years St_year Av annual rainfall\n')
        f.write('\n')
        f.write(' target Station\n')
        f.write(f"     0 {target['index']}  1.000     0    -1   {target['annualRainDepth']}\n")
        f.write('\n')
        f.write(' Nearby Stations\n')
        for loopwrite in range(param['nNearStns']):
            f.write('%6d%6d%7.3f%6d%6d%10.2f\n' % 
                (loopwrite+1,
                nearbyStn['stnIndex'][loopwrite],
                nearbyStn['weight'][loopwrite],
                nearbyStn['nYears'][loopwrite],
                nearbyStn['startYear'][loopwrite],
                nearbyStn['avAnRain'][loopwrite])
        )
        
    ## Write nearby station details file
    print('    No Index Weight Years St_year Av annual rainfall')
    print(' target Station')
    print(f"     0 {target['index']}  1.000     0    -1   {target['annualRainDepth']}")
    print(' Nearby Stations')
    for loopprint in range(param['nNearStns']):
        print('%6d%6d%7.3f%6d%6d%10.2f' % 
            (loopprint+1,
            nearbyStn['stnIndex'][loopprint],
            nearbyStn['weight'][loopprint],
            nearbyStn['nYears'][loopprint],
            nearbyStn['startYear'][loopprint],
            nearbyStn['avAnRain'][loopprint])
    )
    print()