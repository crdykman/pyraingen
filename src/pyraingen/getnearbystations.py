import numpy as np
import pandas as pd
import xarray as xr
import math

def station(param, target, nAttributes=33, fout='nearby_station_details.out'):
    ## Read the Station Data
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
            stnData['LON'][i] - target['lon'] < 0.001 and
            stnData['LAT'][i] - target['lat'] < 0.001
        )]
        canUseStation[Targetidx] = 0
    
    ## Reduce Our Station List
    # That is, work with only those that we can use.
    stnToUse = stnData.sel(index = canUseStation)

    ## Read the Coefficients Data
    modelCoeffs = np.transpose(pd.read_csv(param['pathModelCoeffs']).values)

    ## Loop over the Stations
    invPredictor = np.zeros((len(stnToUse['INDEX']), nAttributes))
    #sortIdxInvPredictor = np.zeros(np.shape(invPredictor))
    invPredMax = np.zeros((nAttributes,1))
    stnWeight = np.zeros((param['nNearStns'],1))

    for loopStation in np.arange(0,len(stnToUse['INDEX']),1):
        deltaLat = abs(target['lat'] - stnToUse['LAT'][loopStation])
        deltaLon = abs(target['lon'] - stnToUse['LON'][loopStation])
        if deltaLat < 0.001 and deltaLon < 0.001:
            print('warnings.warn(''additional dump'')')

        deltaDistToCoast = (abs(target['distToCoast'] 
            - stnToUse['DIST_COAST'][loopStation]) 
            / ((target['distToCoast'] 
            + stnToUse['DIST_COAST'][loopStation])/2)
            )
        deltaElevation = (abs(target['elevation']
            - stnToUse['ELEVATION'][loopStation])
            / ((target['elevation'] 
            + stnToUse['ELEVATION'][loopStation])/2)
            )
        deltaLatLon = deltaLat * deltaLon
        deltaTemp = (abs(target['temp']
            - stnToUse['av_an_tmax'][loopStation])
        )

        for loopAttr in range(nAttributes):
            # Extract a subset array to make the code simpler to read:
            modelCoeffsSubSet = modelCoeffs[loopAttr, :]

            invPredictor[loopStation, loopAttr] = (1.0 / (1.0 + math.exp(-1 * (
                modelCoeffsSubSet[0]
                + modelCoeffsSubSet[1] * deltaLat
                + modelCoeffsSubSet[2] * deltaLon
                + modelCoeffsSubSet[3] * deltaLatLon
                + modelCoeffsSubSet[4] * deltaDistToCoast
                + modelCoeffsSubSet[5] * deltaElevation 
                + modelCoeffsSubSet[6] * deltaTemp)))
                )
            if invPredictor[loopStation, loopAttr] > invPredMax[loopAttr]:
                invPredMax[loopAttr] = invPredictor[loopStation, loopAttr]
        
    ## Now, for this season compute the combined predictor values
    # This is based on the F77 implementation and is a multi-step operation:
    #   -) for each attribute
    #       -) normalise by the maximum
    #       -) sort descending
    #       -) store the rank INDEX, not the value
    #   -) add the

    pValue = np.zeros((len(stnToUse['INDEX']),1)) 

    for loopStation in range(len(stnToUse['INDEX'])):
            for loopAttr in range(nAttributes):
                pValue[loopStation] = (pValue[loopStation] + 
                    invPredictor[loopStation, loopAttr] / invPredMax[loopAttr]
                )

    pValue = pValue / nAttributes

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