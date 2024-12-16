import numpy as np
import math

def targetStations(param, param_path, genSeqOption3, stnDetails):
    """Computes the nearest stations
    This seeks to implement the nearby station routines of
    Westra et al 2011, section 3.1 and 3.3.

    Parameters
    ----------
    param : dict where keys are words and values are float
        Dictionary of the run parameters. 
    param_path : dict where keys are words and values are str
        Dictionary of the necessary paths. 
    genSeqOption3 : dict where keys are words and values are float
        Dictioary of the parameters in case the user
        chooses genSeqOption==3.
    stnDetails : dict where keys are words and values are float
        The station details array structure.\n

    Returns
    ----------
    stnDetails : dict where keys are words and values are float
        The record of stations from index.txt.  This is
        passed back as it may have changed shape to add according to
        genSeqOption3.
    nearStationIdx : array
        A set of indices into station data for the
        seasonal nearest like stations. This array is nominally two
        dimensional with the first dimension being the seasons and the
        second the number of stations returned.  This is because there
        are possible different correlations between stations over
        different seasons.
    """
    #   Dr Peter Brady <peter.brady@wmawater.com.au>
    #   2016-08-10=
    #
    #   Caleb Dykman
    #   2021-09-15

    ## Import the Globals
    from .global_ import nSeasons, nAttributes
    from .readcoefficients import readCoefficients

    ## Read the station details
    # Check if our target station, param.targetIndex, is within our data, if so
    # store the index.  Otherwise force it to the number of stations plus 1 and
    # read additional data from our genSeqOption3 data to become our
    # reference.
    if param['targetIndex'] in stnDetails['stnIndex']:
        idxTarget = np.where(stnDetails['stnIndex'] == param['targetIndex'])
    else:
        idxTarget = np.size(stnDetails['stnIndex'], axis=0)
        #assuming all arrays are of same length so indicies match
        stnDetails['stnIndex']     = np.append(stnDetails['stnIndex'], idxTarget)
        stnDetails['stnLat']       = np.append(stnDetails['stnLat'], genSeqOption3['lat'])
        stnDetails['stnLon']       = np.append(stnDetails['stnLon'], genSeqOption3['lon'])
        stnDetails['stnElevation'] = np.append(stnDetails['stnElevation'], genSeqOption3['elevation'])
        stnDetails['stnDistCoast'] = np.append(stnDetails['stnDistCoast'], genSeqOption3['distToCoast'])
        stnDetails['stnTemp']      = np.append(stnDetails['stnTemp'], genSeqOption3['temperature'])

    ## Allocate Ram
    nearStationIdx = np.zeros((nSeasons,(np.size(stnDetails['stnIndex'], axis=0))))

    ## Read the Coefficients Data
    # Check directory for importing function
    modelCoeffs = readCoefficients(param_path['pathCoeff'])

    for loopSeason in range(nSeasons):
        # Allocate and zero arrays
        invPredictor = np.zeros(((np.size(stnDetails['stnIndex'], axis=0)),nAttributes))
        invPredMax = np.zeros((nAttributes,1))

        for loopStation in range((np.size(stnDetails['stnIndex'], axis=0))):
            # For this station compute the model in equations 9 and 10,
            # although I think the algorithm actually runs the reverse of the
            # description in the paper.  That is, this is the inverse of
            # Equation 9.  See:
            #
            #   https://en.wikipedia.org/wiki/Logit (retrieved 2016-08-10)
            #
            # for a more detailed discussion of the inverse logistic
            # relationships.
            #
            # First up extract some conveniance variables.  These computations
            # are roughly equivalent to step (2) in Section 4.1 of the paper.
            deltaLat = (np.abs(stnDetails['stnLat'][idxTarget] 
                - stnDetails['stnLat'][loopStation])
                )
            deltaLon = (np.abs(stnDetails['stnLon'][idxTarget] 
                - stnDetails['stnLon'][loopStation])
                )
            deltaDistToCoast = (np.abs(stnDetails['stnDistCoast'][idxTarget] 
                - stnDetails['stnDistCoast'][loopStation]) 
                / ((np.abs(stnDetails['stnDistCoast'][idxTarget] 
                + stnDetails['stnDistCoast'][loopStation]))/2)
                )
            deltaElevation = (np.abs(stnDetails['stnElevation'][idxTarget]
                - stnDetails['stnElevation'][loopStation])
                / ((np.abs(stnDetails['stnElevation'][idxTarget] 
                + stnDetails['stnElevation'][loopStation]))/2)                                                                 
                )
            deltaLatLon = deltaLat * deltaLon
            deltaTemp = (np.abs(stnDetails['stnTemp'][idxTarget]
                - stnDetails['stnTemp'][loopStation])
                )

            # Now loop over the combination of attributes for this station and
            # season:
            for loopAttr in range(nAttributes):
                #Extract a subset array to make the code simpler to read:
                modelCoeffsSubSet = np.squeeze(modelCoeffs[loopSeason,:,loopAttr])
                
                # Now we compute the inverse of equation 9.  However,
                # exponetiation is computationally expensive so only do it if
                # we have to, that is for loopStation != idxTarget.
                if loopStation != idxTarget:
                    invPredictor[loopStation, loopAttr] = (
                        1.0 / (1.0 + math.exp(-1 * (
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
                else:
                    invPredictor[loopStation, loopAttr] = 0.0

        # Now, for this season compute the combined predictor values
        pValue = np.zeros(((np.size(stnDetails['stnIndex'], axis=0)),1)) #too many brackets?
        for loopStation in range((np.size(stnDetails['stnIndex'], axis=0))):
            for loopAttr in range(nAttributes):
                pValue[loopStation] = (pValue[loopStation] +
                    invPredictor[loopStation,loopAttr] / invPredMax[loopAttr]
                    )

        # Now sort in descending and get the rank index
        pValueSortIdx = np.argsort(pValue, axis=0)[::-1]

        # With these ranked stations loop through until we get at least the
        # number of years in the record asked for
        yearSum = 0.0
        stationCounter = 0

        for loopProc in np.arange(0,(np.size(stnDetails['stnIndex'], axis=0)),1):
            # Extract some conveniance variables for readability:
            currIPer = stnDetails['stnIPer'][pValueSortIdx[loopProc][0]]
            currDataL = stnDetails['stnDataL'][pValueSortIdx[loopProc][0]]

            # Compute the record length for this station:
            currLen = currDataL * float(currIPer) / 100.0

            # Check that this station contains the minimum length or record.
            # If so store it and continue.
            if currLen > param['minYears']:
                # Store this station:
                nearStationIdx[loopSeason, stationCounter] = (
                    pValueSortIdx[loopProc][0]
                    )

                # Increment our counters
                stationCounter += 1
                yearSum += currLen
                
                # Do we have enough data for this season?
                if yearSum > param['nYearsPool']:
                    # Yes! WOO HOO! So break out.
                    break
    
    return stnDetails, nearStationIdx