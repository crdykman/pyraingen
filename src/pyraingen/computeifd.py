import numpy as np

def computeIFD(rainfallSeries, yearsVector, ifdDurations, minsPerSample=6):
    """This function computes an IFD table from the input rainfall series
    based on a set of standard durations and recurrence intervals.  Please
    Note: this may not be as per a full IFD computation but is designed to
    be used inside the IFD Conditioning code.

    Parameters
    ----------
    rainfallSeries : array
        A matrix of the rain fall time series in six
        minute increments (240 per day), with dimenions (240, nDays,
        nSimulations).
    yearsVector : array  
        A vector of length nDays filled with integers of
        the years.  It is used to slice into the rainfallSeries.  For
        example, if the first year is 1942 then the first 365 records are
        1942, followed by 365 records of 1943, etc.\n

    Returns
    ----------
    array
        The resultant IFD table with dimensions (nYears,
        nSimulations, nIFDDurations)
    """
    #   Caleb Dykman
    #   2022-05-23

    from .aggregaterainfall import aggregateRainfall
    
    # Define some constants:
    # Standard duration of IFD used
    # NB: ifdDurations are in mintes for compatability with further calculations


    #Extract a conveniance vector of years
    yearsUnique = np.unique(yearsVector)

    #Allocate RAM for storage arrays:
    IFD = np.zeros((len(yearsUnique),
                    np.size(rainfallSeries, axis=0),
                    len(ifdDurations)))
    
    # The general process for each simulation is:
    #   -) Aggregate up from 6 minute if required.
    #   -) Extract the annual maximum series
    #   -) ?

    for loopDuration in np.arange(0,len(ifdDurations),1):
        # Aggregate if required
        if ifdDurations[loopDuration] == minsPerSample:
            # no aggregation required as we know that we are at the six minute
            # interval:
            rainfallAggregated = rainfallSeries
        else:
            # Well yes we do have to aggregate up
            # First: compute the number of records we will be aggregating
            # together for the required duration:
            nRecordsToAggregate = int(ifdDurations[loopDuration] / minsPerSample)

            # Then hand off to a dedicated function now we know the level that
            # we wish to aggregate to.
            rainfallAggregated = aggregateRainfall(
                rainfallSeries, nRecordsToAggregate
            )
        
        # Extract annual maximum series
        for loopYears in np.arange(0,np.size(IFD,axis=0),1):
            IFD[loopYears, :, loopDuration] = (
                np.max(np.max(rainfallAggregated[:,yearsVector == yearsUnique[loopYears],:], axis=2),axis=1)
            )
        
    # Finally, before we hand back, sort these
    IFD = np.sort(IFD, axis=0)

    return IFD



