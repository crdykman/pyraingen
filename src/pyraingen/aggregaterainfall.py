import numpy as np

def aggregateRainfall(rainfallSeries, nRecordsToAggregate):
    """This function is designed to aggregate, via a summation, six minute
    recorded rainfalls into longer durations.  The longer durations are
    specified as an integer number of records to aggregate and must divide
    evenly into 240 records, which is the number of six minute records per
    day.  That is mod(nRecordsToAggregate, 6) == 0.

    Parameters
    ----------
    rainfallSeries : array
        Array of the rain fall time series in six
        minute increments (240 per day), with dimenions (240, nDays,
        nSimulations).
    nRecordsToAggregate : int
        Records to aggregate.\n

    Returns
    ----------
    array
        The aggregated rainfall series.
    """
    #   Caleb Dykman
    #   2022-05-23

    # Define some utility constants:
    nRecordsToAggregate = int(nRecordsToAggregate)
    # Derive the records per day from the data rather than assuming six-minute
    # sampling, so non-default sampling intervals aggregate correctly.
    nRecordsPerDay = np.size(rainfallSeries, axis=2)
    if nRecordsPerDay % nRecordsToAggregate != 0:
        raise ValueError(
            f"{nRecordsToAggregate} records per aggregation does not divide "
            f"evenly into {nRecordsPerDay} records per day; the remainder "
            "would be silently discarded"
        )

    rainfallAggregated = np.zeros((np.size(rainfallSeries, axis=0),
                                    np.size(rainfallSeries, axis=1),
                                    int(nRecordsPerDay/nRecordsToAggregate)))
    
    # Now, loop through the day and sum up
    loopCounter = 0
    for loopAggregate in np.arange(0,nRecordsPerDay, nRecordsToAggregate):
        rainfallAggregated[:, :, loopCounter] = (
            np.sum(rainfallSeries[:,:,loopAggregate:loopAggregate+nRecordsToAggregate], axis=2)
        )
        loopCounter = loopCounter + 1

    return rainfallAggregated

