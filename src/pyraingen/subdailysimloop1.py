from numba import njit
from numba.typed import List
import numpy as np
import random

from .subdailydisaggregation import subDailyDisaggregation

@njit
def subDailySimLoop1(nSims, nDailySims, targetDailyRain, param,
                    nGoodDays, fragments, fragmentsState, fragmentsDailyDepth):
    """Loops through nSims performing sub-daily disaggregation.
    for genSeqOption >=3 and <= 4.

    Parameters
    ----------
    nSims : float
        Number of simulations to be performed.
    targetDailyRain : array
        Array of daily rainfall series for target station.
    param : dict where keys are words and values are float
        Dictionary of the run parameters.
    nGoodDays : array
        Number of good days per year per day per season.
        A "good day" is one that is not of state bad and has
        a depth greater than the dryWetCutoff. 
    fragments : list
        A list of arrays (one for each season) containing a count  
        of the number of good fragments available for sampling.
    fragmentsState : list 
        A list of arrays (one for each season) containing the daily wetState sequences
        of the good fragments available for sampling.
    fragmentsDailyDepth : list 
        A list of arrays (one for each season) containing the daily depth sequences
        of the good fragments available for sampling.\n

    Returns
    ----------
    list
        List of arrays (one for each simulation) of size (recordsPerDay, nDays).
    """
    # Check the number of simulations.  This is to handle recycling if
    # the user has asked for 100 simulations, say, but the source only
    # has 50.
    subDailySims = List()
    for loopSim in range(nSims):
        if loopSim+1 > nDailySims:
            simCounter = random.randrange(nDailySims)
        else:
            simCounter = loopSim
        # Read the daily simulation
        targetDailyRainfall = targetDailyRain[:, int(simCounter)].astype(np.float64)
        subDailySim = subDailyDisaggregation(targetDailyRainfall, param,
                nGoodDays, fragments, fragmentsState, fragmentsDailyDepth)
        subDailySims.append(subDailySim)
    return subDailySims
