import numpy as np

def ifdcondobjfun(targetIFD, simulatedIFD, Freq):
    """Calculates Objective function (relative mean absolute error (RMAE))
    between target and simulated IFD.

    Parameters
    ----------    
    targetIFD : array
        Target IFD table with dimensions (nYears, nSimulations, nIFDDurations).
    simulatedIFD : array
        Simulated IFD table with dimensions (nYears, nSimulations, nIFDDurations).
    Freq : array
        Vector of recurrence intervals.\n

    Returns
    ----------
    array
        Vector of relative mean absolute error for nSimulations.
    """
    # Obtain exceedence probabilites of simuated IFD corresponding (closer) to
    # the probabilites of target IFD

    # Standard Bom AEP
    #AEP = np.array([63.20, 50, 20, 10, 5, 2]).astype(int)

    # Frequency
    #Freq = 100 - np.append(1,np.append(np.arange(5,100,5), 99))


    # Compute the length of the annual series once and then reuse it:
    lengthTS = np.size(simulatedIFD, axis=0)

    # What we do first is go through and match the ranks 
    # (required for the exceedence probability Equation P = m/(n+1)) 
    # in our computed IFDs to the standard frequencies.
    # i.e. The Ranks that provides the exceedence probability closest to the standard frequency (smallest difference)

    indexFrequency = np.zeros((len(Freq), 1))
    for loopFreq in np.arange(0,len(Freq),1):
        indexFrequency[loopFreq] = np.argmin(
            abs(Freq[loopFreq] - 100 * np.arange(1,lengthTS+1,1) / (lengthTS + 1))
        )
    indexFrequency = indexFrequency.astype(int)
    
    # Calculate Objective function (relative mean absolute error (RMAE)).
    # This works because of two things:
    #   1) we have already sorted our simulated IFDs elsewhere, and;
    #   2) we have computed the index we need above.
    # We can, therefore, use the indices to slice through the IFD arrays and
    # compute the averages.  In sequence we:
    #   1) use indexFrequency to slice the years record into our standard
    #   frequencies.
    #   2) then take the mean across the second dimension - that is the number
    #   of simulations.
    #   3) the outer loop of durations loops through the third dimension of the
    #   simulatedIFD matrix.
    rmae = np.zeros((np.size(targetIFD, axis=1)))
    for loopDuration in np.arange(0,np.size(targetIFD, axis=1),1):
        # Index explicitly rather than squeezing: np.squeeze also collapses the
        # simulation axis when nSims == 1, which made the mean below raise.
        sims = simulatedIFD[indexFrequency[:, 0], :, loopDuration]
        meanSims = np.mean(sims, axis=1)
        rmae[loopDuration] = (
            np.mean(abs(meanSims - targetIFD[:, loopDuration]))
            / np.mean(targetIFD[:, loopDuration])
        )
    
    return rmae