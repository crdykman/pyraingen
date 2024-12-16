## IFDCond
# Port of Fitsum's R code to MATLAB so that I can understand the method.
# It will most likely be further ported then to Fortran 2008 for deployment
# after the fact.
#
# Dr Peter Brady <peter.brady@wmawater.com.au>
# 2015-06-11
#
# This MATLAB port is based on a sample R script provided by Fitsum
# Woldemeskel <f.woldemeskel@unsw.edu.au>.

# Caleb Dykman
# 2022-05-24

# clear all %#ok<CLSCR>
# close all

# tic


# addpath('../../utils', '-end')

## Packages & Libraries
import numpy as np
from datetime import date
import random
import scipy.optimize
import warnings
import matplotlib.pyplot as plt
warnings.simplefilter('ignore', np.RankWarning)

## Defined Functions
from .readsynthrainnetcdf import readSynthRainNetCDF
from .datevectojd import datevecToJD
from .computeifd import computeIFD
from .ifdcondobjfun import ifdcondobjfun
from .aggregaterainfall import aggregateRainfall
from .ifdcorrection import correction
from .producesubdailynetcdf import produceSubDailyNetCDF

## Global constants
#from .global_ import missingDay
#from .global_ import recordsPerDay


def ifdcond(fileNameInput, fileNameOutput, fileNameTargetIFD, 
            nSims, nRecursions = 2,
            TargetIFDdurationsEst = [30, 60, 360, 720],
            TargetIFDdurations = [30, 60, 360, 720], 
            AEP = [63.20, 50, 20, 10, 5, 2],
            nRecordsPerDay = 240, minsPerSample = 6,
            massScale = 1, plot=True):
    """Algorithm from Fitsum et al (2016) for Constraining continuous 
    rainfall simulations for derived design flood estimation.

    Parameters
    ----------
    fileNameInput : str
        File name and location for raw input data.
    fileNameOutput : str
        File name and location to store corrected output data.
    fileNameTargetIFD : str 
        File name and location for target IFD data. Data should be csv with size
        rows = durations, cols = Frequencies (AEP). No headers or row labels.
    nSims : int
        Number of data simulations in input data that are to  be corrected. 
    nRecursions : int
        Number of recursions (repetitions) to be used in correction procedure
        Default is 2.
    TargetIFDdurationsEst : list
        IFD durations present in input target IFD data. 
        Defailt is [30, 60, 360, 720].
    TargetIFDdurations : list
        IFD durations in input target IFD data to be corrected.
        (i.e. subset of TargetIFDdurationsEst)
        Defailt is [30, 60, 360, 720]. 
    AEP : list
        List of return periods (frequency)  in input target IFD data to be corrected.
        Default is [63.20, 50, 20, 10, 5, 2].
    nRecordsPerDay : int
        Number of subdaily records per day. (Mins per day / minsPerSample)  
        Default is 240.
    minsPerSample : int
        Base timestep of input data in mins.
        Default is 6.
    massScale : float
        Scaling factor for adjusting rainfall mass (~mean annual rainfall)
        Default is 1.
    plot : bool
        True to output Raw vs Conditioned comparison plots for each duration.
        False to suppress.
        Default is True

    Returns
    ----------
    netCDF
        Saves netCDF of corrected sub-daily simulations
        to specified file path.
    """

    ## Now Read The Data
    yearStart, yearEnd, dataRaw = readSynthRainNetCDF(fileNameInput)
    dataRaw = np.swapaxes(dataRaw,0,2)   # for unravelling 

    # Create the time vector
    #dateStart = date(yearStart,1,1)
    #dateEnd =  date(yearEnd,12,31)
    nDays = (date(yearEnd,12,31) - date(yearStart,1,1)).days +1
    #rainfall = np.ones((recordsPerDay, nDays)) * missingDay
    dailyData = np.zeros((nDays,1))
    dayVector = np.arange(0,np.size(dailyData, axis=0))
    dayVector = dayVector + datevecToJD(date(int(yearStart),1,1))

    if nSims < 100:
        SimsIndex = random.sample(range(0, np.size(dataRaw, axis=0)), nSims)
        dataRaw = dataRaw[SimsIndex,:,:]


    # Target IFD
    # The columns of this data field match the IFD durations, which are stored
    # in TargetIFDdurationEst.  The rows match the Annual Exceedence
    # Probability stored in Freq.
    # NB: they are sorted in an alternate order though.  For example the
    # correct plot would be: plot(Freq, targetIFD(end:-1:1, 1))
    targetIFD = np.loadtxt(fileNameTargetIFD,delimiter=',')

    ## General Constants and Preparation
    # Are we dubugging:
    debugMode = False

    # Standard duration of IFD used
    # NB: these are in minutes for compatability with further calculations
    TargetIFDdurationEst = np.array(TargetIFDdurationsEst)

    # Duration of TargetTFD considered for bias correction, change these
    # numbers as required!!!
    TargetIFDduration = np.array(TargetIFDdurations)
    #TargetIFDduration = TargetIFDdurationEst[[0, 2, 4]]

    # Frequency
    # Can only constrain to recurrence intervals that would be present in the record!
    Freq = np.array(AEP).astype(int)

    # Vector of years we are simulating over:
    years = np.arange(yearStart, yearEnd+1, 1)

    #Check
    if Freq.min() < 1/(len(years)+1)*100:
        warnings.warn('''Can only constrain to recurrence intervals 
        that would be present in the length of the record. Check no. 
        of years and defined frequencies.''')

    # Allocate some space to store the objective function results:
    objectiveResults = np.zeros((
        nRecursions,
        len(TargetIFDduration),
        len(TargetIFDdurationEst)
    ))

    # Compute the data vector in MATLAB epoch for the known time series start
    # and end years:
    dateVector = np.arange(date(yearStart,1,1), date(yearEnd+1,1,1), 1)

    # Years Vector
    # The following vector is a conveniance to use when accumulating data up to
    # yearly scale.  It is simply the year component only of the date vector
    yearsVector = dateVector.astype('datetime64[Y]').astype(int) + 1970

    # Number of days in the data set:
    NoDays = len(dateVector)

    ## Reference IFD Computation
    # This is the first time through and this IFD for the raw data
    # set is not computed.  So do it
    IFDRaw = computeIFD(dataRaw, yearsVector, TargetIFDduration)

    # Also compute the reference objective function value from the
    # first time through.
    objectiveReference = ifdcondobjfun(targetIFD, IFDRaw, Freq)

    ## Allocate Additional RAM For Reference Variables
    # These are really only useful for debugging though so they will be
    # disabled by default

    if debugMode == True:
        tmp = np.reshape(dataRaw, (nSims, NoDays * nRecordsPerDay), 
                            order='C' )

        # The number of rainfall records in the series
        nRainDaysRaw = np.sum(tmp.astype(bool), axis=1)
        nRainDaysCor = np.zeros((nRecursions, len(TargetIFDduration), nSims))

        # Total RainFall Depth
        rainfallDepthRaw = np.sum(tmp, axis=1)
        rainfallDepthCor = np.zeros((nRecursions, len(TargetIFDduration), nSims))

        # Mean annual maximum
        meanAnnualMaxRaw = np.mean(max(np.reshape(dataRaw, (nSims, nRecordsPerDay*NoDays), 
                                        order='C')), axis=1)
        meanAnnualMaxCor = np.zeros((nRecursions, len(TargetIFDduration)))

        # Clean up:
        del tmp

    ## Main Computation Body
    # repeat correction for all the recursion (iteration) considered
    for currentRec in range(nRecursions):
        # repeat correction for all the duration of interest
        for currentDur in range(len(TargetIFDduration)):   
            print(
                f'Recursion {currentRec+1} of {nRecursions}, '
                f'Duration {currentDur+1} of {len(TargetIFDduration)}'
            )

            print(f'Step 1.{currentRec+1}.{currentDur+1}: Compute the Annual Maximum and Non-Maximum Series')
            ###################################################################
            # Step 1: Compute the Annual Maximum and Non-Maximum Series
            # Now we have to collapse the data set to the requested duration.
            # We do this in a few stages:
            #  -) allocate zero sized arrays
            #  -) accumulate up

            # I received a further clarification from Fitsum on 2015-06-24.
            # Namely, for the aggregated durations the annual maximum series is
            # computed from the AGGREGATED data but the correction is applied
            # to the six minute data.  In contrast, ALL STATISTICS are computed
            # and applied to the non-maximum series.  In this case the values
            # that contribute to the maximum series are nulled out.

            # First step is to concatenate the 6 minute time steps to
            # higher if required.  This is a generally simple process in
            # that we always start from the first data measurement in a
            # given day and sum from there.  For example: at 6 minute
            # duration we don't do a thing as we are using all the data
            # points.  For the 60 minute duration there are 10 points to
            # make 60 minutes so we split the 240 data sets into 24 records
            # comprised of indices:
            #   Record  1: Indices 1 - 10
            #   Record  2: Indices 11 - 20
            #   ...
            #   Record 24: Indices 231 - 240

            if TargetIFDduration[currentDur]/minsPerSample != 1:
                # This is NOT the six minute time step so we have to
                # accumulate
                # First: compute the number of records we will be
                # aggregating together for the required duration:
                nRecordsToAggregate = TargetIFDduration[currentDur]/minsPerSample
                
                # Then hand off to a dedicated function now we know the level that
                # we wish to aggregate to
                dataAccumulated = aggregateRainfall(dataRaw, nRecordsToAggregate)
            else:
                dataAccumulated = np.copy(dataRaw)

            # Allocate arrays.
            # The array annualMaxValue is independent of both the chosen
            # durations and the number of simulations but it should be rezeroed
            # at each simulation to ensure a clean slate before the next
            # simulation.
            annualMaxValue = np.zeros((nSims, len(years)))
            annualMaxIndexLinear = np.zeros((nSims, len(years)))

            # We are also cutting down the number of arrays and not holding the
            # non-extreme data in a second array but rather we are keeping a
            # logical index into the main data array that points to the extreme
            # values.
            annualMaxIndexLogical = np.zeros((
                nSims,
                NoDays,
                nRecordsPerDay), 
                dtype=bool)

            # This loop goes through each simulation.  It may be possible to
            # further vectorise this but remember that it will most likely have
            # to be ported to Fortran 2008.
            for loopSim in range(nSims):
                # Explanatory note from Fitsum:  all the correction factors are
                # computed using the AGGREGATED series.  After computation of
                # the correction factors we then ignore the aggregated series
                # and apply the correction factors to the underlying six minute
                # data.  Then for all following calculations, such as the IFD
                # and objective functions, we use only the re-corrected six
                # minute series.
                #
                # Steps:
                #   1) allocate RAM for the ENTIRE data set
                #   2) reshape the array as required to linear
                #   3) fix the max indices to -1, done below inside the year
                #   loop
                #   4) sort and drop the -1 values, do this outside of the loop
                #      below in one hit
                
                # With the accumulated data we can now go through each year and
                # extract the maximum and its index.
                for loopYear in range(len(years)):
                    # Get The EXTREME Series (Annual Maximum)
                    # Here I'm using a logical index to extract only the rows
                    # for a given year.
                    #
                    # I'm doing this in one hit and only extracting the first
                    # value in the case that the maximum value occurs more than
                    # once.
                    #
                    # What the following lines of code do is to take the
                    # accumulation array into a temporary array and then zero
                    # out the years that we are NOT interested in.  That way
                    # when we ask for the index of the maximum element it is
                    # the correct linear index into the ENTIRE data set for
                    # that simulation at the accumulation level.  For durations
                    # not equal to six minutes we will have to multiply that
                    # index by the number of elements we are accumulating at to
                    # get the correct index into the full array.
                    tmpAcummulatedData = np.copy(dataAccumulated[loopSim, :, :])
                    tmpAcummulatedData[years[loopYear] != yearsVector, :] = 0
                    tmpMaxValue = tmpAcummulatedData.max()
                    tmpMaxIndex = tmpAcummulatedData.ravel(order='C').argmax() #index into flattened array

                    annualMaxValue[loopSim, loopYear] = tmpMaxValue

                    # The index needs to be adjusted here to take into account
                    # the removal of the accumulated data.
                    #
                    # Also as we are only dealing with a single simulation
                    # slice of the dataRaw array we need to add additional
                    # lengths to the linear index to take account of the loop
                    # through the simulations.
                    #
                    # NB: these indices are into the series that they were
                    # extracted from.  So for the duration ~= six minutes we
                    # adjust the logical indices to include a stride.  That way
                    # we can zero out the values that contribute to the maximum
                    # series below.  We are ONLY doing the logical index here.
                    # The corrections to the linear index for the annual
                    # maximum stride will be done below at the correction
                    # stage.
                    additionalLength = int(loopSim * nRecordsPerDay * NoDays
                        / (TargetIFDduration[currentDur]/minsPerSample)
                    )
                    annualMaxIndexLinear[loopSim, loopYear] = (
                        tmpMaxIndex + additionalLength
                    )
                    if TargetIFDduration[currentDur]/minsPerSample == 1:
                        # We do not correct for stride and simply add the
                        # wrapping term
                        annualMaxIndexLogical.ravel(order='C')[tmpMaxIndex + additionalLength] = True
                    else:
                        # Here we have to "expand" the single linear index to
                        # include a stride term and the additional length to
                        # account for the aggregated to non-aggregated series
                        correctionStride = int(TargetIFDduration[currentDur]/minsPerSample)

                        correctionIndexLinear = int(
                            (annualMaxIndexLinear[loopSim, loopYear]) *
                            TargetIFDduration[currentDur]/minsPerSample 
                        )
                        # Now convert it to a stride
                        correctionIndexLinear = (
                            np.arange(correctionIndexLinear,correctionIndexLinear+correctionStride, 1)
                        )
                        annualMaxIndexLogical.ravel(order='C')[correctionIndexLinear] = True

                        del additionalLength, correctionStride
                        del correctionIndexLinear

                    del tmpAcummulatedData, tmpMaxIndex, tmpMaxValue

            ##############################################################################
            print(f'Step 2.{currentRec+1}.{currentDur+1}: Sort and Rank Values')
            # Step 2a/b
            # Do the extreme first
            annualMaxSortedValue = np.flip(np.sort(annualMaxValue, axis=1),axis=1)
            annualMaxSortedIndex = np.flip(np.argsort(annualMaxValue, axis=1),axis=1)
            meanExtremeRank = np.mean(annualMaxSortedValue, axis=0) #why calculate the mean across sims?? (caleb)

            # Now the non-extreme
            # Remember this is now based on a logical array.  So we reshape the
            # raw data array and then use the logical mask to drop the
            # extreme values from the non-extreme before we take the mean.  We
            # also drop out the data > realiz as we may not be doing every
            # simulation.
            #
            # NB: this ALWAYS works on the raw data set and NEVER on the
            # aggregated data.  As per email from Fitsum dated 2015-06-24.
            #
            # I'm looping over the number of realiz here but Fortran should be
            # able to vectorize using the logical mask.
            #
            # Fitsum has a rounding function in here that is changing the
            # results if I drop it off.  Following conversation with Fitsum we
            # can drop that rounding function and carry the machine precision
            # through until the final output.
            #
            # Because we ignore the zero rainfall days later, to keep
            # consistent indices rolling through all we have to do is inside
            # the tmpData array is zero out annual maximum values.  We then
            # compute as normal and the indices and sizes will work out.
            tmpData = np.copy(dataRaw).reshape(
                (nSims, nRecordsPerDay*NoDays), order='C')
            tmpLogi = np.copy(annualMaxIndexLogical).reshape(
                (nSims, nRecordsPerDay*NoDays), order='C')

            for loopSim in range(nSims):
                tmpData[loopSim, tmpLogi[loopSim, :]] = 0

            nonExtremeSortedValue = np.sort(tmpData, axis=1)
            nonExtremeSortedIndex = np.argsort(tmpData, axis=1)
            meanNonExtremeRank = np.mean(nonExtremeSortedValue, axis=0)  #why calculate the mean across sims?? (caleb) # to dermine elements that are zero across all sims?
            del tmpData, tmpLogi

            # Before we go on, collapse the arrays to dump the zero elements
            # The maximum value indices still exist in the array they have just
            # been zeroed out. When you take out the indices that have a mean of 0 
            # across all sims you may still leave behind indices of max values if 
            # a non zero value exist for that indice in a different sim. This can
            # cause issues when scaling, as a large scaling factor meant for a 
            # non extreme value can be inadvetently applied to a maximum value.
            nonExtremeSortedIndex = nonExtremeSortedIndex[:, meanNonExtremeRank != 0,].astype(int)
            nonExtremeSortedValue = nonExtremeSortedValue[:, meanNonExtremeRank != 0] 
            meanNonExtremeRank = meanNonExtremeRank[meanNonExtremeRank != 0] 
            meanNonExtremeRankLength = len(meanNonExtremeRank)

            # And add in the additional index so that the linear index stored
            # in nonExtremeSortedIndex is a true linear index into the entire
            # array and not just the slice that was sorted.
            for loopSim in range(1,nSims,1):
                additionalLength = int(loopSim * (nRecordsPerDay * NoDays))
                nonExtremeSortedIndex[loopSim, :] = (
                    nonExtremeSortedIndex[loopSim, :] + additionalLength
                )
            del loopSim

            ##############################################################################
            print(f'Step 3.{currentRec+1}.{currentDur+1}: Compute the Empirical ARI')
            # Step 3a/b
            # Compute the empirical ARI in order to compare it with.
            # I'm not 100% sure what they are doing here but the find
            # operations are relatively slow.  Could this just be done with a
            # direct computation?

            # Extreme First
            indexFrequency = np.zeros((len(Freq), 1))
            for loopFreq in range(len(Freq)):
                indexFrequency[loopFreq] = np.argmin(
                    abs(Freq[loopFreq] - 
                    100 * np.arange(1,len(meanExtremeRank)+1,1) / (len(meanExtremeRank) + 1))
                )
            indexFrequency = indexFrequency.astype(int)

            # Force no correction when there is zero rainfall:
            # not likely for annual maxima series? (caleb)
            meanExtremeRank[meanExtremeRank == 0] = 1

            # Make sure target IFD is ordered correctly (caleb)
            # Why use mean extreme value for each rank? (caleb)                                                           
            fparam = (targetIFD[:, TargetIFDduration[currentDur]==TargetIFDdurationEst]
                / meanExtremeRank[indexFrequency[::-1]]
            )

            # Fitsum's code does a polynomial fit here
            # why to power 3? has this been sensitivty tested? (Caleb)
            # fitting a polynomial by only several points? (Caleb)                                                                                
            extremePoly = np.polyfit(np.squeeze(meanExtremeRank[indexFrequency[::-1]]),
                np.squeeze(fparam), 3
            )
            extremeFit = (
                extremePoly[0] * meanExtremeRank**3 +
                extremePoly[1] * meanExtremeRank**2 +
                extremePoly[2] * meanExtremeRank +
                extremePoly[3]
            )
            extremeFit[extremeFit <= 0] = 0

            # Non-Extreme
            # The next block of equations compute the the five "terms" in
            # Equation 8 of the draft paper.  Fitsum's original code contained
            # an extra find loop to only compute the non-zero terms in the
            # non-extreme series.  This works because despite the rank, any
            # term multiplied by zero rainfall results in a zero contribution
            # to the cumulative sum.  However, in MATLAB and Fortran the find
            # operation is an order or more of magnitude slower than
            # vectorising the entire series and computing the sum.  Hence, I'm
            # going to vectorise and sum.
            #
            # Define: Ri = non-extreme rainfall; Rj = extreme rainfall
            #
            # Term 1: sumRiRj in R
            # First compute the total depth and hold it to maintain the overall
            # average after the correction.
            if currentRec == 0 and currentDur == 0:
                meanNonExtremeSum = np.sum(meanNonExtremeRank)
                meanExtremeSum = np.sum(meanExtremeRank)

            termOne = (meanNonExtremeSum + meanExtremeSum) * massScale

            # Term 2: sumRiWeightedNom in R
            termTwo = (extremeFit[-1] *
                np.sum((np.arange(0,meanNonExtremeRankLength,1)) / 
                meanNonExtremeRankLength * meanNonExtremeRank
                ))

            # Term 3: sumRjcor in R
            # This following line is my original implementation
            # termThree = sum(meanExtremeRank .* extremeFit);
            # This is from Fitsum's R
            termThree = 0
            #for loopSim in np.arange(0,np.size(dataRaw,2),1):
            for loopSim in range(nSims):    
                termThree = (termThree + np.sum(annualMaxValue[loopSim, :] *
                    extremeFit[annualMaxSortedIndex[loopSim, :]])
                )
            termThree = termThree / nSims

            # Term 4 (denominator): sumRiWeightedDenom in R
            termFour = np.sum(
                ((meanNonExtremeRankLength+1) - np.arange(1,meanNonExtremeRankLength+1,1))
                / meanNonExtremeRankLength * meanNonExtremeRank
            )

            # Compute Equation 8
            fOneNoEx = (termOne - termTwo - termThree) / termFour

            # Now correct fOneNoEx incase it is negative
            if fOneNoEx >= 0:
                # This is Equations 5 and 6:
                delta = (extremeFit[-1] - fOneNoEx) / meanNonExtremeRankLength
                fnNonEx = fOneNoEx + np.arange(0,meanNonExtremeRankLength,1) * delta
            else:
                # It is negative.  This method is based on equations 9 through
                # 11 in the draft paper.
                #
                # First compute the LHS of equation 11.  We know all these
                # values as they are the pre-correction values.
                lhs = (meanNonExtremeSum + meanExtremeSum) * massScale

                # Next, compute the second term on the RHS of the equation 11
                # as this is the corrected extreme values, which we did above.
                rhsExtreme = np.sum(meanExtremeRank * extremeFit)

                # Now do a minimisation operation to estimate "b".
                # changed absolute value brackets (Caleb 01/09/2022)
                def eqll(x):
                    return (
                        abs(lhs - extremeFit[-1]/(meanNonExtremeRankLength**x)  
                        * np.sum((np.arange(1,meanNonExtremeRankLength+1)**x)
                        * meanNonExtremeRank) - rhsExtreme)
                    )
                # Then minimise:
                b, fval,_ ,_ , flag = scipy.optimize.fmin(func=eqll, x0=0, full_output=True)
                # Now back substitute into Equation 12 to directly compute "a":
                a = extremeFit[-1] / (meanNonExtremeRankLength ** b)  

                # Finally, with a and b compute Equation 9:
                fnNonEx = a * np.arange(1,meanNonExtremeRankLength+1,1) ** b

            if fOneNoEx > 0:
                rhsExtreme = np.sum(meanExtremeRank * extremeFit)
            nonExtremeSum_scaled = np.sum(meanNonExtremeRank * fnNonEx)

            del nonExtremeSum_scaled
            ##############################################################################
            print(f'Step 4.{currentRec+1}.{currentDur+1}: Correct the Rainfall Sequences')
            # Step 4: Correct The Rainfall Sequences
            #
            # The process is, conceptually, quite simple: multiply rain by
            # its correction factor.
            #
            # However, the difficulty is that the correction factors must be
            # multiplied based on the rank order of the rainfall but must then
            # be inserted back into the time sequence.  In order to do this
            # efficiently for every sort operation performed above we have kept
            # the indices of the sort operation.  Therefore, while it may look
            # complex we are slicing back through the sort operations to
            # multiply the data in the original array.
            #
            # The general steps are:
            #   -) loop over simulations:
            #       -) apply corrections to extreme data
            #       -) apply corrections to the non-extreme data
            #
            # NB: I'm repeating the comments here so that I REALLY don't forget
            # that the linear indexes are into aggregated time series.  For six
            # minute durations no further action is required.  However, for
            # other durations where we aggregate up we must adjust the linear
            # indices so that we correct not the aggregated series but the
            # UNDERLYING six minute series.

            # Now Correct

            dataRaw_flat = dataRaw.flatten(order='C')
            dataRaw_nomax = np.copy(dataRaw)
            dataRaw_nomax[annualMaxIndexLogical] = 0
            dataRaw_nomax_flat = dataRaw_nomax.flatten(order='C')

            dataRaw_flat = correction(minsPerSample,
                            nSims,
                            currentDur,
                            years,
                            TargetIFDduration,
                            annualMaxIndexLinear,
                            annualMaxSortedIndex,
                            dataRaw_flat,
                            dataRaw_nomax_flat,
                            extremeFit,
                            nonExtremeSortedIndex,
                            fnNonEx
                            )
            dataRaw = dataRaw_flat.reshape((nSims, nDays, nRecordsPerDay), order='C')
            del dataRaw_flat

            ##############################################################################
            print(f'Step 5.{currentRec+1}.{currentDur+1}: Test the Correction Skill')
            # Step 5: Test the correction skill
            # The process is based on comparing our computed IFD from the
            # corrected data against a known target IFD.
            IFDCorrected = computeIFD(dataRaw, yearsVector, TargetIFDduration)

            # Now compare our corrected IFD with that of the target via the
            # objective function
            objectiveResults[currentRec, currentDur, :] = ifdcondobjfun(targetIFD, IFDCorrected, Freq)

            # Extract debugging data sets:
            if debugMode == True:
                tmp = np.reshape(dataRaw, (len(dateVector) * nRecordsPerDay, nSims), order='F')
                nRainDaysCor[currentRec, currentDur, :] = np.sum(tmp.astype(bool), axis=0)
                rainfallDepthCor[currentRec, currentDur, :] = np.sum(tmp, axis=0)
                meanAnnualMaxCor[currentRec, currentDur] = np.mean(max(np.reshape(
                    dataRaw, (nRecordsPerDay*len(dateVector), nSims), order='F')),axis=0)
                del tmp

    if plot == True:
        fig, ax = plt.subplots(nrows=len(TargetIFDduration), ncols=1, figsize=(10,len(TargetIFDduration)*5))
        for i in range(len(TargetIFDduration)):
            # Raw IFD
            if i == len(TargetIFDduration)-1:
                ax[i].set_xlabel('Frequency (years)', fontsize=14)
            else:
                ax[i].set_xlabel('')
            ax[i].set_ylabel('Depth (mm)', fontsize=14)
            ax[i].fill_between(np.arange(len(years)), 
                                IFDRaw[:,:,i].max(axis=1),
                                IFDRaw[:,:,i].min(axis=1),
                                alpha=0.1, label='Raw Simulations Range')
            ax[i].plot(np.median(IFDRaw[:,:,i], axis=1),
                        linewidth=2, alpha=0.8, 
                        label=f'Raw Median Depth-Frequency {TargetIFDduration[i]} min')
            #Conditioned IFD
            ax[i].fill_between(np.arange(len(years)), 
                                IFDCorrected[:,:,i].max(axis=1),
                                IFDCorrected[:,:,i].min(axis=1),
                                alpha=0.1, label='Cond Simulations Range')
            ax[i].plot(np.median(IFDCorrected[:,:,i], axis=1),
                        linewidth=2, alpha=0.8, 
                        label=f'Cond Median Depth-Frequency {TargetIFDduration[i]} min')
            ax[i].set_xscale('log')
            if i == 0:
                ax[i].legend()
        fig.savefig('ifdcond_plot.png')
    
    ## Dump the Processed Time Series to a NetCDF
    print('Saving conditioned data')
    dataRaw = np.swapaxes(dataRaw,0,2)
    produceSubDailyNetCDF(
        fileNameOutput,
        dataRaw,
        dayVector,
        title = 'Example Conditioned Simulations for Perth',
        institution = 'UNSW Water Research Centre',
        nRecordsPerDay = nRecordsPerDay
    )
    print('Saving conditioned data: done')