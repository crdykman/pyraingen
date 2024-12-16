from numba import njit

@njit
def getDaysPerYear(year):
    """This function is designed to return the number of days in a given year
    to handle leap years.  It has been validated for the years 1800 to
    2100.  Dates outside that range should be treated with caution.
    
    Parameters
    ----------
    year : int
        Year to get the dates from in the Gregorian calender,
        preferably between 1800 and 2100.\n

    Returns
    ----------
    int
        The days in the requested year.
    """
    #   Dr Peter Brady <peter.brady@wmawater.com.au>
    #   2016-08-16

    #   Caleb Dykman
    #   2021-09-16

    ## Local Parameters
    # Days in the year
    daysYearLeap = 366
    daysYear = 365

    # Validated year window
    validYearStart = 1800
    validYearEnd   = 2100

    ## Input Check
    if year < validYearStart  or year > validYearEnd:
        print(
            f'getDaysPerYear has not been validated for the year {year}. Only over the year range ' +
            str(validYearStart) +
            ' - ' +
            str(validYearEnd)
        )


    ## Logic
    # From: https://en.wikipedia.org/wiki/Leap_year#Algorithm (2016-08-14)
    if year % 4 != 0:
        daysPerYear = daysYear
    elif year % 100 != 0:
        daysPerYear = daysYearLeap
    elif year % 400 != 0:
        daysPerYear = daysYear
    else:
        daysPerYear = daysYearLeap

    return daysPerYear


        
        


