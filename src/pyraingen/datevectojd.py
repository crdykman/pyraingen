import numpy as np
from datetime import date

def datevecToJD(dateVector):
    """This function takes a formatted dateVector, as outputed from datevec or
    similar, and converts it to a Julian Date.  There is no timezone
    support.
    
     NB: Julian days start not at midnight but midday.
    
    This implementation is based on the integer maths algorithm of Fliegel,
    H. F. & van Flandern, T. C. 1968, Communications of the ACM, 11, 657
    with modifications from:
           http://aa.usno.navy.mil/faq/docs/JD_Formula.php
    It is extended to handle fractional days with hours, minutes and
    seconds.
    
    Parameters
    ----------
    dateVector : array
        A vector whos elements are:\n
        -) year     (integer)\n
        -) month    (integer)\n
        -) day      (integer)\n
        -) hour     (integer)\n
        -) minute   (integer)\n
        -) second   (double)\n

    Returns
    ----------
    float
        Resultant Julian Date.
    """
    #  Dr Peter Brady <peter.brady@wmawater.com.au>
    #  2016-08-26

    #  Caleb Dykman
    #  2021-09-27

    # Check necessary for the USNO algorithm
    if dateVector.year < 1801 or dateVector.year > 2099:
        raise ValueError('DTJD:YearRange'
        'Out of year range, this algorithm is restricted to the '
        'range: 1801 <= year <= 2099')
    
    year   = dateVector.year #dateVector[0]
    month  = dateVector.month #dateVector[1]
    day    = dateVector.day #dateVector[2]
    hour   = 0 #dateVector[3]
    minute = 0 #dateVector[4]
    second = 0 #dateVector[5]
    subDay = hour / 24 + minute / (24 * 60) + second / (24 * 60 *60)

    # Fliegel and van Flandern does not work for MATLAB due to the way MATLAB
    # handles integer truncation.  It also does not handle fractional days.
    # julianDate = DAY - 32075 + 1461 * (Year + 4800 + (Month - 14) / 12) / 4 +...
    #     367 * (Month - 2 - (Month - 14) / 12 * 12) / 12 - 3 * ...
    #     (( YEAR + 4900 + (Month - 14) / 12) / 100) / 4;

    # USNO Implementation valid for 1801 <= year <= 2099 (modified)
    julianDate = (367 * year - np.fix(7 * (year + np.fix((month + 9) / 12)) / 4)
        + np.fix((275 * month)/9) + day + 1721013.5 + subDay
        - 0.5 * np.sign(100 * year + month - 190002.5) + 0.5
    )

    return julianDate