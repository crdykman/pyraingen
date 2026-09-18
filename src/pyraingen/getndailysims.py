import netCDF4 as nc

def getNDailySims(fnameNC):
    """Function to get the number of daily simulations 
    Get the number of daily simulations ouf of a specifically formatted
    NetCDF file.  The NetCDF would adhere to the structure defined in:
    utils/convert_national_daily_to_NetCDF.m.
     
    Parameters
    ----------
    fnameNC : str
        The full path to the file to read data from.\n

    Returns
    ----------
    int 
        Number of daily simulations within the file fnameNC.
    """
    #    Dr Peter Brady <peter.brady@wmawater.com.au>
    #    2016-08-21
    #
    #    Caleb Dykman
    #    2021-09-14

    ## Set default
    nDailySims = 0

    ## Local Parameters
    targetVarName = 'simulation'

    ## Working
    ds=nc.Dataset(fnameNC)

    for key in ds.dimensions.keys():
        if targetVarName == key:
            nDailySims = len(ds.dimensions[key])
    ds.close()

    return int(nDailySims)
