import os
import numpy as np
import netCDF4 as nc
import datetime

def produceDailyNetCDF(fnameNC, DailyData, dayVector,**kwargs):
    """Creates from scratch the daily netCDF as per the schema.
    This function is intended to be single write point for daily data so
    that all data written uses the same documented netCDF schema.

    Parameters
    ----------
    fnameNC : str
        The name of the netCDF file to dump data into.  If the
        file exists, a warning will be thrown and the fuction return an
        error.
    dailyData : array
        An array of the daily read data with the first
        dimension being the reading and the second the simulation.  For
        BoM or other "real" data, as opposed to simulated data, the
        simulation dimension will == 1.  Data may be padded with BoM or
        similar negative numbers.
    dayVector : array
        a vector of Julian dates for the day of the
        measurement recording.  For daily data this is set to midnight
        and remember that Julian days start at midday.  No timezone
        information is stored.
    title : str 
        (optional) Data set title to be inserted into the global attribute
        "title".  Default: "Daily Rainfall".
    institution : str 
        a text field to identify the institution that
        created this data set in the global attribute "institution".
        Default:
            'UNSW Water Research Centre'
    
    Returns
    ----------
    Bool
        status: true for OK, false for an error.
    """
    #   Dr Peter Brady <peter.brady@wmawater.com.au>
    #   2016-08-26

    #  Caleb Dykman
    #  2022-07-15

    ## Check of the file exists
    if os.path.exists(fnameNC):
        raise Warning('PSDNC:Exists ' + 'The file {}, exists please check the '
        'name and/or delete the file and try again.'.format(fnameNC))
        status = False
        return status
    
    ## Set some defaults
    status = True
    dataSetTitle = 'Daily Rainfall'
    dataSetInstitution = 'UNSW Water Research Centre'
    
    ## Get Some Metadata
    # First find the revision.  NB: this uses Subversion keywords.  Do not
    # change the following two lines.
    # revText = '$Revision:: 188   $';
    # headURLText = '$HeadURL: https://develop.wmawater.com.au/svn/arr_p4/branches/pb_working/utils/produceDailyNetCDF.m $';
    # revision = sscanf(revText(12:1:18), '%d');
    # headURL = sscanf(headURLText, '$HeadURL: %s$');
    # [~, username] = system('whoami');
    # username = username(1:1:(end-1));
    # pcinfo=java.net.InetAddress.getLocalHost;
    # hostname = char(pcinfo.getHostName);

    # Unpack the kwargs
    if len(kwargs) > 0:
        for key in kwargs.keys():
            if key.lower() == 'title':
               dataSetTitle = kwargs[key] 

            if key.lower() == 'institution':
               dataSetInstitution = kwargs[key]

    ## Now Write
    # Create NetCDF
    Dailync = nc.Dataset(fnameNC, 'w', format='NETCDF4')
   
    # Define Dimensions
    Dailync.createDimension('day',len(dayVector))
    Dailync.createDimension('simulation',
        np.size(DailyData,axis=1)
        )
    
    
    # Create Variables
    days = Dailync.createVariable('day','f8',('day',),zlib=True,
        complevel=9, shuffle=True
    )
    rainfalls = Dailync.createVariable('rainfall','f4',
        ('day','simulation',), zlib=True, complevel=9, 
        shuffle=True, fill_value=-999.9
    )

    ## Writing Data
    days[:] = dayVector
    rainfalls[:,:] = DailyData 

    ## Add Attributes
    # Global
    Dailync.title = dataSetTitle
    Dailync.history = 'Generated {} by {}@{}'.format(datetime.datetime.now()
        ,'username','hostname'
    )
    Dailync.source = 'Python function: {}@{}'.format(
        'headURL','revision'
    )
    Dailync.institution = dataSetInstitution
    Dailync.conventions = 'ARR Project 4'
    
    # Variable
    days.units = 'Julian date, no timezone, 1900-01-01 00:00:00 == 2415020.5'
    days.long_name = 'Julian Date'
    rainfalls.units = 'mm'
    rainfalls.long_name = 'Daily Rainfall'

    Dailync.close()

    return status

    
