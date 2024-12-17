import shutil
from importlib import resources
import os

def get_example_data(dname=None, daily=True, subdaily=True, ifd=True):
    """Copy example data to chosen directory.

    Parameters
    ----------
    dname : str
        directory name to copy data to.
        Default is None.
    daily : bool
        True to copy daily example data to specified directory.
        Default is True.
    subdaily : bool
        True to copy subdaily example data to specified directory.
        Default is True.
    ifd : bool
        True to copy ifd example data to specified directory. 
        Default is True.       

    """
    if dname == None:
        dname = os.getcwd()
    
    if daily == True:
        with resources.path("pyraingen.data.example.daily", "rev_dr061003.txt") as f:
            src_pth = str(f)
        source_dir = src_pth[:-16]
        shutil.copytree(source_dir, dname+"/daily")
    
    if subdaily == True:
        with resources.path("pyraingen.data.example.subdaily", "daily.nc") as f:
            src_pth = str(f)
        source_dir = src_pth[:-8]
        shutil.copytree(source_dir, dname+"/subdaily")
    
    if ifd == True:
        with resources.path("pyraingen.data.example.ifd", "subdaily.nc") as f:
            src_pth = str(f)
        source_dir = src_pth[:-11]
        shutil.copytree(source_dir, dname+"/ifd")

    print("Copied")