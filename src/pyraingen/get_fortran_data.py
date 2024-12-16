import shutil
from importlib import resources
import os

def copy_fortran_data():
    """Copy necessary fortran input data to working directory..
    """
    with resources.path("pyraingen.fortran_daily", "data.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "data.inc")

    with resources.path("pyraingen.fortran_daily", "occr.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "occr.inc")

    with resources.path("pyraingen.fortran_daily", "rf_amt.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "rf_amt.inc")
    
    with resources.path("pyraingen.fortran_daily", "pcond.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "pcond.inc")

    with resources.path("pyraingen.fortran_daily", "pminmax.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "pminmax.inc")

    with resources.path("pyraingen.fortran_daily", "occrcor.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "occrcor.inc")

    with resources.path("pyraingen.fortran_daily", "p30avsd.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "p30avsd.inc")

    with resources.path("pyraingen.fortran_daily", "stat.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "stat.inc")

    with resources.path("pyraingen.fortran_daily", "para.inc") as f:
        src_pth = f
    shutil.copy(src_pth, "para.inc")

    if not os.path.exists("stn_record.dat"):
        with resources.path("pyraingen.fortran_daily", "stn_record.dat") as f:
            src_pth = f
        shutil.copy(src_pth, "stn_record.dat")

    if not os.path.exists("data_r.dat"):
        with resources.path("pyraingen.fortran_daily", "data_r.dat") as f:
            src_pth = f
        shutil.copy(src_pth, "data_r.dat")

    # print("Copied")