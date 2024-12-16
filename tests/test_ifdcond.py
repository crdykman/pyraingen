from pyraingen.ifdcond import ifdcond

fileNameInput ="subdaily.nc"
fileNameOutput = "condsubdaily.nc"
fileNameTargetIFD = "targetifds.csv"
TargetIFDdurationsEst = [30, 60, 360, 720]
TargetIFDdurations = [30, 60, 360, 720]
AEP = [63.20, 50, 20, 10, 5, 2]

ifdcond(fileNameInput, fileNameOutput, fileNameTargetIFD, 
        nSims=10, TargetIFDdurationsEst=TargetIFDdurationsEst,
       TargetIFDdurations=TargetIFDdurations, AEP=AEP, plot=True)