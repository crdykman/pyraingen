c	nstnmax-maximum number of stations
c	nyrmax=maximum number of years of data
c	monmax=maximum number of months
c	ndmax=maximum number of days in a month
c	ndatmax=total number of data(days)
c	isnmax -maximum number of seasons in a year
c	ndsmax- max number of observations including station values=ndatmax*nstnmax
c	nvmax- maximum no of parameters or variables
c	ndaymax - total number of days in a year/season
c	nearmax-maximum number of possible patterns considered (number of years*length of window)
c	nlagmax-maximum number of lags
c	lonmax -maximum number of previous days considered

	integer nyrmax,monmax,ndmax,nvmax,nstnmax
	integer ndatmax,nearmax,nvctmax,nlagmax,ndaymax,nearmax1
	parameter (nyrmax=200,monmax=12,ndmax=31)
	parameter (nstnmax=6,isnmax=5,itrmax=100)
	parameter (mvmax=62,nearmax=mvmax*nyrmax,nvctmax=10)
	parameter (ndatmax=nyrmax*monmax*ndmax/3,nlagmax=3)
	parameter (ndsmax=ndatmax*nstnmax/2,ndaymax=ndmax*monmax)
	parameter (nvrmax=5,nvvrmax=nvrmax*nvrmax,lonmax=900)
	parameter (ijlmax=nstnmax,nlmax=nstnmax*nvrmax,nvmax=nvrmax)
	parameter (nearmax1=int(nearmax/2),nvvmax=nvmax*nvmax)
	