	real*4   rf(nstnmax,nyrmax,monmax,ndmax)
	real*4   rain
	integer  nday(monmax)
	integer  nstn,nout,iband,ndata,nvar,nsim
	integer  nrain
	integer  monst,moned,iprf,ng,nsg,llag,lav,iwet
	integer lon(2,nvrmax),nlon,nprev,nprev1,ijl,iyrst,iseed,isn
	integer npp(5,nstnmax),npr(5,nstnmax),iamt,nlav(nstnmax)
	integer nlonr,ispl,nlocal,nyrs(nstnmax),nstrt(nstnmax)

	common /data/ rf,nstn,nout,iband,ndata,iprf,
     1          nrain,rain,nday,nvar,nsim,monst,moned,
     1          nprev,ng,nsg,gend,llag,lav,iwet,lon,nlon,nlav,
     1          nprev1,ijl,iyrst,iseed,isn,npp,iat,npr,
     1          iamt,nlonr,ispl,nlocal,nyrs,nstrt