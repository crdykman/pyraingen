	real*4 rfcur(nearmax1,monmax,ndmax,nstnmax)
	real*4 rfprv(nearmax1,monmax,ndmax,nstnmax)
	real*4 rflcl(nearmax1,monmax,ndmax,nstnmax)
	real*4 stnwt1(nstnmax,monmax,ndmax,nvmax)
	real*4 stnwt2(nstnmax,monmax,ndmax,nvmax)
	real*4 stnwt3(nstnmax,monmax,ndmax,nvmax)
	real*4 stnwt4(nstnmax,monmax,ndmax,nvmax)
	real*4 cov_amt(monmax,ndmax,nstnmax,nstnmax)
	integer irfnxt(nearmax1,monmax,ndmax,nstnmax)
	common /rfamount/ rfcur,rfprv,stnwt1,stnwt2,
	1     stnwt3,stnwt4,cov_amt,irfnxt,rflcl