	real*4  pp1lav(monmax,ndmax,nvmax,nstnmax),
     +          pp2lav(monmax,ndmax,nvmax,nstnmax),
     +        pp1lcov(monmax,ndmax,nvvmax,nstnmax),
     +        pp2lcov(monmax,ndmax,nvvmax,nstnmax)
  
	real*4  p1lav(monmax,ndmax,nvmax,nstnmax),
     +          p2lav(monmax,ndmax,nvmax,nstnmax),
     +        p1lcov(monmax,ndmax,nvvmax,nstnmax),
     +        p2lcov(monmax,ndmax,nvvmax,nstnmax)

	common /p_cond/ pp1lav,pp2lav,pp1lcov,pp2lcov,
     +	                  p1lav,p2lav,p1lcov,p2lcov
