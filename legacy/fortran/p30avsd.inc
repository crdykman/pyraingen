	real*4  pp1av(monmax,ndmax,nvrmax),
     +          pp2av(monmax,ndmax,nvrmax),
     +          pp3av(monmax,ndmax,nvrmax),
     +          pp4av(monmax,ndmax,nvrmax),
     +        pp1cov(monmax,ndmax,nvrmax,nvrmax),
     +        pp2cov(monmax,ndmax,nvrmax,nvrmax),
     +        pp3cov(monmax,ndmax,nvrmax,nvrmax),
     +        pp4cov(monmax,ndmax,nvrmax,nvrmax)
  
	real*4  p1av(monmax,ndmax,nvrmax),
     +          p2av(monmax,ndmax,nvrmax),
     +          p3av(monmax,ndmax,nvrmax),
     +          p4av(monmax,ndmax,nvrmax),
     +        p1cov(monmax,ndmax,nvrmax,nvrmax),
     +        p2cov(monmax,ndmax,nvrmax,nvrmax),
     +        p3cov(monmax,ndmax,nvrmax,nvrmax),
     +        p4cov(monmax,ndmax,nvrmax,nvrmax),
     +        almr(5,monmax,ndmax)
	real*4 pprev1(monmax,ndmax),
     +         pprev2(monmax,ndmax),
     +         pprev3(monmax,ndmax),
     +         pprev4(monmax,ndmax),
     +         prev1(monmax,ndmax),
     +         prev2(monmax,ndmax),
     +         prev3(monmax,ndmax),
     +         prev4(monmax,ndmax)

	real*4  acl(nstnmax)

	real*4  pp1lav(monmax,ndmax,nvrmax),
     +          pp2lav(monmax,ndmax,nvrmax),
     +          pp3lav(monmax,ndmax,nvrmax),
     +          pp4lav(monmax,ndmax,nvrmax),
     +        pp1lcov(monmax,ndmax,nvrmax,nvrmax),
     +        pp2lcov(monmax,ndmax,nvrmax,nvrmax),
     +        pp3lcov(monmax,ndmax,nvrmax,nvrmax),
     +        pp4lcov(monmax,ndmax,nvrmax,nvrmax)
  
	real*4  p1lav(monmax,ndmax,nvrmax),
     +          p2lav(monmax,ndmax,nvrmax),
     +          p3lav(monmax,ndmax,nvrmax),
     +          p4lav(monmax,ndmax,nvrmax),
     +        p1lcov(monmax,ndmax,nvrmax,nvrmax),
     +        p2lcov(monmax,ndmax,nvrmax,nvrmax),
     +        p3lcov(monmax,ndmax,nvrmax,nvrmax),
     +        p4lcov(monmax,ndmax,nvrmax,nvrmax)


	common /p30_avsd/ pp1av,pp2av,pp3av,pp4av,
     +                    pp1cov,pp2cov,pp3cov,pp4cov,
     +	                  p1av,p2av,p3av,p4av,pprev1,
     +                    pprev2,pprev3,pprev4,prev1,prev2,
     +                    prev3,prev4,acl,
     +                    p1cov,p2cov,p3cov,p4cov,almr,
     +	                  pp1lav,pp2lav,pp3lav,pp4lav,
     +                    pp1lcov,pp2lcov,pp3lcov,pp4lcov,
     +	                  p1lav,p2lav,p3lav,p4lav,
     +                    p1lcov,p2lcov,p3lcov,p4lcov
