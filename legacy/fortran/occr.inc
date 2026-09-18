	real*4   covar(monmax,ndmax,nstnmax,nstnmax),
     +  pro1(monmax,ndmax,nstnmax),pro2(monmax,ndmax,nstnmax),
     +  pro(monmax,ndmax,nstnmax),high(nlmax,nyrmax,monmax,ndmax)
	real*4 avl(nvmax,monmax,ndmax,nstnmax)
	real*4 sdl(nvmax,monmax,ndmax,nstnmax)
  
	common /occ_cor/ covar,pro1,pro2,pro,high,avl,sdl

