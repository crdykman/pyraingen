	real*4 solh(nstnmax),solg(itrmax,nstnmax),sol(nstnmax)
	real*4 spellh(nstnmax,15,nyrmax),spell(nstnmax,15,nyrmax),
	1          spellg(itrmax,nstnmax,15,nyrmax),wmaxh(nstnmax),
     1         dmaxh(nstnmax),wmax(nstnmax),dmax(nstnmax),
	1          wmaxg(itrmax,nstnmax),dmaxg(itrmax,nstnmax),
	1          rfavm(nstnmax,monmax),rfsdm(nstnmax,monmax),
     1          roavm(nstnmax,monmax),rosdm(nstnmax,monmax),
     1          rwavm(nstnmax,monmax),rwsdm(nstnmax,monmax),
     1          rfavs(nstnmax,isnmax),rfsds(nstnmax,isnmax),
     1          roavs(nstnmax,isnmax),rosds(nstnmax,isnmax),
     1          rwavs(nstnmax,isnmax),rwsds(nstnmax,isnmax),
	1          rfavmh(nstnmax,monmax),rfsdmh(nstnmax,monmax),
	1          roavmh(nstnmax,monmax),rosdmh(nstnmax,monmax),
	1          rwavmh(nstnmax,monmax),rwsdmh(nstnmax,monmax),
     1          rfavsh(nstnmax,isnmax),rfsdsh(nstnmax,isnmax),
     1          roavsh(nstnmax,isnmax),rosdsh(nstnmax,isnmax),
     1          rwavsh(nstnmax,isnmax),rwsdsh(nstnmax,isnmax),
	1 rfavmg(itrmax,nstnmax,monmax),rfsdmg(itrmax,nstnmax,monmax),
     1 roavmg(itrmax,nstnmax,monmax),rosdmg(itrmax,nstnmax,monmax),
     1 rwavmg(itrmax,nstnmax,monmax),rwsdmg(itrmax,nstnmax,monmax),
     1 rfavsg(itrmax,nstnmax,isnmax),rfsdsg(itrmax,nstnmax,isnmax),
     1 roavsg(itrmax,nstnmax,isnmax),rosdsg(itrmax,nstnmax,isnmax),
     1 rwavsg(itrmax,nstnmax,isnmax),rwsdsg(itrmax,nstnmax,isnmax),
     1 hper(nstnmax,5),hperh(nstnmax,5),hperg(itrmax,nstnmax,5)

	real*4 rfspl(nstnmax,15),rfspn(nstnmax,15),rfspnd(nstnmax,15),
     1      rfsplh(nstnmax,15),rfspnh(nstnmax,15),rfspndh(nstnmax,15),
     1      rfsplg(itrmax,nstnmax,15),rfspng(itrmax,nstnmax,15),
     1      rfspndg(itrmax,nstnmax,15)
	real*4 crsodh(isnmax,nstnmax,nstnmax),
     1    crsoh(isnmax,nstnmax,nstnmax),
     1    crsrh(isnmax,nstnmax,nstnmax),
     1    crsrdwh(isnmax,nstnmax,nstnmax),
     1    crsrdah(isnmax,nstnmax,nstnmax)
	real*4 crsod(isnmax,nstnmax,nstnmax),crso(isnmax,nstnmax,nstnmax),
     1     crsr(isnmax,nstnmax,nstnmax),crsrdw(isnmax,nstnmax,nstnmax),
     1      crsrda(isnmax,nstnmax,nstnmax)
	real*4 crsodg(itrmax,isnmax,nstnmax,nstnmax),
     1         crsog(itrmax,isnmax,nstnmax,nstnmax),
     1         crsrg(itrmax,isnmax,nstnmax,nstnmax),
     1         crsrdwg(itrmax,isnmax,nstnmax,nstnmax),
     1         crsrdag(itrmax,isnmax,nstnmax,nstnmax)


	common/stat_hg/solh,sol,solg,spellh,spell,spellg,dmaxh,dmax,dmaxg,
     1      wmaxh,wmax,wmaxg,rfavmh,rfsdmh,roavmh,rosdmh,rwavmh,
     1      rwsdmh,rfavsh,rfsdsh,roavsh,rosdsh,rwavsh,rwsdsh,rfavm,
     1      rfsdm,roavm,rosdm,rwavm,rwsdm,rfavs,rfsds,roavs,rosds,
     1      rwavs,rwsds,rfavmg,rfsdmg,roavmg,rosdmg,rwavmg,rwsdmg,
     1      rfavsg,rfsdsg,roavsg,rosdsg,rwavsg,rwsdsg,hperh,hper,
     1      hperg,rfsplh,rfspl,rfsplg,rfspnh,rfspn,rfspng,rfspndh,
     1      rfspnd,rfspndg,crsodh,crsoh,crsrh,crsrdwh,crsrdah,crsodg,
     1      crsog,crsrg,crsrdwg,crsrdag,crsod,crso,crsr,crsrdw,crsrda