      subroutine regionalised_daily(idrop)
c	if idrop=-1, program is called first time, return nearby station information
c	if idrop >0, modify nearby station information
c	if idrop=0, start generation
c	generate multi-site rf occurrences using wilks approach
c	considers previous specified day(s) average wetness state
c	requires correlation matrix of normal deviates
c	considers rainfall series at multiple stations
C	Based on today's rainfall conditioned on previous day's rainfall
c	programme is based on Wilks approach
c
c	if rainfall<0.3 mm its a dry day otherwise a wet day
c

c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      parameter (nsmax=5000)
      character ch*2,file_rf_out*80,cc*6,nearby*80,filet*80
      character fpath*80
      real*4    w(nstnmax),rfh(nstnmax,nyrmax,monmax,ndmax)
      real*4    avrf(nstnmax)
      integer   idx(nstnmax),idr(nsmax)
c	supply values of some fixed parameters here
      isuccess=0
      nout=12
      iseed=131
      isn=4
      
      do i=1,500
      ca= gasdev(iseed)
      enddo
      
      open(1, file='data_r.dat',status='old')
      read(1,*)
      read(1,*)
      read(1,*)rain,iband,nstn,nsim,nlon,llag,ng,nsgst,iamt,ival,irf
      read(1,*)
      read(1,*)(lon(1,nk),nk=1,nlon)
      read(1,*)(lon(2,nk),nk=1,nlon)
      read(1,*)
      Read(1,'(a)')nearby
      Read(1,*)
      Read(1,*)
      read(1,*)idxt,tlat,tlon,tele,tdcoast,tanrf,ttemp
      read(1,*)
      read(1,'(a)')fpath
c     if(ival.gt.0)then
      read(1,*)
      Read(1,'(a)')filet
c     endif
      close(unit=1)
      nsg=nsgst-1
      
      if(nlon.eq.0)then
      iyrst=1
      else
      iyrst=0
      do nk=1,nlon
      if(lon(2,nk).gt.iyrst)iyrst=lon(2,nk)
      enddo
      iyrst=int(float(iyrst)/365.0+2.001)
      endif

      open(13,file='t1.out')
c	define number of days in a month
      call day(nday,monmax)
c	read target location lat, lon, elevation, dist to coast and average annual rainfall
c	As an example following stations are considered
c     66037 SYDNEY AIRPORT AMO
c     46106 TINAPAGEE
c     68085 NERRIGA (TOLWONG)      

c	idxt=046106               
c	tlat=-29.4670
c	tlon=144.3830
c	tele=125.00
c	tdcoast=715.34
c	tanrf=308.91	

c	idxt=066037
c	tlat=-33.9410
c	tlon=151.1730
c	tele=6.00
c	tdcoast=0.27
c	tanrf=1086.71	

c	idxt=068085               
c	tlat=-34.8480
c	tlon=150.1350
c	tele=600.00
c	tdcoast=42.80
c	tanrf=823.37	

c	if annual rainfall at target station is not available read from the Bom grided file
c      if(irf.eq.0.or.tanrf.eq.0.0)call find_anrf(tlat,tlon,tanrf,nsmax)
c	check if any nearby stations are to be ignored
      if(idrop.gt.0)then
      open(98,file='drop.out',status='old')
      do kk=1,idrop
      read(98,*)i1,idr(kk)
      enddo
      close(unit=98)
      kk=idrop
      endif

 17   open(21,file=nearby)
c	find nearest 'nstn' no of stations from the target location
      call station(kk,idr,tlat,tlon,tele,tdcoast,ttemp,tanrf,
     1             idx,w,nstn,nsmax,nstnmax)

c	read data of nearby stations from the given file
      write(*,*)'Reading data'
      do jj=1,nstn
      call readdata (fpath,jj,idx(jj),nyrs(jj),nstrt(jj),avrf(jj))
      enddo

c	if data at target station is available (validation), read the data
c      if(ival.eq.1)call readdatah(fpath,rfh,filet,nday,nyh,nst,
c     1             nout,tavrf)
      nsh=nst-1
      w(nstn+1)=1.0
      i=0
c	write(*,21)
      write(21,21)
c	write(*,*)
      write(21,*)
 21   format('    No Index Weight Years St_year Av annual rainfall')
c	write(*,*)'target Station'
      write(21,*)'target Station'
c	write(*,200)I,idxt,w(nstn+1),nyh,nsh,tanrf
      write(21,200)I,idxt,w(nstn+1),nyh,nsh,tanrf
c	write(*,*)
      write(21,*)
c	write(*,*)'Nearby Stations'
      write(21,*)'Nearby Stations'
      DO I=1,nstn
      write(21,200)I,idx(i),w(i),nyrs(i),nstrt(i),avrf(i)
c	write(*,200)I,idx(i),w(i),nyrs(i),nstrt(i),avrf(i)
      if(i.gt.1)w(i)=w(i)+w(i-1)
      enddo
      w(nstn)=1.0
      close(unit=21)
      
      if(idrop.ne.0)return

c	write(*,*)' Do you want to drop any nearby station?'
c	write(*,*)' If yes, enter station sr no otherwise enter zero'
c	read(*,*)idrop
c	if(idrop.gt.0)then
c	kk=kk+1
c	idr(kk)=idx(idrop)
c	goto 17
c	endif

 200  format(2i6,f7.3,2i6,2x,f8.2,5x,f8.2)
 300  format(i6,a6,f7.3,2i6,2f8.2)

C      file_rf_out='MMM_'
C      ii=4
C      write(cc,'(i6.6)')idxt
C      do i=1,6
C      if(cc(i:i).ne.' ')then
C      ii=ii+1
C      file_rf_out(ii:ii)=cc(i:i)
C      endif
C      enddo
C      file_rf_out(ii+1:ii+1)='.'
C      file_rf_out(ii+2:ii+2)='o'
C      file_rf_out(ii+3:ii+3)='u'
C      file_rf_out(ii+4:ii+4)='t'



c	form a moving window for each day and calculate probabilities, mean sd and serial corrl
      call smoothprob()
      if(nlon.gt.0)then
      call long_store()
      call av_sd_lon()
      call smoothavcov()
      endif
      if(iamt.gt.0)call rf_amt_store()
      open(2,file=filet)
      if(ival.gt.0)call stat_h(rfh,1,nsh,nyh,rain,nout,nday,isn)
      call simulate (w,avrf,tanrf)
      close(unit=2)
      call result(isn,1,nsim,nout,ng,ng)
c	isuccess=1
      return 
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine day(nday,monmax)
c      Implicit double precision (A-H,O-Z) 
      integer nday(monmax)
      nday(1)=31
      nday(2)=29
      nday(3)=31
      nday(4)=30
      nday(5)=31
      nday(6)=30
      nday(7)=31
      nday(8)=31
      nday(9)=30
      nday(10)=31
      nday(11)=30
      nday(12)=31
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine daycount(ns,i,id)
c Implicit double precision (A-H,O-Z) 
      integer ns,i,id 
      id=28
      if (mod(ns+i,400).eq.0)id=29
      if (mod(ns+i,100).ne.0.and.mod(ns+i,4).eq.0)id=29
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine day_neg(i,j,l,ip,nout,ns,lag,nday,monmax,ind)
c check the days of moving window for negative values of day
c Implicit double precision (A-H,O-Z)
      integer nday(monmax)
      ind=0
      if(i.le.ip.and.j.eq.1.and.i.le.ip.and.l-lag.lt.1)then
      ind=1
      return
      endif
 10   if(l.gt.0)return
      j=j-1
      if(j.lt.1)then
      j=nout
      i=i-1
      if(i.lt.ip)then
      ind=1
      return
      endif
      endif
       nd=nday(j)
      if(j.eq.2)call daycount(ns,i,nd)
      l=l+nd
      goto 10
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine day_pos(i,j,l,nout,ns,nday,ny,monmax,ind)
c	check the day for moving window for positive values of day
c      Implicit double precision (A-H,O-Z)
      integer nday(monmax)
      ind=0
      if(i.eq.ny.and.j.eq.nout.and.i.eq.ny.and.l.gt.nday(nout))then
      ind=1
      return
      endif
 10   nd=nday(j)
      if(j.eq.2)call daycount(ns,i,nd)
       if(l.le.nd)then
      return
      endif
      l=l-nd
      j=j+1
      if(j.gt.nout)then
      i=i+1
      j=1
      if(i.gt.ny)then
      ind=1
      return
      endif
      endif
      goto 10
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine to identify the nearest stations for parameter estimation at given ungauged location
c	read nearby station coordinates
      subroutine station(kk,idr,tlat,tlon,tele,tdcoast,ttemp,
     1               tanrf,idx,w,nstn,nsmax,nstnmax)
      real*4 lat(nsmax),lon(nsmax),ele(nsmax),dcoast(nsmax),
     1       temp(nsmax),avrf(nsmax),w(nstnmax),
     2       aa(50),air(nsmax),dlatlon(nsmax)
      integer indx(nsmax),ir(nsmax),ilat(nsmax),ilon(nsmax),
     1     iele(nsmax),idcoast(nsmax),iavrf(nsmax),idx(nstnmax),
     2     idlatlon(nsmax),idr(nsmax),itemp(nsmax)
       
      open(1,file='stn_record.dat',status='old')
      read(1,*)
      read(1,*)
      ii=0
      do i=1,50000
      read(1,10,end=11)i1,iy,id,tt,alt,aln,aele,adcoast,atemp,
     1    (aa(j),j=1,16),anrf
 10	  format(2i5,i10,a50,2f10.4,20f10.2)

      if(kk.gt.0)then
      do k1=1,kk
      if(idr(k1).eq.id)goto 12
      enddo
      endif
c	exclude the target location and the station to be dropped
      if(alt.ne.tlat.and.aln.ne.tlon)then
      ii=ii+1
      ele(ii)=abs(aele-tele)
      lat(ii)=abs(alt-tlat)
      lon(ii)=abs(aln-tlon)
      avrf(ii)=abs(anrf-tanrf)
      temp(ii)=abs(atemp-ttemp)
      dcoast(ii)=abs(tdcoast-adcoast)*2/(tdcoast+adcoast)
      dlatlon(ii)=lat(ii)*lon(ii)
      indx(ii)=id
      endif
 12	  continue
      enddo
 11	  continue

c	rank nearest station for each physiographic attribute
      call ranks(ele,iele,ii,nsmax)
      call ranks(lat,ilat,ii,nsmax)
      call ranks(lon,ilon,ii,nsmax)
      call ranks(avrf,iavrf,ii,nsmax)
      call ranks(dcoast,idcoast,ii,nsmax)
      call ranks(dlatlon,idlatlon,ii,nsmax)
      call ranks(temp,itemp,ii,nsmax)
      air=0.0
c combine nearest neighbours under each attribute
c	lat 1.0, lon 0.25, ele 1.0, avrf 0.5, dcoast 0.5
      do i=1,ii
      air(ilat(i))=air(ilat(i))+float(i)
      air(ilon(i))=air(ilon(i))+float(i)*0.5
      air(iele(i))=air(iele(i))+float(i)*0.1
      air(iavrf(i))=air(iavrf(i))+float(i)*0.5
      air(idcoast(i))=air(idcoast(i))+float(i)*0.1
      air(idlatlon(i))=air(idlatlon(i))+float(i)
	  air(itemp(i))=air(itemp(i))+float(i)*0.25
      enddo
      
      call ranks(air,ir,ii,nsmax)

c	calculate stations weights on the basis of total ranking
      sum=0.0
      do i=1,nstn
      w(i)=1.0/air(i)
      sum=sum+w(i)
      enddo
      do i=1,nstn
      w(i)=w(i)/sum
      idx(i)=indx(ir(i))
      enddo
c	write(*,*)(idx(j),j=1,nstn)
c	write(*,*)(w(j),j=1,nstn)
c	pause

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ranks(x,ix,n,nmax)
      real*4 x(nmax)
      integer ix(nmax)
      do i=1,n
      ix(i)=i
      enddo
      do i=1,n-1
      do j=i+1,n
      if(x(j).lt.x(i))then
      a=x(i)
      x(i)=x(j)
      x(j)=a
      ia=ix(i)
      ix(i)=ix(j)
      ix(j)=ia
      endif
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine readdata (fpath,jj,idx,nyear,nstart,avr)
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include 'data.inc'
      include 'occrcor.inc'
      include 'p30avsd.inc'
c	subroutine to read the daily rainfall data from a file
      character*80 ff*5,hh,fpath,tt*16
      
      write(ff,'(i5.5)')idx
      hh=' '
      ii=0
      do i=1,80
      if(fpath(i:i).ne.' ')then
      ii=ii+1
      hh(ii:ii)=fpath(i:i)
      endif
      enddo
      tt='rev_dr046115.txt'
      
      do i=1,5
      tt(7+i:7+i)=ff(i:i)
      enddo
      
      do i=1,16
      hh(ii+i:ii+i)=tt(i:i)
      enddo


c	read daily rainfall

      open(1,file=hh,status='old')
      read(1,100)head
      
      i1=0
      ind=0
      do i=1,1000000000
      read(1,*,err=107,end=302)iy,i2,i3,bb
      if(i2.eq.1.and.i3.eq.1.and.i2.eq.1.and.ind.eq.0)then
      nstart=iy
      ind=1
      ny=nstart-1
      endif
      if(iy.ge.nstart)then
c	check for assigning the sequence of a year
      if(iy.gt.ny)then
c	write(*,*)'....... year ',ny, i1
      i1=i1+1
      ny=iy
      endif
      rf(jj,i1,i2,i3)=bb
      endif
 399  continue
      enddo
 302  close(unit=1)
      if(i2.ne.nout.and.i3.ne.31)i1=i1-1
      nyear=i1

c	calculate average annual rainfall
      avr=0.0
      do i=1,nyear
      do j=1,nout
      nd=nday(j)
      if(j.eq.2)call daycount(nstart-1,i,nd)
      do l=1,nd
      avr=avr+rf(jj,i,j,l)
      enddo
      enddo
      enddo
      avr=avr/float(nyear)
      
      return
 100  format(a)
 107  write(*,*)'Error in data file'
      write(*,*)'year month day ',iy,i2,i3
c	stop
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smoothprob()
c find probability of a wet day given rf of previous day
c Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      include  'occr.inc'
      include 'rf_amt.inc'
      real*4 xc(nearmax,nstnmax)
           real*4 prow(5,nstnmax),pr(nstnmax)
           real*4 xp1(nearmax,nstnmax)
           real*4 xn(nearmax,nstnmax)
      integer i,j,l,kk(nstnmax),nk
      integer ip,jp,lp
       
       
      write(*,*)'Calculating probabilities ..... '
       
       
      do j=1,nout
      is=iseas(j,isn)
      do l=1,nday(j)
      
      do nk=1,nstn
      do i=1,nearmax
      xc(i,nk)=0.0
      xp1(i,nk)=0.0
      enddo
      
      kk(nk)=0
      nyear=nyrs(nk)
      ns=nstrt(nk)-1
      
      do 11 i=1,nyear,1
      if(j.eq.2)then
      call daycount(ns,i,lc)
      if(l.gt.lc)goto 11
      endif
      
      li=l-iband-1
      
      do 10 jh=1,iband*2+1
      li=li+1
      ic=i
      jc=j
      lc=li
      if(li.le.2)call day_neg(ic,jc,lc,1,nout,ns,2,nday,monmax,indx)
      if(li.gt.2)call day_pos(ic,jc,lc,nout,ns,nday,nyear,monmax,indx)
      if(indx.eq.1.)goto 10
c check for prev day
      lp=lc-1
      jp=jc
      ip=ic
      if(lp.lt.1)then
      if(ip.le.1.and.jp.eq.1)goto 10
      jp=jp-1
      if(jp.lt.1)then
      ip=ip-1
      jp=nout
      endif
      lp=nday(jp)
      if(jp.eq.2)call daycount(ns,ip,lp)
      endif
c check for the next day

      in=ic
      jn=jc
      ln=lc+1
      nd=nday(jn)
      if(jn.eq.2)call daycount(ns,in,nd)
c check for windows of last few days of the last year
      if(in.eq.nyear.and.jn.eq.nout.and.jn.eq.nout.and.ln.gt.nd)goto 10
      if(ln.gt.nd)then
      jn=jn+1
      ln=1
      if(jn.gt.nout)then
      in=in+1
      jn=1
      endif
      nd=nday(jn)
      if(in.gt.nyear)goto 10
      endif
      
      kk(nk)=kk(nk)+1
      xc(kk(nk),nk)=rf(nk,ic,jc,lc)
      xp1(kk(nk),nk)=rf(nk,ip,jp,lp)
      xn(kk(nk),nk)=rf(nk,in,jn,ln)
	    
 10	  continue
 11	  continue
      enddo

c calculate conditional probabilies
      call prob1(xc,xp1,nstn,kk,rain,prow,nearmax,nstnmax)
      do nk=1,nstn
      pro1(j,l,nk)=prow(1,nk)
      pro2(j,l,nk)=prow(2,nk)
      pro(j,l,nk)=prow(3,nk)
      pr(nk)=prow(3,nk)
      write(13,60)j,l,nk,(prow(k1,nk),k1=1,3)
      enddo
      enddo
      enddo
 60	  format(3i5,5f8.3)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prob1(xc,xp1,nstn,kp,rain,prow,nearmax,nstnmax)
c      Implicit double precision (A-H,O-Z) 
        real*4 xc(nearmax,nstnmax),xp1(nearmax,nstnmax)
        real*4 prow(5,nstnmax),sum1,sum2,sum5,rain
      integer i,kp(nstnmax),nk
      do nk=1,nstn
      prow(1,nk)=0.0
      prow(2,nk)=0.0
      prow(3,nk)=0.0
      prow(4,nk)=0.0
      prow(5,nk)=0.0
      sum1=0.0
      sum2=0.0
      sum5=0.0
      k1=0
      k2=0
      
      do i=1,kp(nk)
      if(xc(i,nk).ge.rain.and.xp1(i,nk).lt.rain)sum1=sum1+1.0
      if(xc(i,nk).ge.rain.and.xp1(i,nk).ge.rain)sum2=sum2+1.0
      if(xc(i,nk).ge.rain)sum5=sum5+1.0
       if(xp1(i,nk).lt.rain)k1=k1+1
      if(xp1(i,nk).ge.rain)k2=k2+1
      enddo
      if(k1.gt.0)prow(1,nk)=sum1/float(k1)
      if(k2.gt.0)prow(2,nk)=sum2/float(k2)
      if(kp(nk).gt.0)prow(3,nk)=sum5/float(kp(nk))
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine long_store()
c	subroutine to calculate and store higher time scale variables
      include  'para.inc'
      include 'data.inc'
      include 'occr.inc'
      real*4 prev(nvrmax,nstnmax),dsum
      
      do nk=1,nstn
      nyear=nyrs(nk)
      ns=nstrt(nk)-1
      do i=1,nyear
      do j=1,nout
      lc=nday(j)
      if(j.eq.2)call daycount(ns,i,lc)
      do l=1,lc
      do jj=1,nlon
      jk=(nk-1)*nlon+jj
      high(jk,i,j,l)=-9999.0
      enddo
      enddo
      enddo
      enddo
      enddo
      
      do nk=1,nstn
      nyear=nyrs(nk)
      ns=nstrt(nk)-1
      do i=iyrst,nyear,1
      do j=1,nout
      lc=nday(j)
      if(j.eq.2)call daycount(ns,i,lc)
      do l=1,lc
      do jj=1,nlon
      call prevdays(i,j,l,lon(1,jj),lon(2,jj),dsum,ind,nk,ns)
      if(ind.eq.1)then
      write(*,*)'check higher time scale routine'
      goto 110
      endif
      jk=(nk-1)*nlon+jj
      high(jk,i,j,l)=dsum
      enddo
 110  continue

      enddo
      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prevdays(i1,j1,l1,long1,long2,sum,ind,nk,ns)
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      real*4 sum  
      
      ind=0
      sum=0.0
      lp=l1
      jp=j1
      ip=i1
      do l=1,long2
      lp=lp-1
      if(lp.lt.1)then
      if(ip.eq.1.and.jp.eq.1)goto 10
      jp=jp-1
      if(jp.lt.1)then
      ip=ip-1
      if(ip.lt.1)goto 10
      jp=nout
      endif
      lp=nday(jp)
           if(jp.eq.2)call daycount(ns,ip,lp)
      endif
      if(l.ge.long1)then
      if(rf(nk,ip,jp,lp).ge.rain)sum=sum+1.0
      endif
       enddo
      goto 20
 10	  ind=1
      return
 20   al=float(long2-long1)+1.0
      return
      write(*,30)i1,j1,l1,sum
 30   format(3i4,20f6.1)
	
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	Calculate mean and SD of higher time scale data
      subroutine av_sd_lon()
      include  'para.inc'
      include 'data.inc'
      include 'occr.inc'
      real*4 x(nearmax),prev(nvrmax,nstnmax),
     1         x1(nearmax,nstnmax,nvmax),asd(nstnmax,nstnmax),
     2         pp(nvrmax,nstnmax)
      integer nn(nstnmax)
      
      do j=1,nout
      do l=1,nday(j)
      
      do nk=1,nstn
      nn(nk)=0
      nyear=nyrs(nk)
      ns=nstrt(nk)
      do 11 i=iyrst,nyear,1
      if(j.eq.2)then
      call daycount(ns,i,lc)
      if(l.gt.lc)goto 11
      endif
      
      li=l-iband-1
      
      do 10 jh=1,iband*2+1
      li=li+1
      ic=i
      jc=j
      lc=li
      if(li.le.2)call day_neg(ic,jc,lc,1,nout,ns,2,nday,monmax,indx)
      if(li.gt.2)call day_pos(ic,jc,lc,nout,ns,nday,nyear,monmax,indx)
      if(indx.eq.1.)goto 10
c	check for prev day
      lp=lc-1
      jp=jc
      ip=ic
      if(lp.lt.1)then
      if(ip.le.1.and.jp.eq.1)goto 10
      jp=jp-1
      if(jp.lt.1)then
      ip=ip-1
      jp=nout
      endif
      lp=nday(jp)
      if(jp.eq.2)call daycount(ns,ip,lp)
      endif
      
      do jj=1,nlon
      jk=(nk-1)*nlon+jj
      if(high(jk,ic,jc,lc).lt.0)goto 10
      prev(jj,nk)=high(jk,ic,jc,lc)
      enddo
      nn(nk)=nn(nk)+1
      do jj=1,nlon
      x1(nn(nk),nk,jj)=prev(jj,nk)
      enddo

 10   continue
 11   continue
      enddo

      do nk=1,nstn
      do jj=1,nlon
      do i=1,nn(nk)
      x(i)=x1(i,nk,jj)
      enddo
      call basic(x,avl(jj,j,l,nk),sdl(jj,j,l,nk),nn(nk))
      if(j.eq.2.and.l.eq.29)avl(jj,j,l,nk)=avl(jj,j,l-1,nk)
      if(j.eq.2.and.l.eq.29)sdl(jj,j,l,nk)=sdl(jj,j,l-1,nk)
      if(nk.eq.1)write(13,101)nk,jj,j,l,avl(jj,j,l,nk),sdl(jj,j,l,nk)
      enddo
      enddo
      
      enddo
      enddo

 101  format(4i5,7f10.2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE BASIC (data,AVE,SD,N)
c      Implicit double precision (A-H,O-Z) 
      INTEGER n
      REAL ave,var,sd,data(n)
      INTEGER j
      REAL s,ep
      ave=0.0
	  sd=1.0
	  if(n.le.1)return
      do 11 j=1,n
      ave=ave+data(j)
 11   continue
      if(n.gt.0)ave=ave/n
      var=0.0
      ep=0.0
      do 12 j=1,n
      s=data(j)-ave
      ep=ep+s
      var=var+s*s
 12   continue
      sd=sqrt((var-ep**2/n)/(n-1))
      return
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prevavsd0(xc,xp3,nstn,nv,kp,rain,
     1                    ppav,ppcov,pav,pcov,iflag,
     2                    nearmax,nstnmax,nvrmax)

c      Implicit double precision (A-H,O-Z) 

      real*4 xc(nearmax,nstnmax)
      real*4 xp3(nvrmax,nearmax,nstnmax)
      real*4 ppav(2,nvrmax,nstnmax)
      real*4 pav(2,nvrmax,nstnmax)
      real*4 sum(nvrmax)
      real*4 xbar(2,nvrmax)
      real*4 xss(2,nvrmax,nvrmax)
      real*4 ppcov(2,nvrmax,nvrmax,nstnmax),
     1       pcov(2,nvrmax,nvrmax,nstnmax)
      real*4 sighat(nvrmax,nvrmax),df,df1,df2
      real*4 sighat1(nvrmax,nvrmax),sighat2(nvrmax,nvrmax)
      integer nnn(2)
      integer i,kp(nstnmax),nk,nv

      do nk=1,nstn
      do l=1,2
      do i=1,nv+1
      do j=1,nv+1
      pcov(l,i,j,nk)=0.0
      ppcov(l,i,j,nk)=0.0
      enddo
      pav(l,i,nk)=0.0
      ppav(l,i,nk)=0.0
      enddo
      enddo
      xbar=0.0
      nnn=0
      xss=0.0
      det=0.0
      det1=0.0
      det2=0.0
      do i=1,kp(nk)
      ixc=1
      if(xc(i,nk).ge.rain)ixc=2
      nnn(ixc) = nnn(ixc)+1
      
      do kk=1,nv
      sum(kk)=xp3(kk,i,nk)
      enddo


      do 62 k=1,nv
      xbar(ixc,k) = xbar(ixc,k) + sum(k)
      do 62 ll=1,nv
 62   xss(ixc,k,ll) = xss(ixc,k,ll) + sum(k)*sum(ll)
 60   continue
      enddo

      do 70 i=1,2
      do 70 k=1,nv
      if (nnn(i).gt.0) xbar(i,k) = xbar(i,k)/nnn(i)
c	write(*,*) k,i,xbar(i,k),xbar(i,k),nnn(i)
 70   continue
c	
      df = kp(nk)-2
      df1 = nnn(1)-1
      df2 = nnn(2)-1

      do 71 k=1,nv
      do 71 l=1,nv
      sighat(k,l)=0.0
      sighat1(k,l)=0.0
      sighat2(k,l)=0.0
      do 72 i=1,2
      if(i.eq.1)sighat1(k,l) = sighat1(k,l) +
     1      (xss(i,k,l) - nnn(i)*xbar(i,k)*xbar(i,l))
      if(i.eq.2)sighat2(k,l) = sighat2(k,l) +
     1      (xss(i,k,l) - nnn(i)*xbar(i,k)*xbar(i,l))
 72   sighat(k,l) = sighat(k,l) +
     1      (xss(i,k,l) - nnn(i)*xbar(i,k)*xbar(i,l))
      if(df1.gt.0.0)sighat1(k,l) = sighat1(k,l)/df1
      if(df2.gt.0.0)sighat2(k,l) = sighat2(k,l)/df1
 71    if(df.gt.0.0)sighat(k,l) = sighat(k,l)/df

      do 73 k=1,nv
      do 73 l=1,nv
      if(df1.lt.20)sighat1(k,l)=sighat(k,l)
 73   if(df2.lt.20)sighat2(k,l)=sighat(k,l)

c	call subroutine to compute inverse of covar matrix
c	do i=1,nv
c	write(*,*)(sighat1(i,j),j=1,nv),nv
c	write(*,*)(sighat2(i,j),j=1,nv)
c	write(*,*)(sighat(i,j),j=1,nv)
c	enddo
      call solve(sighat,nv,nvrmax,det)
      if(iflag.eq.0)call solve(sighat1,nv,nvrmax,det1)
      if(iflag.eq.0)call solve(sighat2,nv,nvrmax,det2)
c	check for no of observations
      if(nnn(1).lt.10)det1=0.0
      if(nnn(2).lt.10)det2=0.0
      if(kp(nk).lt.10)det=0.0
      do j1=1,nv
      ppav(1,j1,nk)=xbar(2,j1)
      pav(1,j1,nk)=xbar(1,j1)
      enddo
      if(iflag.eq.0)then
      pav(1,nv+1,nk)=det1
      ppav(1,nv+1,nk)=det2
      else
      ppav(1,nv+1,nk)=det
      pav(1,nv+1,nk)=det
      endif
      do j=1,nv
      do j1=1,nv
      if(iflag.eq.0)then
      pcov(1,j,j1,nk)=(sighat1(j,j1))
      ppcov(1,j,j1,nk)=(sighat2(j,j1))
      else
      pcov(1,j,j1,nk)=(sighat(j,j1))
      ppcov(1,j,j1,nk)=(sighat(j,j1))
      endif
      enddo
      enddo
      enddo
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prevavsd1(xc,xp1,xp3,nstn,nv,kp,rain,
     1                    ppav,ppcov,pav,pcov,iflag,
     2                    nearmax,nstnmax,nvrmax)

c      Implicit double precision (A-H,O-Z) 

      real*4 xc(nearmax,nstnmax),xp1(nearmax,nstnmax)
      real*4 xp3(nvrmax,nearmax,nstnmax)
      real*4 ppav(2,nvrmax,nstnmax)
      real*4 pav(2,nvrmax,nstnmax),sum(nvrmax)
      real*4 xbar(2,2,nvrmax)
      real*4 xss(2,2,nvrmax,nvrmax)
      real*4 ppcov(2,nvrmax,nvrmax,nstnmax),
     1       pcov(2,nvrmax,nvrmax,nstnmax)
      real*4 sighat(nvrmax,nvrmax),df,df1,df2
      real*4 sighat1(nvrmax,nvrmax),sighat2(nvrmax,nvrmax)
      integer nnn(2,2)
      integer i,kp(nstnmax),nk,nv
      
      do nk=1,nstn
      do l=1,2
      do i=1,nv+1
      do j=1,nv+1
      pcov(l,i,j,nk)=0.0
      ppcov(l,i,j,nk)=0.0
      enddo
      pav(l,i,nk)=0.0
      ppav(l,i,nk)=0.0
      enddo
      enddo
      xbar=0.0
      nnn=0
      xss=0.0
      det=0.0
      det1=0.0
      det2=0.0
      do i=1,kp(nk)
      ixc=1
      if(xc(i,nk).ge.rain)ixc=2
      ixp=1
      if(xp1(i,nk).ge.rain)ixp=2
      
      nnn(ixc,ixp) = nnn(ixc,ixp)+1
      
      do kk=1,nv
      sum(kk)=xp3(kk,i,nk)
      enddo

      do 62 k=1,nv
      xbar(ixc,ixp,k) = xbar(ixc,ixp,k) + sum(k)
      do 62 ll=1,nv
 62   xss(ixc,ixp,k,ll) = xss(ixc,ixp,k,ll) + sum(k)*sum(ll)
 60   continue
      enddo

      do 70 i=1,2
      do 70 j=1,2 
      do 70 k=1,nv
      if (nnn(i,j).gt.0) xbar(i,j,k) = xbar(i,j,k)/nnn(i,j)
c	write(*,*) k,i,j,xbar(i,j,k),xbar(i,j,k),nnn(i,j)
 70   continue
c	
      df = kp(nk)-2*2
      df1=nnn(1,1)+nnn(2,1)-2
      df2=nnn(1,2)+nnn(2,2)-2

      do 71 k=1,nv
      do 71 l=1,nv
      sighat(k,l)=0.0
      sighat1(k,l)=0.0
      sighat2(k,l)=0.0
      do 72 i=1,2
      do 72 j=1,2 
      if(j.eq.1)sighat1(k,l) = sighat1(k,l) +
     1      (xss(i,j,k,l) - nnn(i,j)*xbar(i,j,k)*xbar(i,j,l))
      if(j.eq.2)sighat2(k,l) = sighat2(k,l) +
     1      (xss(i,j,k,l) - nnn(i,j)*xbar(i,j,k)*xbar(i,j,l))
 72   sighat(k,l) = sighat(k,l) +
     1      (xss(i,j,k,l) - nnn(i,j)*xbar(i,j,k)*xbar(i,j,l))
      if(df1.gt.0.0)sighat1(k,l) = sighat1(k,l)/df1
      if(df2.gt.0.0)sighat2(k,l) = sighat2(k,l)/df2
 71    if(df.gt.0.0)sighat(k,l) = sighat(k,l)/df

      do 73 k=1,nv
      do 73 l=1,nv
      if(df1.lt.20)sighat1(k,l)=sighat(k,l)
 73   if(df2.lt.20)sighat2(k,l)=sighat(k,l)
c	if(df1.lt.20)pro1(jk,lk,nk)=pro(jk,lk,nk)
c	if(df2.lt.20)pro2(jk,lk,nk)=pro(jk,lk,nk)

c	call subroutine to compute inverse of covar matrix
      call solve(sighat,nv,nvrmax,det)
      if(iflag.eq.0)call solve(sighat1,nv,nvrmax,det1)
      if(iflag.eq.0)call solve(sighat2,nv,nvrmax,det2)

c	do i=1,nv
c	write(*,*)(sighat1(i,j),j=1,nv),nv
c	write(*,*)(sighat2(i,j),j=1,nv)
c	write(*,*)(sighat(i,j),j=1,nv)
c	enddo
c	do i=1,nv
c	write(*,*)xbar(1,1,i),xbar(1,2,i),xbar(2,1,i),xbar(2,2,i)
c	enddo
c	

c	check for no of observations
      if(nnn(2,1).lt.10)det1=0.0
      if(nnn(1,1).lt.10)det1=0.0
      if(nnn(2,2).lt.10)det2=0.0
      if(nnn(1,2).lt.10)det2=0.0
      if(kp(nk).lt.10)det=0.0
      do j1=1,nv
      ppav(1,j1,nk)=xbar(2,1,j1)
      pav(1,j1,nk)=xbar(1,1,j1)
      ppav(2,j1,nk)=xbar(2,2,j1)
      pav(2,j1,nk)=xbar(1,2,j1)
      enddo
      if(iflag.eq.0)then
      pav(1,nv+1,nk)=det1
      pav(2,nv+1,nk)=det2
      ppav(1,nv+1,nk)=det1
      ppav(2,nv+1,nk)=det2
      else
      ppav(1,nv+1,nk)=det
      ppav(2,nv+1,nk)=det
      pav(1,nv+1,nk)=det
      pav(2,nv+1,nk)=det
      endif
      do j=1,nv
      do j1=1,nv
      if(iflag.eq.0)then
      pcov(1,j,j1,nk)=(sighat1(j,j1))
      pcov(2,j,j1,nk)=(sighat2(j,j1))
      ppcov(1,j,j1,nk)=(sighat1(j,j1))
      ppcov(2,j,j1,nk)=(sighat2(j,j1))
      else
      pcov(1,j,j1,nk)=(sighat(j,j1))
      pcov(2,j,j1,nk)=(sighat(j,j1))
      ppcov(1,j,j1,nk)=(sighat(j,j1))
      ppcov(2,j,j1,nk)=(sighat(j,j1))
      endif
      enddo
      enddo
      enddo
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine solve(ss,nv,nvrmax,det)
c      Implicit double precision (A-H,O-Z) 
      real*4 det
      real*4 ss(nvrmax,nvrmax)        
      real*4 A(nvrmax,nvrmax),W(nvrmax),V(nvrmax,nvrmax),
     &        temp(nvrmax,nvrmax)


      tol = 1.E-7
      if(nv.eq.1)then
      det=ss(1,1)
      if(ss(1,1).ne.0.0e0)ss(1,1)=1.0e0/ss(1,1)
      return
      endif		 
      do 100 i=1,nv
      do 110 j=1,nv
      A(i,j)=ss(i,j)
 110  continue
 100  continue
      call SVDCMP(a,nv,nv,nvrmax,nvrmax,w,v)

c this bit is calculation of 1/wj*U^T - see eqn 2.9.5 in NR
	  temp=0.0
      do 200 i=1,nv
      do 210 j=1,nv
      if(w(i).gt.tol)temp(i,j)=A(j,i)/w(i)
 210  continue
 200  continue

c this bit is calculation of 2.9.5 in NR
      call matmul(v,temp,ss,nv,nvrmax)

c	calculate determinant
      det=1.0
      do 220 i=1,nv
      det=det*w(i)
 220  continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matmul(A,B,C,n,nmax)
c      Implicit double precision (A-H,O-Z) 
      real*4 A(nmax,nmax), B(nmax,nmax), C(nmax,nmax)

      do 100 i=1,n
      do 110 j=1,n
 110  C(i,j)=0.0
 100  continue
      do 200 i=1,n
      do 210 j=1,n
      do 220 k=1,n
      C(i,j)=C(i,j)+A(i,k)*B(k,j)
 220  continue
 210  continue
 200  continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c subroutine svdcmp
c Extracted from: Numerical Recipes in Fortran 77
c purpose: to carry out singular value decomposition of a matrix
c       given a matrix A(m,n), this routine computes its singular value
c       decomposition, A=U*W*V(T). The matrix U replaces A on output.
c       The diagonal matrix of singular values W is output as a vector
c       W(n). The matrix V is output as V(n,n).
c Input: a(m,n) - matrix A
c        w(n) - dummy variable.
c        v(n,n) - dummy variable.
c        m - actual dimension of matrix
c        n - actual dimension of matrix
c        mp - max. dimension of matrix
c        np - max. dimension of matrix
c Output: a(m,n) - the U matrix.
c         w(n) - vector of diagonal matrix of singular values.
c         v(n,n) - the V matrix.
c
c Calling subroutines: nil.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      SUBROUTINE svdcmp(a,m,n,mp,np,w,v)
c      Implicit double precision (A-H,O-Z)
      INTEGER m,mp,n,np,NMAX
      dimension a(mp,np),v(np,np),w(np)
      PARAMETER (NMAX=100)
      dimension rv1(NMAX)
C    DOES NOT USE pythag - modification by AS
      INTEGER i,its,j,jj,k,l,nm
      real*4 anorm,c,f,g,h,s,scale,x,y,z
      g=0.0
      scale=0.0
      anorm=0.0
      do 25 i=1,n
      l=i+1
      rv1(i)=scale*g
      g=0.0
      s=0.0
      scale=0.0
      if(i.le.m)then
      do 11 k=i,m
      scale=scale+abs(a(k,i))
 11   continue
      if(scale.ne.0.0)then
      do 12 k=i,m
      a(k,i)=a(k,i)/scale
      s=s+a(k,i)*a(k,i)
 12   continue
      f=a(i,i)
      g=-sign(sqrt(s),f)
      h=f*g-s
      a(i,i)=f-g
      do 15 j=l,n
      s=0.0
      do 13 k=i,m
      s=s+a(k,i)*a(k,j)
 13   continue
      f=s/h
      do 14 k=i,m
      a(k,j)=a(k,j)+f*a(k,i)
 14   continue
 15   continue
      do 16 k=i,m
      a(k,i)=scale*a(k,i)
 16   continue
      endif
      endif
      w(i)=scale *g
      g=0.0
      s=0.0
      scale=0.0
      if((i.le.m).and.(i.ne.n))then
      do 17 k=l,n
      scale=scale+abs(a(i,k))
 17   continue
      if(scale.ne.0.0)then
      do 18 k=l,n
      a(i,k)=a(i,k)/scale
      s=s+a(i,k)*a(i,k)
 18   continue
      f=a(i,l)
      g=-sign(sqrt(s),f)
      h=f*g-s
      a(i,l)=f-g
      do 19 k=l,n
      rv1(k)=a(i,k)/h
 19   continue
      do 23 j=l,m
      s=0.0
      do 21 k=l,n
      s=s+a(j,k)*a(i,k)
 21   continue
      do 22 k=l,n
      a(j,k)=a(j,k)+s*rv1(k)
 22   continue
 23   continue
      do 24 k=l,n
      a(i,k)=scale*a(i,k)
 24   continue
      endif
      endif
      anorm1=abs(w(i))+abs(rv1(i))
      if(anorm1.gt.anorm)anorm=anorm1
 25   continue
      do 32 i=n,1,-1
      if(i.lt.n)then
      if(g.ne.0.0)then
      do 26 j=l,n
      v(j,i)=(a(i,j)/a(i,l))/g
 26   continue
      do 29 j=l,n
      s=0.0
      do 27 k=l,n
      s=s+a(i,k)*v(k,j)
 27   continue
      do 28 k=l,n
      v(k,j)=v(k,j)+s*v(k,i)
 28   continue
 29   continue
      endif
      do 31 j=l,n
      v(i,j)=0.0
      v(j,i)=0.0
 31   continue
      endif
      v(i,i)=1.0
      g=rv1(i)
      l=i
 32   continue
      do 39 i=min(m,n),1,-1
      l=i+1
      g=w(i)
      do 33 j=l,n
      a(i,j)=0.0
 33   continue
      if(g.ne.0.0)then
      g=1.0/g
      do 36 j=l,n
      s=0.0
      do 34 k=l,m
      s=s+a(k,i)*a(k,j)
 34   continue
      f=(s/a(i,i))*g
      do 35 k=i,m
      a(k,j)=a(k,j)+f*a(k,i)
 35   continue
 36   continue
      do 37 j=i,m
      a(j,i)=a(j,i)*g
 37   continue
      else
      do 38 j= i,m
      a(j,i)=0.0
 38   continue
      endif
      a(i,i)=a(i,i)+1.0
 39   continue
      do 49 k=n,1,-1
      do 48 its=1,30
      do 41 l=k,1,-1
      nm=l-1
      if((abs(rv1(l))+anorm).eq.anorm)  goto 2
      if((abs(w(nm))+anorm).eq.anorm)  goto 1
 41   continue
 1    c=0.0
      s=1.0
      do 43 i=l,k
      f=s*rv1(i)
      rv1(i)=c*rv1(i)
      if((abs(f)+anorm).eq.anorm) goto 2
      g=w(i)
      h=sqrt(f*f+g*g)
      w(i)=h
      h=1.0/h
      c= (g*h)
      s=-(f*h)
      do 42 j=1,m
      y=a(j,nm)
      z=a(j,i)
      a(j,nm)=(y*c)+(z*s)
      a(j,i)=-(y*s)+(z*c)
 42   continue
 43   continue
 2    z=w(k)
      if(l.eq.k)then
      if(z.lt.0.0)then
      w(k)=-z
      do 44 j=1,n
      v(j,k)=-v(j,k)
 44   continue
      endif
      goto 3
      endif
c     if(its.eq.30)  'no convergence in svdcmp'
      x=w(l)
      nm=k-1
      y=w(nm)
      g=rv1(nm)
      h=rv1(k)
      f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
      g=sqrt(f*f+1.0)
      f=((x-z)*(x+z)+h*((y/(f+sign(g,f)))-h))/x
      c=1.0
      s=1.0
      do 47 j=l,nm
      i=j+1
      g=rv1(i)
      y=w(i)
      h=s*g
      g=c*g
      z=sqrt(f*f+h*h)
      rv1(j)=z
      c=f/z
      s=h/z
      f= (x*c)+(g*s)
      g=-(x*s)+(g*c)
      h=y*s
      y=y*c
      do 45 jj=1,n
      x=v(jj,j)
      z=v(jj,i)
      v(jj,j)= (x*c)+(z*s)
      v(jj,i)=-(x*s)+(z*c)
 45   continue
      z=sqrt(f*f+h*h)
      w(j)=z
      if(z.ne.0.0)then
      z=1.0/z
      c=f*z
      s=h*z
      endif
      f= (c*g)+(s*y)
      x=-(s*g)+(c*y)
      do 46 jj=1,m
      y=a(jj,j)
      z=a(jj,i)
      a(jj,j)= (y*c)+(z*s)
      a(jj,i)=-(y*s)+(z*c)
 46   continue
 47   continue
      rv1(l)=0.0
      rv1(k)=f
      w(k)=x
 48   continue
 3    continue
 49   continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 3#(11,1&#)6UK'VIka5..
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine smoothavcov()
c	find probability of a wet day given rf of previous day
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      include  'occr.inc'
      include  'pcond.inc'
      include  'pminmax.inc'
      real*4 xc(nearmax,nstnmax)
           real*4 xp1(nearmax,nstnmax)
      real*4 xpl(nvmax,nearmax,nstnmax)
      real*4 pprev(nvrmax,nstnmax)
      real*4 pmin(nvrmax,nstnmax),pmax(nvrmax,nstnmax)
      real*4 ppavl(2,nvmax,nstnmax),pavl(2,nvmax,nstnmax)
      real*4 ppcovl(2,nvmax,nvmax,nstnmax),
     1       pcovl(2,nvmax,nvmax,nstnmax)
      integer i,j,l,kk(nstnmax),j1,nk
      integer ip,jp,lp
	
c	write(*,*)'Calculating means and variances of long term variables'

      do j=1,nout
c	write(*,*)'.....for month ',j

      is=iseas(j,isn)
      do l=1,nday(j)
c	write(*,*)'............for day ',l
c	if(j.eq.2.and.l.eq.29)goto 111

      do nk=1,nstn
      do i=1,nearmax
      xc(i,nk)=0.0
      xp1(i,nk)=0.0
      do j1=1,nvrmax
      xpl(j1,i,nk)=0.0
      enddo
      enddo
      do jj=1,nlon
      pmin(jj,nk)=10000.0
      pmax(jj,nk)=-10000.0
      enddo
      kk(nk)=0
      nyear=nyrs(nk)
      ns=nstrt(nk)-1
      do 11 i=iyrst,nyear,1
      if(j.eq.2)then
      call daycount(ns,i,lc)
      if(l.gt.lc)goto 11
      endif
      
      li=l-iband-1
      
      do 10 jh=1,iband*2+1
      li=li+1
      ic=i
      jc=j
      lc=li
      if(li.le.2)call day_neg(ic,jc,lc,1,nout,ns,2,nday,monmax,indx)
      if(li.gt.2)call day_pos(ic,jc,lc,nout,ns,nday,nyear,monmax,indx)
      if(indx.eq.1.)goto 10
c	check for prev day
      lp=lc-1
      jp=jc
      ip=ic
      if(lp.lt.1)then
      if(ip.le.1.and.jp.eq.1)goto 10
      jp=jp-1
      if(jp.lt.1)then
      ip=ip-1
      jp=nout
      endif
      lp=nday(jp)
      if(jp.eq.2)call daycount(ns,ip,lp)
      endif


c	check for previous lon days average wetness state
      if(nlon.gt.0)then
      do jj=1,nlon
      jk=(nk-1)*nlon+jj
      aa=high(jk,ic,jc,lc)
      if(aa.lt.0.0)goto 10
      pprev(jj,nk)=aa
      if(pprev(jj,nk).lt.pmin(jj,nk))pmin(jj,nk)=pprev(jj,nk)
      if(pprev(jj,nk).gt.pmax(jj,nk))pmax(jj,nk)=pprev(jj,nk)
      enddo
      endif


      kk(nk)=kk(nk)+1
c	write(*,12)kk,i,j,l
 12   format(6i5)
c	

      xc(kk(nk),nk)=rf(nk,ic,jc,lc)
      xp1(kk(nk),nk)=rf(nk,ip,jp,lp)
      
      if(nlon.gt.0)then
      do jj=1,nlon
      xpl(jj,kk(nk),nk)=pprev(jj,nk)/sdl(jj,j,l,nk)
      enddo
      endif
 10   continue
 11   continue
      enddo
      if(nlon.gt.0)then
      iflag=1
      if(llag.eq.0)then
      call prevavsd0(xc,xpl,nstn,nlon,kk,rain,
     1     ppavl,ppcovl,pavl,pcovl,iflag,nearmax,nstnmax,nvmax)
      elseif(llag.eq.1)then
      call prevavsd1(xc,xp1,xpl,nstn,nlon,kk,rain,
     1     ppavl,ppcovl,pavl,pcovl,iflag,nearmax,nstnmax,nvmax)
      endif
      
      do nk=1,nstn
      j3=0
      do j1=1,nlon
      do j2=1,nlon
      j3=j3+1
      pp1lcov(j,l,j3,nk)=ppcovl(1,j1,j2,nk)
      pp2lcov(j,l,j3,nk)=ppcovl(2,j1,j2,nk)
      
      p1lcov(j,l,j3,nk)=pcovl(1,j1,j2,nk)
      p2lcov(j,l,j3,nk)=pcovl(2,j1,j2,nk)
      enddo
      enddo
      
      
      do j1=1,nlon
      plmin(j,l,j1,nk)=pmin(j1,nk)
      plmax(j,l,j1,nk)=pmax(j1,nk)
      enddo
      
      do j1=1,nlon+1
      pp1lav(j,l,j1,nk)=ppavl(1,j1,nk)
      pp2lav(j,l,j1,nk)=ppavl(2,j1,nk)
      
      p1lav(j,l,j1,nk)=pavl(1,j1,nk)
      p2lav(j,l,j1,nk)=pavl(2,j1,nk)
      enddo
      enddo
      nk=1
      write(13,*)'inverse covar matrix - longer time scale var'
      do j3=1,2
      do j1=1,nlon
      write(13,67)(pcovl(j3,j1,j2,nk),j2=1,nlon)
      enddo
      do j1=1,nlon
      write(13,67)(ppcovl(j3,j1,j2,nk),j2=1,nlon)
      enddo
      write(13,*)
      enddo
      write(13,*)'average details - longer time scale var'
      do j2=1,2
      write(13,67)(pavl(j2,j1,nk),j1=1,nlon+1)
      write(13,67)(ppavl(j2,j1,nk),j1=1,nlon+1)
      write(13,*)
      enddo
      
      endif
      
      
 67   format(50f15.5)

      enddo
      enddo
c	 correct for Feb 29
      is=iseas(2,isn)
      do nk=1,nstn
      j3=0
      do j1=1,nlon
      do j2=1,nlon
      j3=j3+1
      pp1lcov(2,29,j3,nk)=pp1lcov(2,28,j3,nk)
      pp2lcov(2,29,j3,nk)=pp2lcov(2,28,j3,nk)
      
      p1lcov(2,29,j3,nk)=p1lcov(2,28,j3,nk)
      p2lcov(2,29,j3,nk)=p2lcov(2,28,j3,nk)
      enddo
      enddo
      
      do j3=1,nlon+1
      pp1lav(2,29,j3,nk)=pp1lav(2,28,j3,nk)
      pp2lav(2,29,j3,nk)=pp2lav(2,28,j3,nk)
      
      p1lav(2,29,j3,nk)=p1lav(2,28,j3,nk)
      p2lav(2,29,j3,nk)=p2lav(2,28,j3,nk)
      enddo
      
      do j1=1,nlon
      plmin(2,29,j1,nk)=plmin(2,28,j1,nk)
      plmax(2,29,j1,nk)=plmax(2,28,j1,nk)
      enddo
      
      pro1(2,29,nk)=(pro1(2,28,nk)+pro1(3,1,nk))/2.0
      pro2(2,29,nk)=(pro2(2,28,nk)+pro2(3,1,nk))/2.0
      enddo
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rank1(xc,xp,ixn,nn,nstn,nearmax,nstnmax)
      real*4 xc(nearmax,nstnmax),xp(nearmax,nstnmax)
      integer nn(nstnmax),ixn(nearmax,nstnmax)
      
      do nk=1,nstn
      do i=1,nn(nk)-1
      do j=i+1,nn(nk)
      if(xc(j,nk).gt.xc(i,nk))then
      a=xc(i,nk)
      xc(i,nk)=xc(j,nk)
      xc(j,nk)=a
      
      a=xp(i,nk)
      xp(i,nk)=xp(j,nk)
      xp(j,nk)=a
      
      ia=ixn(i,nk)
      ixn(i,nk)=ixn(j,nk)
      ixn(j,nk)=ia
      
      endif
      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rf_amt_store()
c	Store variables for rf amount downscaling
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      include  'rf_amt.inc'
      real*4 xc(nearmax1,nstnmax),xl(nearmax1,nstnmax)
           real*4 w(nvmax)
           real*4 xp(nearmax1,nstnmax),xrc(nearmax1),yrc(nearmax1,nvmax)
      integer ixn(nearmax1,nstnmax),irhn(nstnmax)
      integer nn(nstnmax)
      
      jk=0
      do j=1,nout
c	write(*,*)'Storing variables for rf_amount .....for month ',j

      is=iseas(j,isn)
      do l=1,nday(j)
      jk=jk+1
      do nk=1,nstn
      do i=1,nearmax1
      xc(i,nk)=0.0
      xp(i,nk)=0.0
      ixn(i,nk)=0
      enddo
      enddo
      nn=0
      do nk=1,nstn
      nyear=nyrs(nk)
      ns=nstrt(nk)
      do 11 i=1,nyear,1
      if(j.eq.2)then
      call daycount(ns,i,lc)
      if(l.gt.lc)goto 11
      endif
      
      li=l-iband-1
      
      do 10 jh=1,iband*2+1
      li=li+1
      ic=i
      jc=j
      lc=li
      if(li.le.2)call day_neg(ic,jc,lc,1,nout,ns,2,nday,monmax,indx)
      if(li.gt.2)call day_pos(ic,jc,lc,nout,ns,nday,nyear,monmax,indx)
      if(indx.eq.1.)goto 10
c	check for prev day
      lp=lc-1
      jp=jc
      ip=ic
      if(lp.lt.1)then
      if(ip.le.1.and.jp.eq.1)goto 10
      jp=jp-1
      if(jp.lt.1)then
      ip=ip-1
      jp=nout
      endif
      lp=nday(jp)
      if(jp.eq.2)call daycount(ns,ip,lp)
      endif

c	check for the next day

      in=ic
      jn=jc
      ln=lc+1
      nd=nday(jn)
      if(jn.eq.2)call daycount(ns,in,nd)
c   check for windows of last few days of the last year
      if(in.eq.nyear.and.jn.eq.nout.and.jn.eq.nout.and.ln.gt.nd)goto 10
      if(ln.gt.nd)then
      jn=jn+1
      ln=1
      if(jn.gt.nout)then
      in=in+1
      jn=1
      endif
      nd=nday(jn)
      if(in.gt.nyear)goto 10
      endif
      irhn(nk)=0
      if(rf(nk,in,jn,ln).ge.rain)irhn(nk)=1

      if(rf(nk,ic,jc,lc).ge.rain)then
      nn(nk)=nn(nk)+1
c	if(nn(nk).gt.nearmax1)write(*,*)' increase dimension nearmax1'
      xc(nn(nk),nk)=rf(nk,ic,jc,lc)
      xp(nn(nk),nk)=rf(nk,ip,jp,lp)
      ixn(nn(nk),nk)=irhn(nk)
      endif

 10	  continue
 11	  continue
      enddo

      call rank1(xc,xp,ixn,nn,nstn,nearmax1,nstnmax)
      do nk=1,nstn
      do ii=1,nn(nk)
      rfcur(ii,j,l,nk)=xc(ii,nk)
      rfprv(ii,j,l,nk)=xp(ii,nk)
      irfnxt(ii,j,l,nk)=ixn(ii,nk)
      enddo
      rfcur(nn(nk)+1,j,l,nk)=-999.0

c	find weights of conditioning variables
c	to normalise rainfall, take square root of values
      stnwt1(nk,j,l,1)=1.0
      stnwt2(nk,j,l,1)=1.0
      stnwt3(nk,j,l,1)=1.0
      stnwt4(nk,j,l,1)=1.0
      enddo
      
      enddo
      enddo
      close(unit=13)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c	subroutine to compute different number of wet spells
c	one day to 10 days, length of maximum wet and dry spells
c	 and solitary wet days (dry days on either side)
      subroutine sur(rf,ny,nstn,spell,rfspl,rfspn,rfspnd,
     1               sol,nout,im,ns,rain,wmax,dmax)
      include  'para.inc'
      dimension rf(nstnmax,nyrmax,monmax,ndmax),sps(nyrmax,15000)
      dimension spw(nyrmax,15000),kmax(nyrmax),spell(nstnmax,15,nyrmax),
     1         amtmax(nyrmax),ar(nyrmax),am(nyrmax),im(monmax),
     2         rfspl(nstnmax,15),rfspn(nstnmax,15),
     3         rfspnd(nstnmax,15),sol(nstnmax),wmax(nstnmax),
     4         dmax(nstnmax)
      
      real rain,sum
      spell=0.0
      do nk=1,nstn
      wmax(nk)=0.0
      dmax(nk)=0.0
      do i=1,15
      rfspl(nk,i)=0.0
      rfspn(nk,i)=0.0
      rfspnd(nk,i)=0.0
      enddo
      do i=1,ny
      do k=1,5000
      sps(i,k)=0.0
      spw(i,k)=0.0
      enddo
      enddo

c	for wet spells ending in each year
      do i=1,ny
      kmax(i)=0
      amtmax(i)=0.0
      ar(i)=0.0
      am(i)=0.0
      enddo
      
      kk=0
      rfs=0.0
      do i=1,ny
      kst=1
      do j=1,nout
      
      md=im(j)
      if(j.eq.2)call daycount(ns,i,md)
      do l=1,md
      if(rf(nk,i,j,l).gt.amtmax(i))amtmax(i)=rf(nk,i,j,l)
      if(rf(nk,i,j,l).ge.rain)ar(i)=ar(i)+1.0
      am(i)=am(i)+rf(nk,i,j,l)
      
      if(rf(nk,i,j,l).lt.rain)then
      if(kk.eq.0)goto 101
      spw(i,kk)=spw(i,kk)+1.0
c	store consecutive per wet day amounts for 1, 2-4,5-8, 8-12 and >12days
      rfs=rfs/float(kk)
      if(kk.eq.1)rfspl(nk,1)=rfspl(nk,1)+rfs
      if(kk.gt.1.and.kk.lt.5)rfspl(nk,2)=rfspl(nk,2)+rfs
      if(kk.gt.4.and.kk.lt.8)rfspl(nk,3)=rfspl(nk,3)+rfs
      if(kk.gt.7)rfspl(nk,4)=rfspl(nk,4)+rfs
      if(kk.eq.1)rfspn(nk,1)=rfspn(nk,1)+1.0
      if(kk.gt.1.and.kk.lt.5)rfspn(nk,2)=rfspn(nk,2)+1.0
      if(kk.gt.4.and.kk.lt.8)rfspn(nk,3)=rfspn(nk,3)+1.0
      if(kk.gt.7)rfspn(nk,4)=rfspn(nk,4)+1.0
      if(kk.gt.kmax(i))kmax(i)=kk
c	endif
      rfs=0.0 
      kk=0
      kst=j
      goto 101
      endif
      kk=kk+1
      rfs=rfs+rf(nk,i,j,l)
 101  continue
      enddo
      enddo
      enddo
      do i=1,ny
      spell(nk,1,i)=amtmax(i)
       spell(nk,2,i)=ar(i)
       spell(nk,3,i)=am(i)
       spell(nk,4,i)=float(kmax(i))
      sum=0.0
      do k=2,kmax(i)
      sum=sum+spw(i,k)
      enddo
      spell(nk,5,i)=sum
      if(float(kmax(i)).gt.wmax(nk))wmax(nk)=float(kmax(i))
      enddo

c	for wet spell amounts
      do i=1,15
      if(rfspn(nk,i).gt.0.0)rfspl(nk,i)=rfspl(nk,i)/rfspn(nk,i)
      rfspn(nk,i)=rfspn(nk,i)/float(ny)
      enddo

c	for dry spell

      do i=1,ny
      do k=1,5000
      sps(i,k)=0.0
      spw(i,k)=0.0
      enddo
      
      kmax(i)=0
      enddo
      kk=0
      do i=1,ny
      kst=1
      do j=1,nout
      md=im(j)
      if(j.eq.2)call daycount(ns,i,md)
      do l=1,md
      if(rf(nk,i,j,l).ge.rain)then
      if(kk.eq.0)goto 102
      spw(i,kk)=spw(i,kk)+1.0
      if(kk.eq.1)rfspnd(nk,1)=rfspnd(nk,1)+1.0
      if(kk.gt.1.and.kk.le.9)rfspnd(nk,2)=rfspnd(nk,2)+1.0
      if(kk.gt.9.and.kk.le.18)rfspnd(nk,3)=rfspnd(nk,3)+1.0
      if(kk.gt.18)rfspnd(nk,4)=rfspnd(nk,4)+1.0
      if(kk.gt.kmax(i))kmax(i)=kk
c	endif
      kk=0
      kst=j
      goto 102
      endif
      kk=kk+1
 102  continue
      enddo
      enddo
      enddo
      do i=1,ny
        spell(nk,6,i)=float(kmax(i))
      sum=0.0
      do k=2,kmax(i)
      sum=sum+spw(i,k)
      enddo
      spell(nk,7,i)=sum
      if(float(kmax(i)).gt.dmax(nk))dmax(nk)=float(kmax(i))
      enddo
      
      do i=1,15
      rfspnd(nk,i)=rfspnd(nk,i)/float(ny)
      enddo


c	for solitary wet days
      sol(nk)=0.0
      in=1
      jn=1
      in=1
      do i=1,ny
      do j=1,nout
      md=im(j)
      if(j.eq.2)call daycount(ns,i,md)			  
      do l=1,md
      if(i.eq.1.and.j.eq.1.and.j.eq.1.and.l.eq.1)goto 301
      if(i.eq.ny.and.j.eq.12.and.j.eq.12.and.md.eq.31)goto 301
      ln=l+1
      if(ln.gt.md)then
      jn=jn+1
      ln=1
      endif
      if(jn.gt.12)then
      in=in+1
      jn=1
      endif
      rff=rf(nk,in,jn,ln)
      if(rf(nk,i,j,l).ge.rain.and.rfp.lt.rain.
     1 and.rf(nk,i,j,l).ge.rain.and.rff.lt.rain)then
      sol(nk)=sol(nk)+1.0
      endif
 301  rfp=rf(nk,i,j,l)
      enddo
      enddo
      enddo
      sol(nk)=sol(nk)/float(ny)
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine BASIC1(XX,AVER,SD,N,mm)
      include  'para.inc'
      real*4 xx(nyrmax,monmax),aver(monmax),sd(monmax),sx,sssx,rn,a
      integer n,i
C       THIS SUBROUTINE CALCULATES THE MEAN AND STANDARD
C       DEVIATION OF AN INPUT SEARIES OF LENGTH N
      do j=1,mm
             SX=0.0
             SSSX=0.0
             RN=FLOAT(N)
c        IF(RN.EQ.0.0)WRITE(*,41)
      DO 21 I=1,N
       A=XX(I,j)
       SX=SX+A
       SSSX=SSSX+A*A
 21   CONTINUE
      AVER(j)=SX/RN
      SD(j)=SSSX-SX*AVER(j)
c	if(sd(j).eq.0.0)write(*,*) 'sd is 0.0'
      if(sd(j).lt.0.0001)sd(j)=0.0001
c        IF((RN-1.0).EQ.0.0)WRITE(*,41)
      SD(j)=sqrt(SD(j)/(RN-1.0))
 41   FORMAT(10X,'DATA LENGTH MUST BE GREATER THAN  2')
      enddo
      RETURN
      END
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rank_stat(xx,np,ny,nstn)
        include  'para.inc'
      dimension xx(nstnmax,15,nyrmax)
      do nk=1,nstn
      do k=1,np
      do i=1,ny-1
      do j=i+1,ny
      if(xx(nk,k,j).gt.xx(nk,k,i))then
      a=xx(nk,k,i)
      xx(nk,k,i)=xx(nk,k,j)
      xx(nk,k,j)=a
      endif
      enddo
      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stat_h(rf,nstn,ns,nyear,rain,nout,nday,isn)
      include 'para.inc'
      include 'stat.inc'
      real*4 rf(nstnmax,nyrmax,monmax,ndmax)
      integer nday(monmax)
c	subroutine to calculate statistics
      nh=7
      call sur(rf,nyear,nstn,spellh,rfsplh,rfspnh,rfspndh,
     1               solh,nout,nday,ns,rain,wmaxh,dmaxh)
      call dmon(rf,nyear,nstn,rwavsh,rwsdsh,rfavmh,rfsdmh,
     1          rfavsh,rfsdsh,roavmh,rosdmh,roavsh,rosdsh,nout,nday,
     2          ns,rain,isn)
      call percentile(rf,nyear,nstn,nout,nday,hperh,rain,isn,ns)
      call rank_stat(spellh,nh,nyear,nstn)
      call logodd(rf,nstn,nyear,crsodh,nday,ns,nout,rain,isn)
      call crscors(rf,nyear,nday,ns,nstn,crsoh,nout,rain,isn,0)
      call crscors(rf,nyear,nday,ns,nstn,crsrh,nout,rain,isn,1)
      call crscor(rf,nyear,nday,ns,nstn,crsrdwh,nout,rain,isn,1)
      call crscor(rf,nyear,nday,ns,nstn,crsrdah,nout,rain,isn,2)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine dmon(hist,ny,nstn,rwavs,rwsds,rfavm,rfsdm,
     1          rfavs,rfsds,roavm,rosdm,roavs,rosds,nout,im,ns,
     2             rain,isn)
      include  'para.inc'
      dimension hist(nstnmax,nyrmax,monmax,ndmax),av1(monmax),
     1          sd1(monmax),xm(nyrmax,monmax),xs(nyrmax,monmax),
     2          xsr(nyrmax,monmax),xso(nyrmax,monmax),im(monmax),
     3          ym(nyrmax,monmax),wm(nyrmax,monmax)
      dimension ss1(isnmax),ss2(isnmax),
     1          rfavm(nstnmax,monmax),rfsdm(nstnmax,monmax),
     2          roavm(nstnmax,monmax),rosdm(nstnmax,monmax),
     3          rwavm(nstnmax,monmax),rwsdm(nstnmax,monmax),
     4          rfavs(nstnmax,isnmax),rfsds(nstnmax,isnmax),
     5          roavs(nstnmax,isnmax),rosds(nstnmax,isnmax),
     6          rwavs(nstnmax,isnmax),rwsds(nstnmax,isnmax)
      real rain
c	calculates mean, sd and no of wet days of daily monthly seasonal and annual rainfall totals
c	calculate wet day rainfall amount per year and compute mean and sd and store
      do nk=1,nstn
      do i=1,ny
      ss1=0.0
      ss2=0.0
      do j=1,nout
      is=iseas(j,isn)
      md=im(j)
      if(j.eq.2)call daycount(ns,i,md)
      sm=0.0
      dm=0.0
      wm=0.0
      do l=1,md
      ss1(is)=ss1(is)+hist(nk,i,j,l)
      if(hist(nk,i,j,l).ge.rain)ss2(is)=ss2(is)+1.0
      ss1(isn+1)=ss1(isn+1)+hist(nk,i,j,l)
      if(hist(nk,i,j,l).ge.rain)ss2(isn+1)=ss2(isn+1)+1.0
      sm=sm+hist(nk,i,j,l)
      if(hist(nk,i,j,l).ge.rain)dm=dm+1.0
      enddo
      xm(i,j)=sm
      ym(i,j)=dm
      if(sm.gt.0)wm(i,j)=dm/sm
      enddo
      do is=1,isn+1
      xsr(i,is)=ss1(is)
      xso(i,is)=ss2(is)
      xs(i,is)=0.0
      if(ss2(is).gt.0.0)xs(i,is)=ss1(is)/ss2(is)
      enddo
      enddo
      call basic1(xm,av1,sd1,ny,12)
      do j=1,12
      rfavm(nk,j)=av1(j)
      rfsdm(nk,j)=sd1(j)
      enddo
      call basic1(ym,av1,sd1,ny,12)
      do j=1,12
      roavm(nk,j)=av1(j)
      rosdm(nk,j)=sd1(j)
      enddo
      call basic1(wm,av1,sd1,ny,12)
      do j=1,12
      rwavm(nk,j)=av1(j)
      rwsdm(nk,j)=sd1(j)
      enddo
      
      call basic1(xsr,av1,sd1,ny,isn+1)
      do is=1,isn+1
      rfavs(nk,is)=av1(is)
      rfsds(nk,is)=sd1(is)
      enddo
      call basic1(xso,av1,sd1,ny,isn+1)
      do is=1,isn+1
      roavs(nk,is)=av1(is)
      rosds(nk,is)=sd1(is)
      enddo
      call basic1(xs,av1,sd1,ny,isn+1)
      do is=1,isn+1
      rwavs(nk,is)=av1(is)
      rwsds(nk,is)=sd1(is)
      enddo
      enddo
      return
      end  
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine percentile(hist,ny,nstn,nout,im,hper,rain,isn,ns)
c	calculates 5th precentile of daily rainfall
      include 'para.inc'
      parameter (ntot=nyrmax*monmax*ndmax)
      dimension hist(nstnmax,nyrmax,monmax,ndmax),x(ntot),
     1         hper(nstnmax,5)
      dimension im(nout)
      
      hper=0.0
      do nk=1,nstn
      do is=1,isn+1
        ir=0
      do i=1,ny
      do j=1,nout
      if(is.ne.iseas(j,isn).and.is.le.isn)goto 10
      md=im(j)
      if(j.eq.2)call daycount(ns,i,md)
      do l=1,md
      if(hist(nk,i,j,l).ge.rain)then
      ir=ir+1
      x(ir)=hist(nk,i,j,l)
      endif
      enddo
 10   continue
      enddo
      enddo
      do i=1,ir-1
      do j=i+1,ir
      if(x(j).gt.x(i))then
      a=x(i)
      x(i)=x(j)
      x(j)=a
      endif
      enddo
      enddo
      i=int(0.05*float(ir)+0.5)
      hper(nk,is)=x(i)
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine logodd(colv,nstn,ny,alg,nday,ns,nout,rain,isn)
      include 'para.inc'
      dimension colv(nstnmax,nyrmax,monmax,ndmax),
     1 p(isnmax,5,nstnmax,nstnmax),alg(isnmax,nstnmax,nstnmax)
      integer nday(monmax)
      
      do is=1,isn+1
      do i=1,5
      do j=1,nstn-1
      do k=j+1,nstn
      p(is,i,j,k)=1.0
      enddo
      enddo
      enddo
      enddo
      do i=1,ny
      do j2=1,nout
      is=iseas(j2,isn)
      nd=nday(j2)
      if(j2.eq.2)call daycount(ns,i,nd)
      do l=1,nd
      do j=1,nstn-1
      do j1=j+1,nstn
      if(colv(j,i,j2,l).lt.rain.and.colv(j1,i,j2,l).lt.rain)
     1  p(is,1,j,j1)=p(is,1,j,j1)+1.0
      if(colv(j,i,j2,l).lt.rain.and.colv(j1,i,j2,l).ge.rain)
     1  p(is,2,j,j1)=p(is,2,j,j1)+1.0
      if(colv(j,i,j2,l).ge.rain.and.colv(j1,i,j2,l).ge.rain)
     1  p(is,3,j,j1)=p(is,3,j,j1)+1.0
      if(colv(j,i,j2,l).ge.rain.and.colv(j1,i,j2,l).lt.rain)
     1  p(is,4,j,j1)=p(is,4,j,j1)+1.0
      
      if(colv(j,i,j2,l).lt.rain.and.colv(j1,i,j2,l).lt.rain)
     1  p(isn+1,1,j,j1)=p(isn+1,1,j,j1)+1.0
      if(colv(j,i,j2,l).lt.rain.and.colv(j1,i,j2,l).ge.rain)
     1  p(isn+1,2,j,j1)=p(isn+1,2,j,j1)+1.0
      if(colv(j,i,j2,l).ge.rain.and.colv(j1,i,j2,l).ge.rain)
     1  p(isn+1,3,j,j1)=p(isn+1,3,j,j1)+1.0
      if(colv(j,i,j2,l).ge.rain.and.colv(j1,i,j2,l).lt.rain)
     1  p(isn+1,4,j,j1)=p(isn+1,4,j,j1)+1.0
      
      enddo
      enddo
      enddo
 100  continue
      enddo
      enddo

      do is=1,isn+1
      do j=1,nstn-1
      do j1=j+1,nstn
      alg(is,j,j1)=0.0
      alg(is,j,j1)=(p(is,1,j,j1)*p(is,3,j,j1))/
     1             (p(is,2,j,j1)*p(is,4,j,j1))
      alg(is,j,j1)=log(alg(is,j,j1))
 20   continue
      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine crscors(hist,ny,nday,ns,nstn,crs,nout,rain,isn,ind)
      include 'para.inc'
      real*4 hist(nstnmax,nyrmax,monmax,ndmax),
     1       crs(isnmax,nstnmax,nstnmax),sumd(isnmax)
      real*4 x(nstnmax,nyrmax,isnmax),xx1(nyrmax),yy1(nyrmax)
      integer nday(nout)
      
      do k=1,nstn
      do i=1,ny
      sumd=0.0
      do jj=1,nout
      is=iseas(jj,isn)
      md=nday(jj)
      if(jj.eq.2)call daycount(ns,i,md)
      do i1=1,md
      if(ind.eq.0)then
      if(hist(k,i,jj,i1).ge.rain)sumd(is)=sumd(is)+1.0
      if(hist(k,i,jj,i1).ge.rain)sumd(isn+1)=sumd(isn+1)+1.0
      else
      sumd(is)=sumd(is)+hist(k,i,jj,i1)
      sumd(isn+1)=sumd(isn+1)+hist(k,i,jj,i1)
      endif
      enddo
      enddo
      do is=1,isn+1
      x(k,i,is)=sumd(is)
      enddo
      enddo
      enddo
      
      do k1=1,nstn-1
      do k2=k1+1,nstn
      do j=1,isn+1
      do i=1,ny
      xx1(i)=x(k1,i,j)
      yy1(i)=x(k2,i,j)
      enddo
      call basic(xx1,ax,sx,ny)
      call basic(yy1,ay,sy,ny)
      if(sx.eq.0.0.or.sy.eq.0.0)then
      crs(j,k1,k2)=0.0
      goto 10
      endif
      s1=0.0
      do i=1,ny
      s1=s1+(xx1(i)-ax)*(yy1(i)-ay)
      enddo
      s1=s1/(float(ny-1)*sx*sy)
      crs(j,k1,k2)=s1

 10   continue
      enddo
      enddo
      enddo
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine crscor(hist,ny,nday,ns,nstn,crs,nout,rain,isn,ind)
c	ind=1 correlation of wet days only else all days
      include 'para.inc'
      parameter(nnmax=monmax*ndmax*nyrmax)
      real*4 hist(nstnmax,nyrmax,monmax,ndmax),
     1       crs(isnmax,nstnmax,nstnmax)
      real*4 x(nnmax),y(nnmax),xx(isnmax,nnmax),yy(isnmax,nnmax)
      integer nday(monmax),ij(isnmax)
      
      
      do k=1,nstn-1
      do j=k+1,nstn
      ij=0
      
      do i=1,ny
      do jj=1,nout
      is=iseas(jj,isn)
      md=nday(jj)
      if(jj.eq.2)call daycount(ns,i,md)
      do l=1,md
      if(ind.eq.1)then
      if(hist(k,i,jj,l).lt.rain.or.hist(j,i,jj,l).lt.rain)goto 10
      endif
      ij(is)=ij(is)+1
      xx(is,ij(is))=hist(k,i,jj,l)
      yy(is,ij(is))=hist(j,i,jj,l) 
      ij(isn+1)=ij(isn+1)+1
      xx(isn+1,ij(isn+1))=hist(k,i,jj,l)
      yy(isn+1,ij(isn+1))=hist(j,i,jj,l) 
 10   continue

      enddo
      enddo
      enddo
      
      do is=1,isn+1
      do i=1,ij(is)
      x(i)=xx(is,i)
      y(i)=yy(is,i)
      enddo
      n=ij(is)
      call basic(x,ax,sx,n)
      call basic(y,ay,sy,n)
      sum=0.0
      do i=1,n
      sum=sum+(x(i)-ax)*(y(i)-ay)
      enddo
      sum=sum/(float(n-1)*sx*sy)
      crs(is,k,j)=sum
      enddo
      
      enddo
      enddo
      
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine firsttimeh(ist,nk,sumlon,ind)
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      real*4 sumlon(nvrmax,lonmax)
      ind=0
      do kk=1,nlon
      ln=lon(2,kk)
      im=ist
      jm=nout
      nd=nday(jm)
      ns=nstrt(nk)-1
      if(jm.eq.2)call daycount(ns,im,nd)
      ls=nd
 10   js=ls-ln
      if(js.le.0)then
      jm=jm-1
      if(jm.lt.1)then
      jm=nout
      im=im-1
      endif
      if(im.lt.1)then
      ind=1
      return
      endif
      nd=nday(jm)
      if(jm.eq.2)call daycount(ns,im,nd)
      ls=ls+nd
      goto 10
      endif
      js=js-1
      nd=nday(jm)
      if(jm.eq.2)call daycount(ns,im,nd)
      do l=1,ln
      js=js+1
      if(js.gt.nd)then
      js=js-nd
      jm=jm+1
      if(jm.gt.nout)then
      jm=1
      im=im+1
      endif
      nd=nday(jm)
      if(jm.eq.2)call daycount(ns,im,nd)
      endif
c	write(*,*)im,jm,js
c	
      sumlon(kk,l)=0.0
      if(rf(nk,im,jm,js).ge.rain)sumlon(kk,l)=1.0
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      Subroutine aloncov(pav,ppav,pcov,ppcov,det,ddet,
     1                 p1,llag,rain,nv,nk,j,l1)
      include  'para.inc'
      include  'pcond.inc'
      real*4 ppav(nvmax),pav(nvmax),p1
      real*4 pcov(nvmax,nvmax),
     +       ppcov(nvmax,nvmax)
       
      rn=rain
      if(llag.eq.0)then
       j3=0
      do j1=1,nv
      do j2=1,nv
      j3=j3+1
      ppcov(j1,j2)=pp1lcov(j,l1,j3,nk)
      pcov(j1,j2)=p1lcov(j,l1,j3,nk)
      enddo
      ppav(j1)=pp1lav(j,l1,j1,nk)
      pav(j1)=p1lav(j,l1,j1,nk)
      enddo
      ddet=pp1lav(j,l1,nv+1,nk)
      det=p1lav(j,l1,nv+1,nk)
      endif
      if(llag.eq.1)then
       j3=0
      do j1=1,nv
      do j2=1,nv
      j3=j3+1
      if(p1.lt.rain)ppcov(j1,j2)=pp1lcov(j,l1,j3,nk)
      if(p1.ge.rain)ppcov(j1,j2)=pp2lcov(j,l1,j3,nk)
      if(p1.lt.rain)pcov(j1,j2)=p1lcov(j,l1,j3,nk)
      if(p1.ge.rain)pcov(j1,j2)=p2lcov(j,l1,j3,nk)
      enddo
      if(p1.lt.rain)ppav(j1)=pp1lav(j,l1,j1,nk)
      if(p1.ge.rain)ppav(j1)=pp2lav(j,l1,j1,nk)
      if(p1.lt.rain)pav(j1)=p1lav(j,l1,j1,nk)
      if(p1.ge.rain)pav(j1)=p2lav(j,l1,j1,nk)
      enddo
      if(p1.lt.rain)ddet=pp1lav(j,l1,nv+1,nk)
      if(p1.ge.rain)ddet=pp2lav(j,l1,nv+1,nk)
      if(p1.lt.rain)det=p1lav(j,l1,nv+1,nk)
      if(p1.ge.rain)det=p2lav(j,l1,nv+1,nk)
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine check_lon(igen,sumlon,lon,nlon)
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      real*4 sumlon(nvrmax,lonmax)
      integer lon(2,nvrmax),igen
      do kk=1,nlon
      nn=lon(2,kk)
      
      if(nn.gt.1)then
      do l=1,nn-1
      sumlon(kk,l)=sumlon(kk,l+1)
      enddo
      endif
      sumlon(kk,nn)=float(igen)
      enddo
      
      do kk=1,nlon
      nn=lon(2,kk)-lon(1,kk)+1
      
      sum=0.0
      do l=1,nn
      sum=sum+sumlon(kk,l)
      enddo
      sumlon(kk,lon(2,kk)+1)=sum
      enddo
      return
c	write(*,20)sumlon(kk,nn+1)
c 20	format(20f7.2)
	
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine multiply(x,av,cov,nv,nvrmax,sss)
c      Implicit double precision (A-H,O-Z) 
      real*4 av(nvrmax),cov(nvrmax,nvrmax),x(nvrmax),sss
      real*8 ss, part, qf
      real*8 dav(nvrmax),dcov(nvrmax,nvrmax),dx(nvrmax)


      do k = 1,nv
      do l = 1,nv
      dcov(k,l)=dble(cov(k,l))
      enddo
      dav(k)=dble(av(k))
      dx(k)=dble(x(k))
      enddo

c Multivariate normal
      qf=0.d0
      do 10 k = 1,nv
      do 10 l = 1,nv

c	write(*,*)k,l,qf,cov(k,l),x(k),av(k),x(l),av(l)
 10   qf = qf + dcov(k,l)*(dx(k)-dav(k))*(dx(l)-dav(l))
c       atmf = (1/sqrt(2*pi))**natm * (1/sqrt(det)) * dexp(-.5*qf)
      part=-0.5d0*qf
      if(part.lt.-15.d0)part=-15.d0
      ss = dexp(part)
c       write(*,*) i,j,atmf
      sss=real(ss)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c subroutine calc_href
c Written by: Kin Luk, university of new south wales
c purpose: an option to calculate the bandwidth
c Input: nx - sample size
c        nv - no. of dimensions.
c
c Output: h - the bandwidth.
c
c           
c Calling subroutines: Nil.
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine calc_href(nx,nv,h)
c      Implicit double precision (A-H,O-Z)
      fac = (4.0/(float(nv)+2.0))**(1.0/(float(nv)+4.0))
      h = fac * float(nx)**(-1.0/(float(nv)+4.0))
c         write(3,*) 'h=', h
      return
      end
c
c  Numerical recipes function RAN1 for uniform random
c  variate generation
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rank_h(y,n,rn,nmax)
c	rank the given series from highest to lowest
c        Implicit double precision (A-H,O-Z)
      real*4 x(nmax),y(nmax)
      do i=1,n
      x(i)=y(i)
      enddo
      do i=1,n-1
      do j=i+1,n
      if(x(j).gt.x(i))then
      a=x(i)
      x(i)=x(j)
      x(j)=a
      endif
      enddo
      enddo
      i1=int(float(n)*0.25+0.05)
      i2=int(float(n)*0.75+0.05)
      rn=abs(x(i1)-x(i2))
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE hfracx(x,hgamma,n,nv,nmax,nvmax)
c      Implicit double precision (A-H,O-Z)
        parameter (pi = 3.1415926535897932385)
      real*4 x(nmax,nvmax),hgamma(nmax,nvmax),xx(nmax),xgamma
      real*4 z(nmax,nvmax),new,lambda
      hupper=3.5
      hlower=1.0/hupper

      fac = (4.0/(float(nv)+2.0))**(1.0/(float(nv)+4.0))
      href = fac * float(n)**(-1.0/(float(nv)+4.0))
c	href=href*0.80

      do j=1,nv
      do i=1,n
      xx(i)=x(i,j)
      hgamma(i,j)=href
      enddo
c	enddo
c	return



c	call basic(xx,av,sd,n)
      call rank_h(xx,n,sd,nmax)
      do i=1,n
      ii=0
      aup=xx(i)+0.5*sd
      alr=xx(i)-0.5*sd
      if(aup.gt.amx)then
      alr=alr-(aup-amx)
      aup=amx
      endif
      if(alr.lt.amn)then
      aup=aup+(amn-alr)
      alr=amn
      endif
      
      do jj=1,n
      if(xx(jj).gt.alr.and.xx(jj).lt.aup)ii=ii+1
      enddo
      hgamma(i,j)=(1.2-float(ii)/float(n))*href
      if(hgamma(i,j).gt.hupper*href)hgamma(i,j)=hupper*href
      if(hgamma(i,j).lt.hlower*href)hgamma(i,j)=hlower*href
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      SUBROUTINE hfracy(z,xgiven,hgamma,n,nmax)
c      Implicit double precision (A-H,O-Z)
      parameter (pi = 3.1415926535897932385)
      real*4 x(nmax),g(nmax),new,lambda,z(nmax)
      hupper=3.50
      hlower=1.0/hupper
      
      nv=1
      fac = (4.0/(float(nv)+2.0))**(1.0/(float(nv)+4.0))
      href = fac * float(n)**(-1.0/(float(nv)+4.0))
c	href=href*0.80
      hgamma=href
c	return
c	prewhiten the data
c	call basic(z,av,sd,n)
      call rank_h(z,n,sd,nmax)
      aup=xgiven+0.5*sd
      alr=xgiven-0.5*sd
      if(aup.gt.amx)then
      alr=alr-(aup-amx)
      aup=amx
      endif
      if(alr.lt.amn)then
      aup=aup+(amn-alr)
      alr=amn
      endif

      ii=0
      do j=1,n
      if(z(j).gt.alr.and.z(j).lt.aup)ii=ii+1
      enddo
      hgamma=(1.2-float(ii)/float(n))*href
      if(hgamma.gt.hupper*href)hgamma=hupper*href
      if(hgamma.lt.hlower*href)hgamma=hlower*href
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  subroutine estxcpdfstats
c
c  to estimate sxx, sxxi and wts for conditional probability density
c
c  written by
c    Ashish Sharma
c    CEE, UNSW
c    21/02/1999
c
c  1 - find covariance of x
c  2 - find sxxi
c  3 - find wts
c  4 - end
c
c  modified 2/10/2001 to consider uncertainty in observations
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine estxcpdfstats(x,nx,nv,xt,nn,weight,cov,
     1        wx,xcond,sxxi,h,wts,nvmax,nmax,hgammax)

c      Implicit double precision (A-H,O-Z)
        
      real*4 x(nmax,nvmax),xx(nmax,nvmax),
     1        xcond(nvmax),xt(nmax,nvmax),
     2        cov(nvmax,nvmax),xmn(nvmax),
     3        wx(nmax,nvmax),wxmn(nvmax),weight(nvmax),
     5        wcov(nvmax,nvmax),hgammax(nmax,nvmax)
      real*4 ratio,wts(nmax),sxxi(nvmax,nvmax),uutemp(nmax,nvmax),
     1    sumwt,uu1,uu

c        write(*,*)'1nx',nx,'nv',nv
      do 100 i=1,nv
      do 110 j=1,nv
      cov(i,j)=0.0
      wcov(i,j)=0.0
 110  continue
      xmn(i)=0.0
      wxmn(i)=0.0
 100  continue

      do 150 i=1,nn
      do 160 j=1,nv
      xmn(j)=xmn(j)+wx(i,j)*xt(i,j)
      wxmn(j)=wxmn(j)+wx(i,j)
 160  continue
 150  continue

c      write(*,*)'2nx',nx,'nv',nv
      do 170 i=1,nv
 170  if(wxmn(i).ne.0.0)xmn(i)=xmn(i)/wxmn(i)

      do 200 i=1,nn
      do 210 j=1,nv
      do 220 k=1,nv
      cov(j,k)=cov(j,k)+wx(i,j)*(xt(i,j)-xmn(j))*
     1            wx(i,k)*(xt(i,k)-xmn(k))
      wcov(j,k)=wcov(j,k)+wx(i,j)*wx(i,k)
 220  continue
 210  continue
 200  continue

c      write(*,*)'3nx',nx,'nv',nv
      do 230 j=1,nv
      do 240 k=1,nv

      if(wcov(j,k).ne.0.0)cov(j,k)=cov(j,k)/wcov(j,k)
      sxxi(j,k) = cov(j,k)
 240  continue
c      write(*,*)'cov',(cov(j,k),k=1,nv),nx,x(1,j),x(nx,j)
 230  continue
c covariance estimated - now to find inverse
c        write(*,*)'4nx',nx,'nv',nv
      call solve(sxxi,nv,nvmax,det)
c	do j=1,nv
c        write(*,*)'sxxi',(sxxi(j,k),k=1,nv),nx,x(1,j),x(nx,j)
c	enddo
c	
c       write(*,*)'5nx',nx,'nv',nv
        
c        write(*,*)'okstats',h*h
      call hfracx(x,hgammax,nx,nv,nmax,nvmax)
c now to find the weights
      hsq = h*h
        
      do 410 j1 = 1,nx
      do 420 j2 = 1,nv
      xx(j1,j2) = (xcond(j2) - x(j1,j2))*sqrt(weight(j2))
      uutemp(j1,j2) = 0.d0
 420  continue
 410  continue

c     write(*,*)'6nx',nx,'nv',nv
	  sumwt = 0.0
      do 500 j1 = 1,nx
      uu = 0.0
      wtmp = 1.0
      
      do 510 j2 = 1,nv
      do 520 j3 = 1,nv
	  hsq = hgammax(j1,j2)*hgammax(j1,j3)
      uutemp(j1,j2)=uutemp(j1,j2)+xx(j1,j3)*sxxi(j2,j3)/hsq

 520  continue
      uu = uu + uutemp(j1,j2)*xx(j1,j2)
      wtmp = wtmp*wx(j1,j2)
 510  continue

      uu1=uu/(2.0)
      if(uu1.gt.20.0)then
      uu1=0.0
      else
      uu1=exp(-uu1)
      endif


      wts(j1) = wtmp*uu1
      sumwt = sumwt + wts(j1)
 500  continue
c	write(*,*)'sumwt=',sumwt
c	
      do 600 j1 = 1,nx
	  ratio=0.0
	  if(sumwt.gt.0.0)ratio=wts(j1)/sumwt
      wts(j1) = ratio
 600  continue

	  if(sumwt.eq.0.0)then
      do 700 j1 = 1,nx
      wts(j1) = 1.0/float(nx)
 700  continue
      endif


c          do j1 = 1,nx
c	write(*,943)j1,wts(j1),xx(j1,nv),x(j1,nv),xcond(nv)
c943	format(i5,4f10.5)
c	enddo
c	


      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matmul1(A,B,C,nv,nvmax)
c      Implicit double precision (A-H,O-Z)
      real*4 A(nvmax,nvmax)
      real*4 B(nvmax), C(nvmax)

      do 100 i=1,nv
      C(i)=0.0
 100  continue
      do 200 i=1,nv
	  sum=0.0
      do 210 j=1,nv
      sum=sum+real(A(i,j))*B(j)
 210  continue
	  C(i)=sum
 200  continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine matmul2(A,B,C,nv,nvmax)
c      Implicit double precision (A-H,O-Z)
      real*4 A(nvmax), B(nvmax), C

      C=0.0

      do 200 i=1,nv
      C=C+A(i)*B(i)
 200  continue

      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  subroutine estscond
c
c to estimate conditional covariance and sxxixy
c
c  written by
c    ashish sharma
c    cee,unsw
c    22/02/1999
c    modified 3/10/2001 to incorporate uncertainty in observations
c
c
c  1 - estimate sxy
c  2 - estimate sxxixy
c  3 - estimate sxyxxixy
c  4 - estimate scond = sy - sxyxxixy
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine estscond(y,n,wy,x,nv,wx,
     1        sxxi,sxxixy,scond,nvmax,nmax)
c      Implicit double precision (A-H,O-Z)
        
       
      real*4 y(nmax), x(nmax,nvmax),
     1        sxxixy(nvmax),
     2        xmn(nvmax),sxy(nvmax),wxmn(nvmax),
     3        wy(nmax),wx(nmax,nvmax),wsxy(nvmax)
      real*4 sxxi(nvmax,nvmax)
      nx = n
      do 100 i=1,nv
      xmn(i)=0.0
      wxmn(i)=0.0
 100  continue
      ymn = 0.0
      wymn = 0.0

      do 150 i=1,nx
      do 160 j=1,nv
      xmn(j)=xmn(j)+wx(i,j)*x(i,j)
      wxmn(j)=wxmn(j)+wx(i,j)
 160  continue
      ymn = ymn+wy(i)*y(i)
      wymn = wymn+wy(i)
 150  continue
      do 170 i=1,nv
      xmn(i)=xmn(i)/wxmn(i)
      sxy(i) = 0.0
      wsxy(i) = 0.0
 170  continue
      ymn = ymn/wymn

      syy = 0.0
      wsy = 0.0
      do 200 i=1,nx
      syy = syy+(wy(i)*(y(i)-ymn))**2
      wsy = wsy+wy(i)**2
      do 210 j=1,nv
      sxy(j)=sxy(j)+wx(i,j)*(x(i,j)-xmn(j))*
     1 wy(i)*(y(i)-ymn)
      wsxy(j)=wsxy(j)+wx(i,j)*wy(i)
 210  continue
 200  continue
c        write(*,*)'sxy subroutine'
c        write(7,*)(xmn(i),i=1,nv),ymn
      do 230 j=1,nv
      if(wsxy(j).ne.0.0)sxy(j)=sxy(j)/wsxy(j)
c          write(*,*)'sxy',sxy(j)
 230  continue
	  syy=syy/wsy
      call matmul1(sxxi,sxy,sxxixy,nv,nvmax)
          
      call matmul2(sxy,sxxixy,scond,nv,nvmax)

      if(scond.eq.0.0)write(*,*)'check covar matrix'
      scond = syy - scond
      if(scond.lt.0.0)then
      write(*,*)'scond is -ive, considering as zero',scond,syy,nx
      scond=0.0
      endif
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c  subroutine psim
c
c  subroutine to simulate the probabilistic
c
c  written by
c        Ashish Sharma
c        CEE, UNSW
c        30/09/2002
c
c  NOTE - this is a modification of an earlier subroutine called egen
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine psim(y,n,x,nv,xcond,sxxixy,h,wp,
     1        cwts,scond,iseed,ygen,vv,xr,rain,nvmax,nmax,hgammax)
      
c      Implicit double precision (A-H,O-Z)
      real*4          y(nmax)
     1                , x(nmax,nvmax)
     2                , xcond(nvmax)
     3                , sxxixy(nvmax)
     4                , wp(nvmax)
     5                , hgammax(nmax,nvmax)
      real*4            cwts(nmax),rnum
             
      hused=h
c          rnum=RAN1(iseed)*cwts(n)


      rnum=vv*cwts(n)
      do 200 j=1,n
      if(cwts(j).ge.rnum)then
      iwt=j
      go to 300
      endif
 200  continue

 300  continue
      call hfracy(y,y(iwt),hgamma,n,nmax)
      hused=hgamma


      b=0.0
      do 310 j=1,nv
      b=b+(xcond(j)-x(iwt,j))*sxxixy(j)*wp(j)*hused/hgammax(iwt,j)
c          b=b+(xcond(j)-x(iwt,j))*sxxixy(j)*wp(j)
 310  continue
c	checking for the conditional mean, imposing limits on it if necessary
c	if(b.gt.1.5*y(1))b=1.5*y(1)
c	if(b.lt.y(n)/2)b=y(n)/2.0

c 3.3 - generating the realisation N(b,h*scond)
      bran=hused*sqrt(scond)*xr
      ygen=y(iwt)+b+bran

c	if(ygen.gt.2.0*y(iwt))ygen=2.0*y(iwt)
c	if(ygen.lt.y(iwt)/2.0)ygen=y(iwt)/2.0


c	checking for the very high value, imposing limits on it if necessary
      if(ygen.gt.y(1)*3.0)then
      if(b.gt.1.5*y(1))b=1.5*y(1)
      ygen=y(iwt)+b+bran
      endif


c 4 - checking for negative realisations
	  ii=0
 410  ii=ii+1
      if(ii.gt.5)return
      if(ygen.lt.0.0)then
      hused=hused/2.0
      bran=hused*sqrt(scond)*xr
      if(b.lt.0.0)b=b/2.0
      ygen=y(iwt)+b+bran
	  goto 410
c	ygen=y(n)
      endif
  
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine psimmain(y,n,x,nv,xcond,ys,iseed,vv,
     1                    xr,rain,weight,nmax,nvmax,ind)
             
c      Implicit double precision (A-H,O-Z)

      real*4          y(nmax)
     1                , x(nmax,nvmax)
     2                , xcond(nvmax)
     5                , sxxixy(nvmax)
     6                , wy(nmax)
     7                , wx(nmax,nvmax),hgammax(nmax,nvmax)
     8                , weight(nvmax),cov(nvmax,nvmax)

      real*4            cwts(nmax), wts(nmax), sxxi(nvmax,nvmax)


      call calc_href(n,nv,h)
      h=0.80*h
c	initialise some variables
      do i=1,n
      do j=1,nvmax
      wx(i,j)=1.0
      enddo
c	write(*,2)i,y(i),(x(i,j),j=1,nv)

      enddo
c	
 2    format(i5,10f7.1)

      do i=1,nmax
      wts(i)=0.0
      wy(i)=1.0
      enddo
 109  if(ind.eq.1)then
c	perform unconditional simulation
 10   iwt=int(vv*float(n)+1.0)
	  if(iwt.lt.1)iwt=1
	  if(iwt.gt.n)iwt=n
      if(n.lt.5)then
      ys=y(iwt)
      return
      endif
      call basic(y,av,sd,n)
      call hfracy(y,y(iwt),h,n,nmax)
      ys = y(iwt)+h*sd*xr
      if(ys.gt.y(1)*3.0)then
      ys=y(iwt)
      endif
      return
      endif

c	write(*,*)'1',n        
      call estxcpdfstats(x,n,nv,x,n,weight,cov,
     1        wx,xcond,sxxi,h,wts,nvmax,nmax,hgammax)

c	write(*,*)'2'     
      call estscond(y,n,wy,x,nv,wx,
     1        sxxi,sxxixy,scond,nvmax,nmax)

c	write(*,*)'3'
      do 100 i1=1,n
      cwts(i1)=wts(i1)*wy(i1)
 100  continue

      do 110 i1=2,n
      cwts(i1)=cwts(i1-1)+cwts(i1)
 110  continue


c	write(*,*)'4'          
      call psim(y,n,x,nv,xcond,sxxixy,h,weight,
     1        cwts,scond,iseed,ys,vv,xr,rain,nvmax,nmax,hgammax)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine stat_g(itr,rf_gen,nstn,nsg,ng,rain,nout,nday,isn)
      include 'para.inc'
      include 'stat.inc'
      real*4 rf_gen(nstnmax,nyrmax,monmax,ndmax)
      integer nday(monmax)
c	subroutine to calculate statistics
      nh=7
      call sur(rf_gen,ng,nstn,spell,rfspl,rfspn,rfspnd,
     1               sol,nout,nday,nsg,rain,wmax,dmax)
      call dmon(rf_gen,ng,nstn,rwavs,rwsds,rfavm,rfsdm,
     1          rfavs,rfsds,roavm,rosdm,roavs,rosds,nout,nday,nsg,
     1             rain,isn)
      call percentile(rf_gen,ng,nstn,nout,nday,hper,rain,isn,nsg)

      goto 11
c	open(9,file='temp.out')
c	do i=1,ng
c	do j=1,nout
c	md=nday(j)
c	if(j.eq.2)call daycount(nsg,i,md)
c	do l=1,md
c	write(9,901)i+nsg,j,l,rf_gen(1,i,j,l)
c	enddo
c	enddo
c	enddo
c	stop
 901  format(3i5,50(1x,f6.1))



 11   do nk=1,nstn
      do i=1,15
      rfsplg(itr,nk,i)=rfspl(nk,i)
      rfspng(itr,nk,i)=rfspn(nk,i)
      rfspndg(itr,nk,i)=rfspnd(nk,i)
      enddo
      solg(itr,nk)=sol(nk)
      enddo
      call rank_stat(spell,nh,ng,nstn)
      
      
      do k=1,nh
      do i=1,ng
      do nk=1,nstn
      spellg(itr,nk,k,i)=spell(nk,k,i)
      enddo	
      enddo
      enddo
      
      do nk=1,nstn
      do j=1,nout
      rfavmg(itr,nk,j)=rfavm(nk,j)
      rfsdmg(itr,nk,j)=rfsdm(nk,j)
      roavmg(itr,nk,j)=roavm(nk,j)
      rosdmg(itr,nk,j)=rosdm(nk,j)
      rwavmg(itr,nk,j)=rwavm(nk,j)
      rwsdmg(itr,nk,j)=rwsdm(nk,j)
      enddo
      do j=1,isn+1
      rfavsg(itr,nk,j)=rfavs(nk,j)
      rfsdsg(itr,nk,j)=rfsds(nk,j)
      roavsg(itr,nk,j)=roavs(nk,j)
      rosdsg(itr,nk,j)=rosds(nk,j)
      rwavsg(itr,nk,j)=rwavs(nk,j)
      rwsdsg(itr,nk,j)=rwsds(nk,j)
      hperg(itr,nk,j)=hper(nk,j)
      enddo
      wmaxg(itr,nk)=wmax(nk)
      dmaxg(itr,nk)=dmax(nk)
      enddo
      
      call logodd(rf_gen,nstn,ng,crsod,nday,nsg,nout,rain,isn)
      call crscors(rf_gen,ng,nday,nsg,nstn,crso,nout,rain,isn,0)
      call crscors(rf_gen,ng,nday,nsg,nstn,crsr,nout,rain,isn,1)
      call crscor(rf_gen,ng,nday,nsg,nstn,crsrdw,nout,rain,isn,1)
      call crscor(rf_gen,ng,nday,nsg,nstn,crsrda,nout,rain,isn,2)
      
      do is=1,isn+1
      do i=1,nstn
      do j=1,nstn
      crsodg(itr,is,i,j)=crsod(is,i,j)
      crsog(itr,is,i,j)=crso(is,i,j)
      crsrg(itr,is,i,j)=crsr(is,i,j)
      crsrdwg(itr,is,i,j)=crsrdw(is,i,j)
      crsrdag(itr,is,i,j)=crsrda(is,i,j)
      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine simulate (ws,avrf,tavrf)
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      include  'occr.inc'
      include  'pcond.inc'
      include  'pminmax.inc'
      real*4 pi
      parameter (pi = 3.1415926535897932385)
           real*4 pw, p1,tavrf
           real*4 prevsum(nvrmax),avrf(nstnmax),
     +       sumlon(nvrmax,lonmax),gen
      real*4 det,ddet,rn,rf_gen(nstnmax,nyrmax,monmax,ndmax)
      real*4 ppavl(nvmax),pavl(nvmax),prevl(nvrmax)
      real*4 wgen(5),whist(5,nstnmax),prevm(nvrmax),
     1       whamt(5,nstnmax),wgenr(5),sumho(5),sumhr(5),
     2       	sumg(5),sumgr(5),prevlong(nvmax),ws(nstnmax)


      real*4 pn,pd
      real*4 pcovl(nvmax,nvmax),ppcovl(nvmax,nvmax)
      integer nk,j1,i,j,l1,nd,jfx
      integer igen,icur,inext
      real*4 ip1,ip2,ip3
      rn=rain


      call wetdays(whist,whamt)

c	pick a station to borrow statistics
      vv=ran1(iseed)
      do kk=1,5
      if(vv.le.ws(kk))goto 531
      enddo
 531  nk=kk
      nyear=nyrs(nk)

 11   ist=ran1(iseed)*nyear+1
      if(ist.lt.iyrst.or.ist.gt.nyear)goto 11
      nd=nday(nout)

      ip1=0.0
      ip2=0.0
      ip3=0.0
      if(rf(nk,ist,nout,nd).ge.rain)ip1=1.0
      if(rf(nk,ist,nout,nd-1).ge.rain)ip2=1.0
      if(rf(nk,ist,nout,nd-2).ge.rain)ip3=1.0

      ind=0
      if(nlon.gt.0)call firsttimeh(ist,nk,sumlon,ind)
      if(ind.eq.1)goto 11
      igen=0
      if(rf(nk,ist,nout,31).ge.rain)igen=1
      p1=rf(nk,ist,nout,31)
      gen=rf(nk,ist,nout,31)


      do it=1,nsim
      iic=0
 

      wgen=0.0
      wgenr=0.0
      do i=1,ng
      sumg=0.0
      sumamt=0.0
      sumgr=0.0
      do j=1,nout
      is=iseas(j,isn)

      nd=nday(j)
      if(j.eq.2)call daycount(nsg,i,nd)
      do l1=1,nd
c	pick a station to borrow statistics
      vv=ran1(iseed)
      do kk=1,5
      if(vv.le.ws(kk))goto 631
      enddo
 631  nk=kk
      nyear=nyrs(nk)
      ns=nstrt(nk)-1

c	check for previous days wetness state
      if(nlon.gt.0)call check_lon(igen,sumlon,lon,nlon)
      
      if(nlon.gt.0)then
      do j1=1,nlon
      nn=lon(2,j1)
      prevl(j1)=sumlon(j1,nn+1)
      enddo
      endif
c	
c	find conditional probability of rf occurrence (based on prev days rf)
c	 and modify the probability based on prev days wetness state and value of atm variables

      jfx=1
      if(llag.eq.0)pw=pro(j,l1,nk)
      if(llag.eq.1)then
      if(ip1.lt.rain)pw=pro1(j,l1,nk)
      if(ip1.ge.rain)pw=pro2(j,l1,nk)
      if(ip1.lt.rain)jfx=1
      if(ip1.ge.rain)jfx=2
      endif
      
      j2=j
      l2=l1
      
      
      if(nlon.gt.0)then
      ind=0
      do j1=1,nlon
      nn=lon(2,j1)
      aa=sumlon(j1,nn+1)
      prevsum(j1)=aa
      prevm(j1)=sumlon(j1,nn+1)
      enddo
      endif
      
      if(nlon.gt.0)then
      do j1=1,nlon
      prevlong(j1)=prevsum(j1)/sdl(j1,j,l1,nk)
      enddo
      endif
      
      anlon=float(nlon)
      
      if(nlon.gt.0)call aloncov(pavl,ppavl,pcovl,ppcovl,detl,ddetl,
     1             ip1,llag,rain,nlon,nk,j,l1)
      
      pn=1.0
      pd=1.0
      pp=pw
      
      if(detl.eq.0.or.ddetl.eq.0)goto 503 
      
           if(nlon.eq.0)goto 503
      call multiply(prevlong,pavl,pcovl,nlon,nvmax,pd)
      call multiply(prevlong,ppavl,ppcovl,nlon,nvmax,pn)
      if(ddet.gt.1.0e-8.and.det.gt.1.0e-8)then
      pn=(1./(2.*pi)**(0.5*anlon)) * (1./sqrt(ddet)) * pn
      pd=(1./(2.*pi)**(0.5*anlon)) * (1./sqrt(det)) * pd
      endif
      
      do jj=1,nlon
      aa=prevm(jj)
      if(aa.lt.(1.0/1.5)*plmin(j,l1,jj,nk).or.aa.gt.1.5*
     1               plmax(j,l1,jj,nk))then
      pd=1.0
      pn=1.0
      endif
      enddo
c	write(*,*)nk,pp,pn,pd
c	pause
c	
      pp=pn*pp/(pp*pn+(1-pp)*pd)
 503  continue

      pw=pp
 22   format(4i4,25f11.7)
 23   format(6f11.5)
 27   format(3i4,33x,8f11.4)
 29   format(4i4,5f8.3,i5)
 31   format(f16.12,8x,10f8.3)
 32   continue

c	generate rf occurrences
      vv=ran1(iseed)
      igen=0
      if(vv.le.pw)igen=1
      
      if(iamt.gt.0)then
      ic=ip
      jc=jp
      lc=lp
      inext=igen
      else
      ic=i
      jc=j
      lc=l1
      endif
      
      if(i.eq.1.and.j.eq.1.and.iamt.gt.0.and.l1.eq.1)then
      icur=igen
      p1=0.0
      if(icur.ge.rain)p1=rain
      goto 902
      endif					   
      
      if(iamt.gt.0)then
           if(icur.gt.0)then
      call rf_amt_gen(icur,gen,nk,ic,jc,lc,p1,inext)
      gen=gen*tavrf/avrf(nk)
      if(gen.lt.rain)gen=rain
      else
      gen=0.0
      endif
      endif
      
      if(iamt.eq.0)gen=float(igen)
      rf_gen(1,ic,jc,lc)=gen
      write(2,901)ic+nsg,jc,lc,gen
 901  format(3i5,50(1x,f6.1))

      p1=gen
      ip3=ip2
      ip2=ip1
      ip1=float(igen)
      icur=igen
      
      wgen(is)=wgen(is)+float(igen)
      sumg(is)=sumg(is)+float(igen)
      wgenr(is)=wgenr(is)+gen
      sumgr(is)=sumgr(is)+gen

 902  continue
      ip=i
      jp=j
      lp=l1
           enddo
           enddo
      
      sumho=0.0
      sumhr=0.0
      do is=1,isn
      sumho(is)=sumho(is)+whist(is,nstn+1)
      sumhr(is)=sumhr(is)+whamt(is,nstn+1)
      enddo
      if(it.le.3)write(*,703)i,(sumho(is),sumg(is),sumhr(is),sumgr(is),
     1           is=1,isn)
           enddo
      if(iamt.eq.1)then
       write(2,901)ng+nsg,nout,nday(nout),gen
      rf_gen(1,ng,nout,nday(nout))=gen
      endif
      
      
      ind=0
      sumho=0.0
      sumhr=0.0
      sumg=0.0
      sumgr=0.0
      
      do is=1,isn
      wgen(is)=wgen(is)/float(ng)
      wgenr(is)=wgenr(is)/float(ng)
      
      if(wgen(is).le.(1./1.5)*whist(is,nstn+1).or.wgen(is).ge.1.5*
     1          whist(is,nstn+1))ind=ind+1
      enddo
c	write(*,702)it,ind,(whist(is,nstn+1),
c	1                       wgen(is),whamt(is,nstn+1),wgenr(is),
c     2                       is=1,isn)

      write(2,*)it
 702  format(2i4,50f7.1)
 703  format('  year ',i4,50f7.1)
c				
      call stat_g(it,rf_gen,1,nsg,ng,rain,nout,nday,isn)
        enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine wetdays (wet,rfwet)
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      real*4 wet(5,nstnmax),rfwet(5,nstnmax)
      real*4 sum(5),rfsum(5)
      
      do is=1,isn
      do nk=1,nstn+1
      wet(is,nk)=0.0
      rfwet(is,nk)=0.0
      enddo
      enddo
      
      do nk=1,nstn
      nyear=nyrs(nk)
      ns=nstrt(nk)
      do i=1,nyear
             do j=1,nout
             nd=nday(j)
      	is=iseas(j,isn)
             if(j.eq.2)call daycount(ns,i,nd)
             do l=1,nd
      rfwet(is,nk)=rfwet(is,nk)+rf(nk,i,j,l)
      if(rf(nk,i,j,l).ge.rain)wet(is,nk)=wet(is,nk)+1.0
      enddo
      enddo
      enddo
      enddo
      
      sum1=0.0
      rfsum=0.0
      do is=1,isn
      do nk=1,nstn
      nyear=nyrs(nk)
      wet(is,nk)=wet(is,nk)/float(nyear)
      rfwet(is,nk)=rfwet(is,nk)/float(nyear)
      sum(is)=sum(is)+wet(is,nk)
      rfsum(is)=rfsum(is)+rfwet(is,nk)
      enddo
      wet(is,nstn+1)=sum(is)/float(nstn)
      rfwet(is,nstn+1)=rfsum(is)/float(nstn)
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rf_amt_gen(icur,gen,nk,ic,jc,lc,p1,inext)
c      Implicit double precision (A-H,O-Z) 
      include  'para.inc'
      include  'data.inc'
      include  'rf_amt.inc'
      real*4 yt(nearmax1),xt(nearmax1,nvmax),gen
      real*4 xm,xn,xcond(nvmax)
      real*4 p1,w(nvmax)
      integer icur,inext
      
      gen=0.0
      if(icur.eq.0)return

c	generate rainfall amounts
      jk=lc
      if(jc.gt.1)then
      do j=1,jc-1
      jk=jk+nday(j)
      enddo
      endif
      
      icp=0			
c	generate two vectors of normally distributed random numbers
 402  nst=1

      vv=ran1(iseed)
      xr=gasdev(iseed)
      
      icc=0

c	form vector of conditioning variables for the given day
      jj=0
      nv=0
      if(p1.ge.rain)then
      nv=nv+1
      xcond(nv)=p1
      endif


c	 find rainfall class
      if(p1.ge.rain.and.inext.gt.0)icg=1
      if(p1.ge.rain.and.inext.eq.0)icg=2
      if(p1.lt.rain.and.inext.gt.0)icg=3
      if(p1.lt.rain.and.inext.eq.0)icg=4



 212  yt=0.0
      xt=0.0
      ii=0
      do i=1,10000000
      if(rfcur(i,jc,lc,nk).lt.-100.0)goto 199
      if(rfprv(i,jc,lc,nk).ge.rain.and.irfnxt(i,jc,lc,nk).gt.0)ih=1
      if(rfprv(i,jc,lc,nk).ge.rain.and.irfnxt(i,jc,lc,nk).eq.0)ih=2
      if(rfprv(i,jc,lc,nk).lt.rain.and.irfnxt(i,jc,lc,nk).gt.0)ih=3
      if(rfprv(i,jc,lc,nk).lt.rain.and.irfnxt(i,jc,lc,nk).eq.0)ih=4
c	  if(icg.eq.1.and.irfnxt(i,jc,lc,nk).gt.0)ih=1
c	  if(icg.eq.3.and.irfnxt(i,jc,lc,nk).gt.0)ih=3
c	  if(icg.eq.2.and.irfnxt(i,jc,lc,nk).eq.0)ih=2
c	  if(icg.eq.4.and.irfnxt(i,jc,lc,nk).eq.0)ih=4

      if(icg.ne.ih)goto 198
      ii=ii+1
      yt(ii)=rfcur(i,jc,lc,nk)
      
      nvv=0
      if(p1.ge.rain)then
      nvv=nvv+1
      xt(ii,nvv)=rfprv(i,jc,lc,nk)
      endif
 198  continue
      enddo
 199  continue
 211  nn=ii


      if(nn.gt.nearmax1)write(*,*)' check dimension nearmax1 line 4357'
      if(nn.le.0)then
      gen=rain
      goto 10
      endif
c	write(*,*)nn,yt(1),xt(1,1),xt(1,2)
c	

      ind=0
      if(nvv.lt.1.or.nn.lt.15)ind=1
      w=1.0
      rn=rain


      call psimmain(yt,nn,xt,nvv,xcond,ys,iseed,vv,xr,rn,w,
     1              nearmax1,nvmax,ind)
c	ys=ys*sd1
      if(ys.gt.600.0)write(*,*)'Daily rainfall at station gt 600 mm'
     1               ,nk,ys
      
      if(ys.gt.600.0)ys=600.0
      if(ys.lt.rain)ys=rain
      
      gen=ys
 10   continue
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rank_r3(x,nstn,isn,nsim,nstnmax,isnmax,itrmax)
      real*4 x(itrmax,isnmax,nstnmax,nstnmax)
      do nk1=1,nstn-1
      do nk2=nk1+1,nstn
      do is=1,isn
      do i=1,nsim-1
      do j=i+1,nsim
      if(x(j,is,nk1,nk2).gt.x(i,is,nk1,nk2))then
      a=x(i,is,nk1,nk2)
      x(i,is,nk1,nk2)=x(j,is,nk1,nk2)
      x(j,is,nk1,nk2)=a
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine rank1_stat(xx,nsim,ny,np,nstn)
      include  'para.inc'
      dimension xx(itrmax,nstnmax,15,nyrmax)
      
      do kk=1,nstn
      do k=1,np
      do nk=1,ny
      do i=1,nsim-1
      do j=i+1,nsim
      if(xx(j,kk,k,nk).gt.xx(i,kk,k,nk))then
      a=xx(i,kk,k,nk)
      xx(i,kk,k,nk)=xx(j,kk,k,nk)
      xx(j,kk,k,nk)=a
      endif
      enddo
      enddo
      enddo
      enddo
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine result(isn,nstn,nsim,nout,nyear,ng)
c	subroutine to write the computed statistics in a file
      include 'para.inc'
      include 'stat.inc'
      
      nh=7
      open(3,file='stat.out')
      ns=isn+1
      write(3,*)'Mean of Seasonal & Annual wet days'
      it=0
      write(3,10)it,((roavsh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((roavsg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'Mean of Seasonal & Annual rf amount'
      it=0
      write(3,10)it,((rfavsh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rfavsg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'Mean of Seasonal & Annual rf amount per wet day'
      it=0
      write(3,10)it,((rwavsh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rwavsg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'SD of Seasonal & Annual wet days'
      it=0
      write(3,10)it,((rosdsh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rosdsg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'SD of Seasonal and Annual rf amount'
      it=0
      write(3,10)it,((rfsdsh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rfsdsg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'SD of Seasonal and Annual rf amount per wet day'
      it=0
      write(3,10)it,((rwsdsh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rwsdsg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      ns=nout
      write(3,*)'Mean of Monthly wet days'
      it=0
      write(3,10)it,((roavmh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((roavmg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'Mean of Monthly rf amount'
      it=0
      write(3,10)it,((rfavmh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rfavmg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'SD of Monthly wet days'
      it=0
      write(3,10)it,((rosdmh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rosdmg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'SD of Monthly rf amount'
      it=0
      write(3,10)it,((rfsdmh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rfsdmg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'Maxm wet & dry spell lengths (days)& 5th% rf amt (mm)'
      it=0
      write(3,10)it,(wmaxh(nk),dmaxh(nk),hperh(nk,isn+1),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,(wmaxg(it,nk),dmaxg(it,nk),hperg(it,nk,isn+1),
     1           nk=1,nstn)
      enddo
      
      ns=4
      write(3,*)'Average occurrence of wet pells of varying lengths'
      write(3,*)'      1-day     2-4 days     5-7 days     >7 days'
      it=0
      write(3,10)it,((rfspnh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rfspng(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'Total rainfall in wet spells of varying lengths'
      write(3,*)'      1-day     2-4 days     5-7 days     >7 days'
      it=0
      write(3,10)it,((rfsplh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rfsplg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      write(3,*)'Average occurrence of dry pells of varying lengths'
      write(3,*)'      1-day     2-9 days    9-18 days    >18 days'
      it=0
      write(3,10)it,((rfspndh(nk,is),is=1,ns),nk=1,nstn)
      do it=1,nsim
      write(3,10)it,((rfspndg(it,nk,is),is=1,ns),nk=1,nstn)
      enddo
      call rank_r3(crsodg,nstn,isn+1,nsim,nstnmax,isnmax,itrmax)
      call rank_r3(crsrdwg,nstn,isn+1,nsim,nstnmax,isnmax,itrmax)
      call rank_r3(crsrdag,nstn,isn+1,nsim,nstnmax,isnmax,itrmax)
      call rank_r3(crsog,nstn,isn+1,nsim,nstnmax,isnmax,itrmax)
      call rank_r3(crsrg,nstn,isn+1,nsim,nstnmax,isnmax,itrmax)
      
      n=isn+1
c	find 5% and 95% levels
      i5=int(float(nsim)*5.0/100.0+0.5)
      i95=int(float(nsim)*95.0/100.0+0.5)
      i50=int(float(nsim)*50.0/100.0+0.5)
      if(i5.lt.1)i5=1
      if(i95.gt.nsim)i95=nsim
      goto 31
      
      write(3,*)'Daily log-odds ratio on seasonal and annual basis'
      write(3,19)
 19   format('  STN1  STN2   obs     5%     50%     95%  and repeat')

      do i=1,nstn-1
      do j=i+1,nstn
      write(3,40)i,j,(crsodh(k,i,j),crsodg(i5,k,i,j),
     1      crsodg(i50,k,i,j),crsodg(i95,k,i,j),k=1,n)
      enddo
      enddo
      write(3,*)'Daily crs-corrl of rf-wet days on seas & annual basis'
      do i=1,nstn-1
      do j=i+1,nstn
      write(3,40)i,j,(crsrdwh(k,i,j),crsrdwg(i5,k,i,j),
     1      crsrdwg(i50,k,i,j),crsrdwg(i95,k,i,j),k=1,n)
      enddo
      enddo
      write(3,*)'Daily crs-corrl of rf-all days on seas & annual basis'
      do i=1,nstn-1
      do j=i+1,nstn
      write(3,40)i,j,(crsrdah(k,i,j),crsrdag(i5,k,i,j),
     1      crsrdag(i50,k,i,j),crsrdag(i95,k,i,j),k=1,n)
      enddo
      enddo
      write(3,*)'Cross-corrl of number wet days in season and year'
      do i=1,nstn-1
      do j=i+1,nstn
      write(3,40)i,j,(crsoh(k,i,j),crsog(i5,k,i,j),
     1      crsog(i50,k,i,j),crsog(i95,k,i,j),k=1,n)
      enddo
      enddo
      write(3,*)'Cross-corrl of seasonal & annual rainfall'
      do i=1,nstn-1
      do j=i+1,nstn
      write(3,40)i,j,(crsrh(k,i,j),crsrg(i5,k,i,j),
     1      crsrg(i50,k,i,j),crsrg(i95,k,i,j),k=1,n)
      enddo
      enddo

 31   call rank1_stat(spellg,nsim,ng,nh,nstn)

      do isw=1,nstn
      do jj=1,5
      if(jj.eq.1)write(3,*)'Maximum rainfall amount in a year'
      if(jj.eq.2)write(3,*)'Number of wet days in a year'
      if(jj.eq.3)write(3,*)'Average rainfall amount in a year'
      if(jj.eq.4)write(3,*)'Maximum wet spell length in a year'
      if(jj.eq.5)write(3,*)'Average number of wet spells>1 in a year'
      write(3,21)
 21   format('  STN    year  ppo    ppg   obs      5%     50%     95% ')

      nn=max(nyear,ng)

      do nk=1,nn
c	plotting position
      aph=(float(nk)-0.44)/(nyear+0.12)
      apg=(float(nk)-0.44)/(ng+0.12)
      write(3,20)isw,nk,aph,apg,spellh(isw,jj,nk),
     1 spellg(i5,isw,jj,nk),spellg(i50,isw,jj,nk),spellg(i95,isw,jj,nk)
      enddo
      enddo
      do jj=1,2
      if(jj.eq.1)write(3,*)'Maximum dry spell length in a year'
      if(jj.eq.2)write(3,*)'Average number of dry spells>1 in a year'
      do nk=1,nn
c	plotting position
      aph=(float(nk)-0.44)/(nyear+0.12)
      apg=(float(nk)-0.44)/(ng+0.12)
      write(3,20)isw,nk,aph,apg,spellh(isw,jj+5,nk),
     1           spellg(i5,isw,jj+5,nk),spellg(i50,isw,jj+5,nk),
     2           spellg(i95,isw,jj+5,nk)
      enddo
      enddo
      enddo


 10   format(i5,1200f7.1)
 20   format(i5,i7,2f7.2,1200f7.1)
 30   format(i5,1200f7.2)
 40   format(i5,i7,1200f7.2)
      close(unit=3)
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      function iseas(j,isn)
      integer iseas,j,isn
      if(isn.eq.4)then
      if(j.ge.3.and.j.le.5)iseas=1
      if(j.ge.6.and.j.le.8)iseas=2
      if(j.ge.9.and.j.le.11)iseas=3
      if(j.le.2.or.j.eq.12)iseas=4
      return
      elseif(isn.eq.2)then
      iseas=2
      if(j.ge.6.and.j.le.10)iseas=1
      return
      elseif(isn.eq.3)then
      if(j.ge.1.and.j.le.5)iseas=1
      if(j.ge.6.and.j.le.9)iseas=2
      if(j.ge.10.and.j.le.12)iseas=3
      return
      endif
      return
      end
	  
      FUNCTION ran1(idum)
c     calculates a random number in the range 0-1       
c      Implicit double precision (A-H,O-Z) 
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV  
      REAL ran1,AM,EPS,RNMX  
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,  
     1   NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)  
      INTEGER j,k,iv(NTAB),iy  
      SAVE iv,iy  
      DATA iv /NTAB*0/, iy /0/  
      if (idum.le.0.or.iy.eq.0) then  
      idum=max(-idum,1)  
      do 11 j=NTAB+8,1,-1  
      k=idum/IQ  
      idum=IA*(idum-k*IQ)-IR*k  
      if (idum.lt.0) idum=idum+IM  
      if (j.le.NTAB) iv(j)=idum  
 11   continue  
      iy=iv(1)  
      endif  
      k=idum/IQ  
      idum=IA*(idum-k*IQ)-IR*k  
      if (idum.lt.0) idum=idum+IM  
      j=1+iy/NDIV  
      iy=iv(j)  
      iv(j)=idum  
      ran1=min(AM*iy,RNMX)  
      return  
      END

      real*4 FUNCTION gasdev(idum)
c	function to generate normal random deviates
c      Implicit double precision (A-H,O-Z) 
      INTEGER idum
c      real*4 gasdev
c    USES ran1
      INTEGER iset
      real*4 fac,gset,rsq,v1,v2,ran1
      SAVE iset,gset
      DATA iset/0/
      if(iset.lt.1)then
 1    v1=2.*ran1(idum)-1.
      v2=2.*ran1(idum)-1.
      rsq=v1**2+v2**2
      if(rsq.ge.1..or.rsq.eq.0.)goto 1
      fac=sqrt(-2.*log(rsq)/rsq)
      gset=v1*fac
      gasdev=v2*fac
      iset=1
      else
      gasdev=gset
      iset=0
      endif
      return
      END