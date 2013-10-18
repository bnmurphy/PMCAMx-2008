      subroutine aqchem(gas,aerosol,rh,press,temp,
     &       lwc_c,t0,t1,dt,ierr,kchm,height,chtype) 
CKF     &       lwc_c,t0,t1,dt,ierr,kchm,height) 
ckf     &                                  lwc_c,t0,t1,dt,ierr)
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
ckf
ckf
      include 'camx.prm'
c      include 'camx.com'
      include 'grid.com'
      include 'flags.com'
      include 'filunit.com'
c      include "ptemiss.com"
c      include "bndary.com"
ckf         
ckf 

      real*4 gas(ngas_aq), aerosol(nsect,naers)   
      real*4 gasav(ngas_aq), aerosav(nsect,naers)   
      real*4 rh,temp,press,lwc_c,t0,t1,dt,p
ckf   
      real*4 hgt, klay
      integer chtype
      common/hgt/klay, hgt
ckf      

      logical lfrst
      data lfrst /.true./
      save sulfbefppm,sulfbef0
c
c      iseed=43
c      do i=1,64
c        call randnumgen(iseed,randx)
c        call randnumgen(iseed,randx)
c      enddo
c
CKF
       hgt = height
       klay = kchm
CKF
c SET MAIN PROGRAM VARIABLES
c
      istep = 1                      ! print counter
      istart = int(t0)                     ! beginning of simulation in min
      iend =  int(t1)                     ! end of simulation in min
      t0_s=t0
      t1_s=t1
      istart_s=istart
      iend_s=iend
      deltat= dt                    ! operator timestep in min
      p= press                         ! pressure in atm
      if ( lfrst ) then
	lfrst = .false.
        iaq = 1                        ! first call flag
        sulfbef0 = 0.0
        do i=1, nsect
          sulfbef0 = sulfbef0 + aerosol(i,na4)*32./96.
     &      + aerosol(i,nahso5)*32./113.+aerosol(i,nahmsa)*32./111.
        enddo
        sulfbef0 = sulfbef0 + gas(ngso2)*32.*p/(8.314e-5*temp)
        sulfbefppm = 0.0
        do i=1, nsect
          sulfbefppm= sulfbefppm+aerosol(i,na4)
        enddo
        sulfbefppm=sulfbefppm*0.0243113/96.+gas(ngso2)
      endif
c
      do isect=1,nsect
        do isp=1,naers
          aerosav(isect,isp) = aerosol(isect,isp)
        enddo
      enddo
      do i=1,ngas_aq
        gasav(i) = gas(i)
      enddo
c
      so2init=gas(ngso2)
      h2o2init=gas(ngh2o2)
c
c     CALCULATION OF TOTAL SULFUR MASS BEFORE THE CALL
c
      sulfbef = 0.0
      do i=1, nsect
      sulfbef = sulfbef + aerosol(i,na4)*32./96.
     &  + aerosol(i,nahso5)*32./113.+aerosol(i,nahmsa)*32./111.
      enddo
      sulfbef = sulfbef + gas(ngso2)*32.*p/(8.314e-5*temp)
c      write(*,*)'sulfbef',sulfbef
c      pause
 560  call vsrm(gas, aerosol, lwc_c, t0, t1, deltat,
     &          temp, iaq, p, rh, chtype)
ckf     &          temp, iaq, p, rh)
c
c
c     CALCULATION OF TOTAL SULFUR MASS AFTER THE CALL
c
      sulfaf = 0.0
      sulfafppm = 0.0
      do i=1, nsect
      sulfafppm= sulfafppm+aerosol(i,na4)
      sulfaf = sulfaf + aerosol(i,na4)* 32./96.
     &  + aerosol(i,nahso5)*32./113.+aerosol(i,nahmsa)*32./111.
      enddo
      sulfafppm=sulfafppm*0.0243113/96.+gas(ngso2)
      sulfaf = sulfaf+gas(ngso2)*32.*p/(8.314e-5*temp)
c      write(*,*)'sulfaf',sulfaf
cgy
c try an S balance patch
c
      sbal = sulfaf/sulfbef
      sbalppm = sulfafppm/sulfbefppm
      if (deltat.gt.1.0d0) then
      if (sulfbef.gt. 0.1 .and.
     &    (sbal.lt.0.99 .or. sbal.gt.1.01) ) then
cgy
        do isect=1,nsect
          do isp=1,naers
            aerosol(isect,isp) = aerosav(isect,isp)
          enddo
        enddo
        do i=1,ngas_aq
          gas(i) = gasav(i)
        enddo
        t0=t0_s
	t1=t1_s
        istart = istart_s
        iend = iend_s
        deltat = deltat/2.
        goto 560
      endif
      endif
ckf        
ckf      call settling(lwc_c, aerosol, chtype, t0, t1)      
ckf
c
      return
      end
