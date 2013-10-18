      subroutine gresgrow(n,dt,rkv,axiszmx,axisymx,o3amb,dumpmass)
c
c-----CAMx v4.02 030709
c
c     GRESGROW performs the following functions:
c             (1) grows puffs
c             (2) calculates the mass to dump if growth restricted
c             (3) entrains grid mass if growth allowed
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        10/3/02     Modified lateral (sigy) puff growth rate
c
c     Input arguments:
c        n                   puff index
c        dt                  time step (s)
c        rkv                 vertical diffusivity (m2/s)
c        axiszmx             maximum allowable vertical depth (m)
c        axisymx             maximum allowable lateral size (m)
c        o3amb               ambient ozone concentrations (umol/m3)
c
c     Output arguments:
c        dumpmass            NO, NO2, HNO3 mass to dump (umol)
c
c     Subroutines called:
c        none
c
c     Called by:
c        GRESDRIVE
c 
      include "camx.prm"
      include "pigsty.com"
      include "filunit.com"
c
      dimension dumpmass(3)
      data pi/3.1415926/
c
c-----Entry point
c
c-----Check for phantom puff
c
      if (fmspig(n) .le. 0.0) then
        ingrd(n)=0
        return
      endif
c
c-----Calculate horizontal diffusivity from Kv, and grow puff in cross-section
c     only
c
      zscore = 1.5
      rkh = amax1(20.*rkv**(2./5.),rkv) + 8.e-5*(sigy(n)**2)
      volold = xlength(n)*axisy(n)*axisz(n)*pi
      sigzo = sigz(n)
      sigyo = sigy(n)
      sigz(n) = sqrt(sigz(n)*sigz(n) + 2.*rkv*dt)
      sigy(n) = sqrt(sigy(n)*sigy(n) + 2.*rkh*dt)
c
      if (lnewg(n)) then
        axisy(n) = zscore*sigy(n)
        axisz(n) = zscore*sigz(n)
        fmspig(n) = 1. - exp(-zscore*zscore/2.)
        goto 900
      endif
c
c-----Limit puff dimensions by axiszmx and axisymx
c
      fmsold = fmspig(n)
      axisz(n) = zscore*sigz(n)
      if (2.*axisz(n).gt.axiszmx) then
        axisz(n) = axiszmx/2.
        zscore = axisz(n)/sigz(n)
        fmspig(n) = 1. - exp(-zscore*zscore/2.)
      endif
      axisy(n) = zscore*sigy(n)
c
      if (2.*axisy(n).gt.axisymx) then
        axisy(n) = axisymx/2.
        zscore = axisy(n)/sigy(n)
        axisz(n) = zscore*sigz(n)
        fmspig(n) = 1. - exp(-zscore*zscore/2.)
      endif
c
c-----Artificially shrink puffs older than agemax so that 10% of puff mass
c     is released each timestep
c
      if (agepig(n).ge.agemax) then
        fmspig(n) = amax1(0.,fmsold - 0.1)
        zscore = sqrt(-2.*alog(1. - fmspig(n)))
        axisy(n) = zscore*sigy(n)
        axisz(n) = zscore*sigz(n)
      endif
c
c-----Pump out puff mass due to expansion
c
      delms = (fmsold - fmspig(n))/fmsold
      if (delms.gt.1.e-5) then
        do l=1,3
          dumpmass(l) = puffmass(l,n)*delms
          puffmass(l,n) = puffmass(l,n) - dumpmass(l)
        enddo
      else
        do l=1,3
          dumpmass(l) = 0.
        enddo
      endif
c
 900  volnew = xlength(n)*axisy(n)*axisz(n)*pi
      if (volnew.lt.0.0) then
        write(iout,'(//,a)') 'ERROR in GRESGROW:'
        write(iout,*) 'Negative puff volume'
        write(iout,*) 'Puff#,length,axis_x/y/z:'
        write(iout,*) n,xlength(n),axisy(n),axisz(n),axiszmx
        write(iout,*) 'Size/mass parameters:'
        write(iout,*) zscore,fmspig(n),fmsold
        call camxerr()
      endif
c
c-----Update ozone puff mass due to entrainment
c
      if (lnewg(n)) then
        puffmass(4,n) = o3amb*volnew
      elseif (volnew.gt.0.) then
        dvol = volold*(sigy(n)*sigz(n)/(sigyo*sigzo) - 1.)
        if (volnew.le.volold) dvol = amin1(dvol,volnew)
        conpig = puffmass(4,n)/volnew
        if (o3amb.gt.conpig) puffmass(4,n) = puffmass(4,n) + o3amb*dvol
      endif
c
      return
      end
