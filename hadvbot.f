      subroutine hadvbot(nn,dt,dx,con,vel,mscl,flxarr,flux1,
     &                              flux2,saflux,fpc,fmc,fc1,fc2)
c  
c-----CAMx v4.02 030709
c  
c     HADVBOT performs horizontal advection using the positive-definite  
c     mass conservative non-monotonic scheme of Bott (1989). Fourth order  
c     area preserving polynomials are used for the interior computational 
c     cells. First and second order polynomials are applied at the boundaries. 
c  
c     The following definitions are used:  
c  
c              |-----------> Positive direction  
c  
c     |Boundary|<-----------------Domain----------------->|Boundary| 
c  
c     | CON(1) | CON(2) |  ...  | CON(I) |  ...  |CON(N-1)| CON(N) |  
c  
c     VEL(1)-->|     VEL(I-1)-->|        |-->VEL(I)       |-->VEL(N-1)  
c  
c      FP(1)-->|      FP(I-1)-->|        |-->FP(I)        |-->FP(N-1)  
c  
c      FM(2)<--|        FM(I)<--|        |<--FM(I+1)      |<--FM(N)  
c             
c                            -->|   DS   |<-- 
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c  
c     Modifications:  
c        Spr 2000  Fluxes extracted for advection of DDM sens coeffs
c        5/17/00   small modification to flux1,2 to improve mass accounting
c  
c     Input arguments:  
c        nn                  Number of cells  
c        dt                  Time step (s)  
c        dx                  Length of cell (m)  
c        con                 Concentration*area vector (umol/m) 
c        vel                 Wind speed vector (m/s)  
c        mscl                Map-scale factor (squared) at cell centroid
c  
c     Output arguments:  
c        con                 Concentration*area vector (umol/m)  
c        flxarr              Interfacial mass flux (umol/s)
c        flux1-2             Boundary fluxes (umol/s)  
c                            (1=west/south, 2=east/north)
c        saflux              Interfacial mass flux times area
c                            (used for tracer transport)
c        fpc, fmc            Positive and negative concentration fluxes before
c                            the flux limiting renormalization step (umol/m)
c                            Used in calculating sensitivities.
c        fc1-2               Conc*area vector change from flux (umol/m)
c                            (1=west/south, 2=east/north)
c                            (used for Process Analysis)
c  
c     Routines called:  
c        none  
c  
c     Called by:  
c        XYADVEC
c        XYADVECS  
c        ZRATES
c        TIMADV
c 
      include "camx.prm"
c
      real con(nn),vel(nn),flxarr(nn),mscl(nn), saflux(nn)
      real fpc(nn), fmc(nn)
      real cn(MX1D),fm(MX1D),fp(MX1D)
      real*8 flux1,flux2
c
c======================== Process Analysis Begin ====================================
c
      real fc1(mx1d), fc2(mx1d)
c
c========================= Process Analysis End =====================================
c
      data small / 1.e-20/
c             
c-----Entry point 
c 
c-----Advect (concentration) X (interfacial area) 
c     Set all fluxes to zero, calculate Courant number vector
c
      do i = 1,nn
        fm(i) = 0.
        fp(i) = 0.
        cn(i) = vel(i)*dt/dx
c
c======================== Process Analysis Begin ====================================
c
        fc1(i) = 0.
        fc2(i) = 0.
c
c========================= Process Analysis End =====================================
c
      enddo
c             
c-----First order polynomials for boundary inflow fluxes 
c             
c-----Left boundary 
c
      cp = amax1(0.,cn(1))
      fp(1) = amin1(con(1),
     &              cp*(con(1)+0.5*(1.-cp)*(con(2)-con(1))))
c             
c-----Right boundary 
c
      cm = -amin1(0.,cn(nn-1))
      fm(nn) = amin1(con(nn),
     &               cm*(con(nn)-0.5*(1.-cm)*(con(nn)-con(nn-1))))
c             
c-----Second order polynomials for fluxes from the first and last  
c     computational cells 
c             
c-----Leftmost computational cell 
c
      a0 = (-con(1) + 26.*con(2) - con(3))/24.
      a1 = (-con(1)              + con(3))/16.
      a2 = ( con(1) -  2.*con(2) + con(3))/48.
      cm = -amin1(0.,cn(1))
      x1 = 1. - 2.*cm
      x2 = x1*x1
      fm(2) = amax1(0., a0*cm-a1*(1.-x2)+a2*(1.-x1*x2))
      cp = amax1(0.,cn(2) )
      x1 = 1. - 2.*cp
      x2 = x1*x1
      fp(2) = amax1(0., a0*cp+a1*(1.-x2)+a2*(1.-x1*x2))
c
c======================== DDM Begin =======================
c
c
c-----Save fluxes for calculating sensitivities
c
      fpc(1) = fp(1)
      fmc(nn) = fm(nn)
      fmc(2) = fm(2)
      fpc(2) = fp(2)
c
c======================== DDM End =======================
c
c  
c-----Flux limiting renormalization
c
      wt =  con(2)/amax1(con(2), fm(2)+fp(2)+small)
      fm(2) = fm(2)*wt
      fp(2) = fp(2)*wt
c             
c-----Rightmost computational cell 
c
      a0 = (-con(nn-2) + 26.*con(nn-1) - con(nn))/24.
      a1 = (-con(nn-2)                 + con(nn))/16.
      a2 = ( con(nn-2) -  2.*con(nn-1) + con(nn))/48.
      cm = -amin1(0.,cn(nn-2))
      x1 = 1. - 2.*cm
      x2 = x1*x1
      fm(nn-1) = amax1(0., a0*cm-a1*(1.-x2)+a2*(1.-x1*x2))
      cp = amax1(0.,cn(nn-1) )
      x1 = 1. - 2.*cp
      x2 = x1*x1
      fp(nn-1) = amax1(0., a0*cp+a1*(1.-x2)+a2*(1.-x1*x2))
c
c======================== DDM Begin =======================
c
c
c-----Save fluxes for calculating sensitivities
c
      fmc(nn-1) = fm(nn-1)
      fpc(nn-1) = fp(nn-1)
c
c======================== DDM Begin =======================
c
c-----Flux limiting renormalization
c
      wt = con(nn-1)/amax1(con(nn-1), fm(nn-1)+fp(nn-1)+small)
      fm(nn-1) = fm(nn-1)*wt
      fp(nn-1) = fp(nn-1)*wt
c             
c-----Fourth order polynomials for fluxes from the interior cells 
c
      do i = 3,nn-2
        a0 = (   9.*(con(i+2) + con(i-2)) - 
     &         116.*(con(i+1) + con(i-1)) + 2134.*con(i))/1920.
        a1 = (  -5.*(con(i+2) - con(i-2)) +  
     &          34.*(con(i+1) - con(i-1)))/384.
        a2 = (      -con(i+2) + 
     &          12.*(con(i+1) + con(i-1)) - 
     &           22.*con(i)   - con(i-2))/384.
        a3 = (       con(i+2) - 
     &           2.*(con(i+1) - con(i-1)) - con(i-2))/768.
        a4 = (       con(i+2) - 
     &           4.*(con(i+1) + con(i-1)) + 
     &            6.*con(i)   + con(i-2))/3840.
        cm = -amin1(0.,cn(i-1))
        x1 = 1. - 2.*cm
        x2 = x1*x1
        x3 = x1*x2
        fm(i) = amax1(0., a0*cm - a1*(1. - x2)    +
     &                            a2*(1. - x3)    -
     &                            a3*(1. - x1*x3) +
     &                            a4*(1. - x2*x3))
        cp = amax1(0.,cn(i))
        x1 = 1. - 2.*cp
        x2 = x1*x1
        x3 = x1*x2
        fp(i) = amax1(0., a0*cp + a1*(1. - x2)    +
     &                            a2*(1. - x3)    +
     &                            a3*(1. - x1*x3) +
     &                            a4*(1. - x2*x3))
c
c======================== DDM Begin =======================
c
c
c------Save fluxes for calculating sensitivities
c
        fmc(i) = fm(i)
        fpc(i) = fp(i)
c
c======================== DDM End =======================
c
c
c-----Flux limiting renormalization
c
        wt = con(i)/amax1(con(i), fm(i)+fp(i)+small)
        fm(i) = fm(i)*wt
        fp(i) = fp(i)*wt
      enddo
c 
c-----Update concentrations in computational cells 
c
      flxarr(1) = (fp(1) - fm(2))*dx/dt
      saflux(1) = flxarr(1)*dt/dx
      do i = 2,nn-1
         flxarr(i) = (fp(i) - fm(i+1))*dx/dt
         con(i) = con(i) - mscl(i)*(flxarr(i) - flxarr(i-1))*dt/dx
         saflux(i) = flxarr(i)*dt/dx
         con(i) = amax1(con(i),1.0e-20)
c
c======================== Process Analysis Begin ====================================
c
         fc1(i) =   mscl(i)*flxarr(i-1)*dt/dx
         fc2(i) = - mscl(i)*flxarr(i)*dt/dx
c
c========================= Process Analysis End =====================================
c
      enddo
      flux1 = mscl(2)*flxarr(1)
      flux2 = mscl(nn-1)*flxarr(nn-1)
c
      return
      end
