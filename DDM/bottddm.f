      subroutine bottddm(nn,dt,dx,sen,nsen,con0,con,fpc,fmc,vel,mscl)
c  
c-----CAMx v4.02 030709
c  
c     HADVBOTS performs horizontal advection of the sensitivity coefficients.
c     The advection algorithm is the positive-definite, mass conservative,
c     non-monotonic scheme of Bott (1989). The sensitivity coefficients are
c     calculated by the decoupled direct method. It is assumed that none of 
c     the sensitivity coefficients represents the sensitivity with respect to
c     a parameter involved in the advection step.  E.g., none of the
c     parameters is the wind velocity. Fourth order area preserving
c     polynomials are used for the interior computational cells. First and
c     second order polynomials are applied at the boundaries.  This
c     subroutine is and must be consistent with HADVBOT.
c  
c     The following definitions are used:  
c  
c              |-----------> Positive direction  
c  
c     | Boundary|<------------------Domain------------------->| Boundary| 
c  
c     |  SEN(1) |  SEN(2) |  ...  |  SEN(I) |  ...  | SEN(N-1)|  SEN(N) |  
c     | CON0(1) | CON0(2) |  ...  | CON0(I) |  ...  |CON0(N-1)| CON0(N) |  
c     |  CON(1) |  CON(2) |  ...  |  CON(I) |  ...  | CON(N-1)|  CON(N) |  
c  
c      VEL(1)-->|      VEL(I-1)-->|         |-->VEL(I)        |-->VEL(N-1)  
c  
c       FP(1)-->|       FP(I-1)-->|         |-->FP(I)         |-->FP(N-1)  
c      FPC(1)-->|      FPC(I-1)-->|         |-->FPC(I)        |-->FPC(N-1) 
c  
c       FM(2)<--|         FM(I)<--|         |<--FM(I+1)       |<--FM(N)  
c      FMC(2)<--|        FMC(I)<--|         |<--FMC(I+1)      |<--FMC(N) 
c             
c                              -->|   DS    |<-- 
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c  
c     Modifications:  
c        None 
c  
c     Input arguments:  
c        nn                  Number of cells  
c        dt                  Time step (s)  
c        dx                  Length of cell (m)  
c        sen                 Sensitivity coefficient*area vector 
c                              (umol/m/parameter unit) 
c        nsen                Number of sensitivity coefficients
c        con0                Concentration*area vector (umol/m) at start of
c                            timestep
c        con                 Concentration*area vector (umol/m) at end of
c                            timestep
c        fpc, fmc            Positive and negative concentration fluxes before
c                            the flux limiting renormalization step (umol/m)
c        vel                 Wind speed vector (m/s)  
c        mscl                Map-scale factor (squared) at cell centroid
c  
c     Output arguments:  
c        sen                 Sensitivity coefficient*area vector 
c                              (umol/m/parameter unit)  
c  
c     Routines called:  
c        none  
c  
c     Called by:  
c        XYADVEC  
c 
      include "camx.prm"
      include "tracer.com"
c
      real con0(nn),con(nn),fpc(nn),fmc(nn),vel(nn),mscl(nn)
      real sen(MX1D,MXTRSP)
      real a0(MXTRSP),a1(MXTRSP),a2(MXTRSP),a3(MXTRSP),a4(MXTRSP)
      real cn(MX1D) 
      real fm(MX1D,MXTRSP),fp(MX1D,MXTRSP)
c
      data small / 1.e-20/
c             
c-----Entry point 
c 
c-----Advect (senstivity coefficient) X (interfacial area) 
c     Calculate Courant number vector
c 
      do i = 1,nn
        cn(i) = vel(i)*dt/dx
      enddo

c     Set all sensitivity fluxes to zero.
c
      do j = 1,nsen
        do i = 1,nn
          fm(i,j) = 0.
          fp(i,j) = 0.
        enddo
      enddo
c             
c-----First order polynomials for boundary inflow fluxes 
c             
c-----Left boundary 
c
      cp = amax1(0.,cn(1))
      if (fpc(1).lt.con0(1)) then
        do j=1,nsen
          fp(1,j) = cp*(sen(1,j)+0.5*(1.-cp)*(sen(2,j)-sen(1,j)))
        enddo
      else
        do j=1,nsen
          fp(1,j) = sen(1,j)
        enddo
      endif
c             
c-----Right boundary 
c  
      cm = -amin1(0.,cn(nn-1))
      if (fmc(nn).lt.con0(nn)) then
        do j=1,nsen
          fm(nn,j) = cm*(sen(nn,j)-0.5*(1.-cm)*(sen(nn,j)-sen(nn-1,j)))
        enddo
      else
        do j=1,nsen
          fm(nn,j) = sen(nn,j)
        enddo
      endif
c             
c-----Second order polynomials for fluxes from the first and last  
c     computational cells 
c             
c-----Leftmost computational cell 
c
      do j=1,nsen
        a0(j) = (-sen(1,j) + 26.*sen(2,j) - sen(3,j))/24.
        a1(j) = (-sen(1,j)                + sen(3,j))/16.
        a2(j) = ( sen(1,j) -  2.*sen(2,j) + sen(3,j))/48.
      enddo
      cm = -amin1(0.,cn(1))
      x1 = 1. - 2.*cm
      x2 = x1*x1
      if (fmc(2).gt.0.) then
        do j=1,nsen
          fm(2,j) = a0(j)*cm-a1(j)*(1.-x2)+a2(j)*(1.-x1*x2)
        enddo
      endif
c  
      cp = amax1(0.,cn(2) )
      x1 = 1. - 2.*cp
      x2 = x1*x1
      if (fpc(2).gt.0.) then
        do j=1,nsen
          fp(2,j) = a0(j)*cp+a1(j)*(1.-x2)+a2(j)*(1.-x1*x2)
        enddo
      endif
c
c-----Flux limiting renormalization
c
      if (fmc(2)+fpc(2)+small.gt.con0(2)) then
        denom = fmc(2)+fpc(2)+small
        wt = con0(2)/denom
        wtm = fmc(2)/denom
        wtp = fpc(2)/denom
        do j=1,nsen
          fm(2,j) = fm(2,j)*wt + sen(2,j)*wtm
     &            -(fm(2,j)+fp(2,j))*wt*wtm
          fp(2,j) = fp(2,j)*wt + sen(2,j)*wtp
     &            -(fm(2,j)+fp(2,j))*wt*wtp
        enddo
      endif
c
c-----Rightmost computational cell 
c
c
      do j=1,nsen
        a0(j) = (-sen(nn-2,j) + 26.*sen(nn-1,j) - sen(nn,j))/24.
        a1(j) = (-sen(nn-2,j)                   + sen(nn,j))/16.
        a2(j) = ( sen(nn-2,j) -  2.*sen(nn-1,j) + sen(nn,j))/48.
      enddo
      cm = -amin1(0.,cn(nn-2))
      x1 = 1. - 2.*cm
      x2 = x1*x1
      if (fmc(nn-1).gt.0.) then
        do j=1,nsen
          fm(nn-1,j) = a0(j)*cm-a1(j)*(1.-x2)+a2(j)*(1.-x1*x2)
        enddo
      endif
c  
      cp = amax1(0.,cn(nn-1) )
      x1 = 1. - 2.*cp
      x2 = x1*x1
      if (fpc(nn-1).gt.0.) then
        do j=1,nsen
          fp(nn-1,j) = a0(j)*cp+a1(j)*(1.-x2)+a2(j)*(1.-x1*x2)
        enddo
      endif
c
c-----Flux limiting renormalization
c
      if (fmc(nn-1)+fpc(nn-1)+small.gt.con0(nn-1)) then
        denom = fmc(nn-1)+fpc(nn-1)+small
        wt = con0(nn-1)/denom
        wtm = fmc(nn-1)/denom
        wtp = fpc(nn-1)/denom
        do j=1,nsen
          fm(nn-1,j) = fm(nn-1,j)*wt + sen(nn-1,j)*wtm
     &             -(fm(nn-1,j)+fp(nn-1,j))*wt*wtm
          fp(nn-1,j) = fp(nn-1,j)*wt + sen(nn-1,j)*wtp
     &             -(fm(nn-1,j)+fp(nn-1,j))*wt*wtp
        enddo
      endif
c
c-----Fourth order polynomials for fluxes from the interior cells 
c
      do 10 i = 3,nn-2
        do j=1,nsen
          a0(j) = (  9.*(sen(i+2,j) + sen(i-2,j)) - 
     &             116.*(sen(i+1,j) + sen(i-1,j)) +
     &             2134.*sen(i,j))/1920.
          a1(j) = ( -5.*(sen(i+2,j) - sen(i-2,j)) +  
     &              34.*(sen(i+1,j) - sen(i-1,j)))/384.
          a2(j) = (     -sen(i+2,j) + 
     &              12.*(sen(i+1,j) + sen(i-1,j)) - 
     &               22.*sen(i,j)   - sen(i-2,j))/384.
          a3(j) = (      sen(i+2,j) - 
     &               2.*(sen(i+1,j) - sen(i-1,j)) - sen(i-2,j))/768.
          a4(j) = (      sen(i+2,j) - 
     &               4.*(sen(i+1,j) + sen(i-1,j)) + 
     &                6.*sen(i,j)   + sen(i-2,j))/3840.
        enddo
        cm = -amin1(0.,cn(i-1))
        x1 = 1. - 2.*cm
        x2 = x1*x1
        x3 = x1*x2
        if (fmc(i).gt.0.) then
          do j=1,nsen
            fm(i,j) = a0(j)*cm - a1(j)*(1. - x2)    +
     &                           a2(j)*(1. - x3)    -
     &                           a3(j)*(1. - x1*x3) +
     &                           a4(j)*(1. - x2*x3)
          enddo
        endif
c 
        cp = amax1(0.,cn(i))
        x1 = 1. - 2.*cp
        x2 = x1*x1
        x3 = x1*x2
        if (fpc(i).gt.0.) then
          do j=1, nsen
            fp(i,j) = a0(j)*cp + a1(j)*(1. - x2)    +
     &                           a2(j)*(1. - x3)    +
     &                           a3(j)*(1. - x1*x3) +
     &                           a4(j)*(1. - x2*x3)
          enddo
        endif
c
c-----Flux limiting renormalization
c
        if (fmc(i)+fpc(i)+small.gt.con0(i)) then
          denom = fmc(i)+fpc(i)+small
          wt = con0(i)/denom
          wtm = fmc(i)/denom
          wtp = fpc(i)/denom
          do j=1,nsen
            fm(i,j) = fm(i,j)*wt + sen(i,j)*wtm
     &             -(fm(i,j)+fp(i,j))*wt*wtm
            fp(i,j) = fp(i,j)*wt + sen(i,j)*wtp
     &             -(fm(i,j)+fp(i,j))*wt*wtp
          enddo
        endif
10    continue
c
c-----Update sensitivities in computational cells 
c
      do i = 2,nn-1
c-----Test if negative or small concentration was reset in hadvbot to 1.0e-20.
        if (con(i).gt.1.0e-20) then
          do j=1,nsen
            sen(i,j) = sen(i,j) - mscl(i)*( fp(i,j) - fm(i+1,j) - 
     &                 (fp(i-1,j) - fm(i,j)))
          enddo
        else
          do j=1,nsen
            sen(i,j) = 0.0
          enddo
        endif
      enddo
c
      return
      end

