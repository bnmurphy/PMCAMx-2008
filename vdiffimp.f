      subroutine vdiffimp(nn,dt,vdep,depth,rho,rkv,rr,nparddm,
     &                    fcup,fcdn,ldoipts)
c
c-----CAMx v4.02 030709
c
c     VDIFFIMP performs vertical diffusion of concentrations using 
c     an implicit method, where a tri-diagonal matrix is solved.
c     This version also performs vertical diffusion of sensitivities
c     if DDM is enabled.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        4/17/00   Revised diffusion equations to weight fluxes by density
c
c     Input arguments:
c        nn                number of layers
c        dt                time step (s)
c        vdep              deposition velocity (m/s)
c        depth             layer depth (m)
c        rho               atmospheric density (mb/K)
c        rkv               vertical diffusion coefficient (mb m2/s/K)
c        rr                species concentrations (umol/m3) followed by 
c                          sensitivities (umol/m3/parameter unit).
c        nparddm           number of parameters for which sensitivities
c                          are calculated.  Should be zero if DDM is not
c                          enabled.
c        ldoipts           flag to calculate data for Process Analysis
c
c     Output arguments:
c        rr                species concentrations (umol/m3) followed by
c                          sensitivities (umol/m3/parameter unit)
c        fcup              change in layer concentration due to flux across
c                          upper interface (umol/m3) -- FOR Process Analysis
c        fcdn              change in layer concentration due to flux across
c                          lower interface (umol/m3) -- FOR Process Analysis
c
c     Routines Called:
c        TRDIAG
c
c     Called by:
c        DIFFUS
c
      include "camx.prm"
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
c
      real aa_old(MXLAYA),cc_old(MXLAYA)
      real fcup(MXLAYA),fcdn(MXLAYA)
      logical ldoipts
c
c========================= Process Analysis End =====================================
c
      dimension rr(nn+nn*nparddm),depth(nn),rkv(nn),rho(nn)
      dimension aa(MXLAYA),bb(MXLAYA),cc(MXLAYA)
c
c-----Entry point
c
c-----Lower boundary condition
c
      aa(1) = 0.
      bb(1) = 1. + dt/depth(1)*
     &             (vdep + 2.*rkv(1)/(depth(2)+depth(1))/rho(1))
c
c-----Upper boundary condition
c
      cc(nn) = 0.
      bb(nn) = 1. + dt/depth(nn)*
     &              2.*rkv(nn-1)/(depth(nn-1)+depth(nn))/rho(nn)
c
      do k = 2,nn
        aa(k) = -dt/depth(k)*2.*rkv(k-1)/(depth(k-1)+depth(k))/rho(k-1)
      enddo
      do k = 1,nn-1
        cc(k) = -dt/depth(k)*2.*rkv(k)/(depth(k+1)+depth(k))/rho(k+1)
      enddo
      do k = 2,nn-1
        bb(k) = 1.
     &        + dt/depth(k)*2.*rkv(k-1)/(depth(k-1)+depth(k))/rho(k)
     &        + dt/depth(k)*2.*rkv(k)/(depth(k+1)+depth(k))/rho(k)
      enddo
c
c======================== Process Analysis Begin ====================================
c
      if (ldoipts) then
        do k = 1,nn
          aa_old(k) = 0.
          cc_old(k) = 0.
          if (k.gt.1)  aa_old(k) = aa(k)*rho(k-1)
          if (k.lt.nn) cc_old(k) = cc(k)*rho(k+1)
        enddo
      endif
c 
c========================= Process Analysis End =====================================
c
c
c-----Solve the equations
c
      call trdiag(aa,bb,cc,rr,nn,1+nparddm)
c
c
c======================== Process Analysis Begin ====================================
c
      if (ldoipts) then
        do k = 2,nn-1
           fcup(k) = (rr(k)/rho(k) - rr(k+1)/rho(k+1))*cc_old(k)
           fcdn(k) = (rr(k)/rho(k) - rr(k-1)/rho(k-1))*aa_old(k)
        enddo
c
c-----Lower boundary
c
        fcup(1) = (rr(1)/rho(1) - rr(2)/rho(2))*cc_old(1)
        fcdn(1) = -rr(1)*dt/depth(1)*vdep
c
c-----Upper boundary
c
        fcup(nn) =  0.0
        fcdn(nn) =  (rr(nn)/rho(nn) - rr(nn-1)/rho(nn-1))*aa_old(nn)
      endif
c 
c========================= Process Analysis End =====================================
c
      return
      end
