c  03/09/03: bkoo
c            - modified to eliminate organics (eqparto now deals with organics)
      subroutine step(n,q)

c     This subroutine calculates the equilibrium concetrations and
c     the partial pressures at each time-step for the new generation
c     multicomponent aerosol dynamics model (MADM)
c *************************************************************************
c               WRITTEN BY DR. CHRISTODOULOS PILINIS
c                          April 1997
c *************************************************************************
c

      include 'dynamic.inc'
      real*8 q(ntotal)
c
c     Perform an equilibration for each Section to find
c     the concentrations of all condensible inorganic species.
c     Also partial pressures of species at the particle surface.
c     Then save these concetrations for future use.
c     ps(Pa) in common block in subroutine

c      ng=n*nsp

c     use reverse mode for "safe" thermodynamic range and forward mode
c     for RH conditions that predict partial delequescence (we don't 
c     trust ISORROPIA in this range anyways)
c      if(rh.le.0.6.or.rh.ge.0.80)then
         iprob=1 ! 0--> forward, 1--> reverse
c      else
c         iprob=0
c      endif
c      tpr=1.01325d-1*pres ! pres (bkoo, 06/09/00)
      do k=1,n  ! tmg (10/22/02)
         call equaer(k,q,iprob)
c     if we are near equilibrium use forward mode to prevent the oscillations
c     caused by the reverse method
c         if(abs(1-ps(k,inh3)/q(ng+inh3)/tpr).lt.0.05)then
c            if(abs(1-ps(k,ihno3)/q(ng+ihno3)/tpr).lt.0.05)then
c               if(abs(1-ps(k,ihcl)/q(ng+ihcl)/tpr).lt.0.05)then
c                  iprob=0
c                  call equaer(k,q,iprob)
c                  iprob=1
c                  write(51,*)tcom,' fwd'
c               endif
c            endif
c         endif
      enddo

c     Find the partial pressures of organic species at the
c     particle surface.
c     ps(Pa) in common block in subroutine
cbk   removed - bkoo (03/09/03)
cbk      do k=1,n
cbk         call organics(k,q)
cbk      enddo

      return
      end



cbk      subroutine organics(k,q)

c     This is a dummy routine that returns partial pressures of organics
c     at the surface of the particles
c *************************************************************************
c               WRITTEN BY DR. CHRISTODOULOS PILINIS
c                          May 1997
c *************************************************************************
c
c
cbk      include 'dynamic.inc'
cbk      real*8 q(ntotal)
c
cbk      do i=1,norg
cbk         ps(k,IHCL+i) = 0.0d0
c         ps(k,IHCL+i) = 0.1d-5 ! in Pa 
cbk      enddo

cbk      return
cbk      end

