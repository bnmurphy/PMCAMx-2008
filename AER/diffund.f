c=========================================================================
c 03/09/03: bkoo
c           - added call to orgps2 to calculate PS for organics
c 03/07/03: bkoo
c           - removed redundant qtemp since step has now parameter for nsec
c           - removed redundant IF statement
c           - commented out icount-initializing code
c             (note that HYBR calls MADM multiple times)
c 11/20/02: tmg
c           - modified timestep check to include factor of qt (and changed
c             error tolerances
c 10/22/02: tmg
c           - added parameter to subroutine step because not all sections
c             are always in use
c 07/02/01: tmg
c           - changed to trajectory grid method for solving dynamic equations
c 01/27/01: bkoo
c           - change the STEP position
c           - water change in HYBR
c 01/26/01: bkoo
c           -
c 01/07/01: bkoo 
c           - try a new scheme to estimate xx 
c 10/08/00: bkoo
c           - move declaration statements before inline function declaration
c 09/19/00: bkoo
c           - combine diffundm & diffundh, dryinm & dryinh
c MAY 2000: bkoo
c           - combine MADM with HYBRID
c           - diffundm & dryinm: diffund & dryin for MADM
c           - diffundh & dryinh: diffund & dryin for HYBRID
cgy write error messages only to unit 6, only before a stop
c
c========================================================================+

      subroutine diffund(neq, t, q, tout)

c
c     This subroutine calculates the derivatives of the concetrations of all
c     the aerosol species, for all sections, as well as the concentrations of 
c     the gases for the new generation
c     multicomponent aerosol dynamics model (MADM)
c *************************************************************************
c               WRITTEN BY DR. CHRISTODOULOS PILINIS
c                          January-February 1998
c *************************************************************************
c

      include 'dynamic.inc'

      real*8 dq(ntotal)
      real*8 q(ntotal), q2(ntotal), qlast(ntotal)
      real*8 dqdt(ntotal)      ! Aerosol module master array for
c     concentration derivatives of aerosol  species in ug/m3/sec and
c     gases in ppm/sec.
      double precision KnD
      dimension dqsum(nsp)    ! dqdt mass balance for each species
      dimension pgas(ngas)    ! partial pressures of condensible gases (Pa)
c     add ht array back in (tmg, 07/02/01)
      dimension ht(nsec)     ! total rate of change
cbk      real*8 qtemp(ntotal) ! bkoo (03/07/03)
      real*8 merr(nsec), htavg(nsec) 
c     parameters for gas "source" term by the equilibrium sections
c     real*8  sdgas(ngas),sgasrate(ngas),residual(ngas)
      integer icase(nsec)
      
c     new variables (tmg, 07/02/01)
      dimension qtlast(nsec), hplast(nsec), hptlast(nsec)
      dimension htlast(nsec), hilast(nsec,ngas)
      real*8 masserr, errtol, ltol

      common / counters / icount ! tmg (01/21/02)
c     common / eqgas / sdgas,sgasrate,residual,t0c
      common /debug/ icnt
      common /drycase/ icase

c     Equation for transition regime (Dahneke 1983)
      f(yKnD, alfa)=(1.0d0+yKnD)/(1.0d0+2.0d0*yKnD*(1.0d0+yKnD)/alfa)

      masserr = 0.0d0
      tmincr=10.0d0
cbk      icount=0 ! bkoo (03/07/03)


      tpres=t
      errtol = 1.0d-3 
      ltol = 1.0d-6

c     tinys=1.0d-9             ! minimum mass in ug/m3
      tiny2= tinys * govra ! 02/14/02

      ng   = neq - ngas

 20   do i=1,neq
        dqdt(i)=0.0d0          ! Initialize derivatives
      enddo

c      if (icount.gt.200005) return ! tmg (12/17/02)

      do i=1,nsecx
         qt(i)=0.0d0            ! Initialize total masses
         ht(i)=0.0d0            ! Initialize total rates
      enddo

c     STEP subroutine updates c, ps, dry, q(H2O), keeping dqdt intact
cbk      if(aerm.eq.'MADM') call step(nsecx,q)  ! (10/22/02)
cbk      if(aerm.eq.'HYBR') then
cbk        do i=1,ntotal
cbk          if(i.le.ng) then
cbk            qtemp(i)=q(i)
cbk          else
cbk            qtemp(i)=0.0d0
cbk          endif
cbk        enddo
      call step(nsecx,q) ! bkoo (03/07/03)
cbk        call step(nsecx,qtemp) ! this will "front end" common variables
c                          ps(i,j) and c(i,j). So call STEP again before
c                          using these variables outside of the dynamics
c                          integrator (for equilibrium code, or printing, etc.)
cbk        do i=1,nsecx
cbk          q((i-1)*nsp+KH2O)=qtemp((i-1)*nsp+KH2O) ! bkoo (01/27/01)
cbk        enddo
cbk      endif

cbk   calculate the saturation concentrations for organics - bkoo (04/18/02) (03/09/03)
      call orgps2(q)

c If an aerosol oscillates between wet and dry states the numerical
c  integration can hardly converge because the PS calculations for
c  the two states are much different from each other. In that case,
c  we consider the aerosol as metastable liquid state to avoid the
c  oscillation.
c  - bkoo (03/05/02)
      if(isfirst) then
        isfirst=.FALSE.
        do i=1,nsecx
          icdry(i)=0
          dold(i)=dry(i)
        enddo
      else
        do i=1,nsecx
          if(ims(i).eq.0) then ! if not metastable
            if(dold(i) .and. .not.dry(i)) then
              if(icdry(i).eq.100) then
                ims(i)=1
c                call errprt('METASTABLE    ',t,0,i,1.d0,0)
              endif
              icdry(i)=icdry(i)+1
            endif
            dold(i)=dry(i)
          endif
        enddo
      endif

      do i=1,nsecx                                            ! dbg 
       do k=1,ngas                                            ! dbg 
        if(ps(i,k).lt.0. .and. ps(i,k).ge.-1.e-15) ps(i,k)=0. ! dbg 
        if(ps(i,k).lt.0.) then                                ! dbg 
c         call errprt('NEG-PS        ',t,k,i,ps(i,k),1)        ! dbg
c         stop                                                 ! dbg 
         ps(i,k)=0.                                           ! dbg
        endif                                                 ! dbg 
       enddo                                                  ! dbg 
       do k=1,ninti                                               ! dbg
        if(c(i,k).lt.0.d0) c(i,k)=0.d0 ! bkoo (02/14/02)            dbg
       enddo                                                      ! dbg
      enddo                                                   ! dbg 

c     do i=1,nsecp1
c      do k=1,nsp
c       intr(i,k)=0.0d0
c      end do
c     end do

      do i=1,nsecx              ! Calculate total mass in each section
         do k=1,nsp
            qt(i)=qt(i)+q((i-1)*nsp+k)	
         enddo
         if(qn(i+nsecx2).gt.0.and.qt(i).gt.tinys) then !avoid empty sections
            dsec(i+nsecx2)=(qt(i)/qn(i+nsecx2))**(0.33333)
cgy         if(dsec(i+nsecx2).lt.0.005) write(90,*) 'at t= ',t,
cgy    &          ' Diameter of section ',i+neqsec,
cgy    &          ' is smaller than 0.005 um at Dp = ',dsec(i+neqsec)
cgy         if(dsec(i+nsecx2).gt.20.0) write(90,*) 'at t= ',t,
cgy    &          ' Diameter of section ',i+neqsec,
cgy    &          ' is larger than 20.0 um at Dp = ',dsec(i+neqsec)
         endif
      enddo

c     calculate the partial pressures of the consensible gases
c     Note that since for gases, q is in ppm, to convert to Pa
c     we multiply by (1.0d-6 * 1.015d5)
c we multiply by (1.01325d-1 * pres) instead of 1.015d-1 (bkoo, 06/09/00)
c
      pgas(INH3)  =q(ng+INH3)   * 1.01325d-1 * pres
      pgas(IHNO3) =q(ng+IHNO3)  * 1.01325d-1 * pres
      pgas(IH2SO4)=q(ng+IH2SO4) * 1.01325d-1 * pres
      pgas(IHCL)  =q(ng+IHCL)   * 1.01325d-1 * pres
c
c     Do now the organic gases (organics follow HCl)
c
      do i=1,norg-1
         kk=IHCL+i
         pgas(kk)=q(ng+kk) * 1.01325d-1 * pres
      enddo

      call flush(6)
      do i=1,nsecx
       if(qt(i).le.tinys) then      ! no aerosol in this section-->
        do k=1,ngas
         hi(i,k)=0.0d0              ! no condensation/evaporation
        enddo
       else
        indx=(i-1)*nsp
        dp = 1.0d-6 * dsec(i+nsecx2) ! representative particle diameter ( m )
        if(.not.dry(i)) then         ! WET AEROSOL
         icase(i)=0                              
         do k=1,IHCL                ! inorganic
          call flush(6)
          KnD=2.0d0*lamda(k)/dp
          hi(i,k)=12.0d0*diffus(k)*(pgas(k)-ps(i,k))*
     &            f(KnD, delta(k))*(gmw(k)*1.0d-3)  /
     &            (1.0d3*dp**2*rgas*temp)
         enddo
	   
c ****  restrict NH3 transfer to force electroneutrality ***********

c ****  modified by bkoo (02/14/02)
c
         GG=2.d0*c(i,6)/intmw(6)+c(i,7)/intmw(7) ! [umol/m3]
     &          +c(i,8)/intmw(8)+c(i,5)/intmw(5)+c(i,19)/intmw(19)
     &          -c(i,3)/intmw(3)-c(i,4)/intmw(4)
c
         alpha=0.1*ABS(GG)           ! in umoles/m3/s 10% change/s in H+ conc.
c
c    First check for electroneutrality of the fluxes 
        elec=qt(i)*(2*hi(i,ih2so4)/gmw(ih2so4)+hi(i,ihno3)/gmw(ihno3)+
     &       hi(i,ihcl)/gmw(ihcl)-hi(i,inh3)/gmw(inh3))
c
c    If not within alpha find xx.  Iteration is occationally required 
c    due to precission issues and the sensitivity of the new elec to xx 
c
         if(elec.lt.-alpha) then ! flux over basic limit
          if(hi(i,ihno3).lt.0.d0.and.hi(i,ihcl).lt.0.d0) then ! both NO3 & Cl evaporate
           x1=-alpha/qt(i)-2.d0*hi(i,ih2so4)/gmw(ih2so4)
     &        +hi(i,inh3)/gmw(inh3)
           if(x1.eq.0.d0) then
            hi(i,ihno3)=0.d0
            hi(i,ihcl)=0.d0
           else
            xx=(hi(i,ihno3)/gmw(ihno3)+hi(i,ihcl)/gmw(ihcl))/x1
            if(xx.ge.1.d0) then
             hi(i,ihno3)=hi(i,ihno3)/xx
             hi(i,ihcl)=hi(i,ihcl)/xx
            else
             hi(i,ihno3)=0.d0
             hi(i,ihcl)=0.d0
            endif
           endif
          elseif(hi(i,ihno3).lt.0.d0) then ! Cl condenses
           x1=-alpha/qt(i)-2.d0*hi(i,ih2so4)/gmw(ih2so4)
     &        -hi(i,ihcl)/gmw(ihcl)+hi(i,inh3)/gmw(inh3)
           if(x1.eq.0.d0) then
            hi(i,ihno3)=0.d0
           else
            xx=hi(i,ihno3)/gmw(ihno3)/x1
            if(xx.ge.1.d0) then
             hi(i,ihno3)=hi(i,ihno3)/xx
            else
             hi(i,ihno3)=0.d0
            endif
           endif
          elseif(hi(i,ihcl).lt.0.d0) then ! NO3 condenses
           x1=-alpha/qt(i)-2.d0*hi(i,ih2so4)/gmw(ih2so4)
     &        -hi(i,ihno3)/gmw(ihno3)+hi(i,inh3)/gmw(inh3)
           if(x1.eq.0.d0) then
            hi(i,ihcl)=0.d0
           else
            xx=hi(i,ihcl)/gmw(ihcl)/x1
            if(xx.ge.1.d0) then
             hi(i,ihcl)=hi(i,ihcl)/xx
            else
             hi(i,ihcl)=0.d0
            endif
           endif
          endif ! both NO3 & Cl evaporate
         elseif(elec.gt.alpha .and. hi(i,inh3).lt.0.d0) then ! flux over acidic limit
          x1=2.d0*hi(i,ih2so4)/gmw(ih2so4)+                   ! & NH4 evaporates
     &       hi(i,ihno3)/gmw(ihno3)+hi(i,ihcl)/gmw(ihcl)-alpha/qt(i)
          if(x1.eq.0.d0) then
           hi(i,inh3)=0.d0
          else
           xx=hi(i,inh3)/gmw(inh3)/x1
           if(xx.ge.1.d0) then
            hi(i,inh3)=hi(i,inh3)/xx
           else
            hi(i,inh3)=0.d0
           endif
          endif
         endif ! flux over basic limit
c
         elec2=qt(i)*(2.d0*hi(i,ih2so4)/gmw(ih2so4)+hi(i,ihno3)/
     &         gmw(ihno3)+hi(i,ihcl)/gmw(ihcl)-hi(i,inh3)/gmw(inh3))
         if(elec2.lt.-alpha) then
          hmax=12.d0*diffus(inh3)*pgas(inh3)
     &        *f(2.d0*lamda(inh3)/dp,delta(inh3))
     &        *(gmw(inh3)*1.d-3)/(1.d3*dp**2*rgas*temp)
          hi(i,inh3)=MIN((2.d0*hi(i,ih2so4)/gmw(ih2so4)
     &              +hi(i,ihno3)/gmw(ihno3)+hi(i,ihcl)/gmw(ihcl)
     &              +alpha/qt(i))*gmw(inh3),hmax)
         elseif(elec2.gt.alpha) then
          hmax1=12.d0*diffus(ihno3)*pgas(ihno3)
     &         *f(2.d0*lamda(ihno3)/dp,delta(ihno3))
     &         *(gmw(ihno3)*1.d-3)/(1.d3*dp**2*rgas*temp)
          hmax2=12.d0*diffus(ihcl)*pgas(ihcl)
     &         *f(2.d0*lamda(ihcl)/dp,delta(ihcl))
     &         *(gmw(ihcl)*1.d-3)/(1.d3*dp**2*rgas*temp)
          if(hi(i,ihno3).gt.0.d0 .and. hi(i,ihcl).gt.0.d0) then
           wf=hi(i,ihcl)/gmw(ihcl) /
     &       (hi(i,ihcl)/gmw(ihcl)+hi(i,ihno3)/gmw(ihno3))
           hi(i,ihno3)=MIN((alpha/qt(i)+hi(i,inh3)/gmw(inh3)
     &                -2.d0*hi(i,ih2so4)/gmw(ih2so4))*gmw(ihno3)
     &                *(1-wf),hmax1)
           hi(i,ihcl)=MIN((alpha/qt(i)+hi(i,inh3)/gmw(inh3)
     &               -2.d0*hi(i,ih2so4)/gmw(ih2so4))*gmw(ihcl)
     &               *wf,hmax2)
          elseif(hi(i,ihno3).gt.0.d0) then
           hi(i,ihno3)=MIN((alpha/qt(i)+hi(i,inh3)/gmw(inh3)
     &                -2.d0*hi(i,ih2so4)/gmw(ih2so4)
     &                -hi(i,ihcl)/gmw(ihcl))*gmw(ihno3),hmax1)
          elseif(hi(i,ihcl).gt.0.d0) then
           hi(i,ihcl)=MIN((alpha/qt(i)+hi(i,inh3)/gmw(inh3)
     &               -2.d0*hi(i,ih2so4)/gmw(ih2so4)
     &               -hi(i,ihno3)/gmw(ihno3))*gmw(ihcl),hmax2)
          endif
         endif                                              ! dbg

c *******************************************************************

       else                        ! DRY AEROSOL
        call dryin(i,dp,ht,pgas)   ! dry inorganics
       endif                       ! if WET or DRY

c
c     if very small quantities are evaporating - inorganics (02/14/02)
c
        if(q(indx+KNO3).lt.tinys .and. hi(i,ihno3).lt.0.d0) then ! KNO3
          tmp=hi(i,inh3)-hi(i,ihno3)/gmw(ihno3)*gmw(inh3)
          if((q(indx+KNH4).ge.tinys .and. tmp.le.0.d0) .or.
     &       (tmp.ge.0.d0 .and. q(ng+inh3).ge.tiny2)) then
           hi(i,inh3)=tmp
          else
           hi(i,inh3)=0.d0
           if(q(indx+KCL).ge.tinys) then
            hi(i,ihcl)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihcl)
           else
            hi(i,ihcl)=0.d0
           endif
          endif
          hi(i,ihno3)=0.d0
        endif ! KNO3
        if(q(indx+KCL) .lt.tinys .and. hi(i,ihcl) .lt.0.d0) then ! KCL
          tmp=hi(i,inh3)-hi(i,ihcl)/gmw(ihcl)*gmw(inh3)
          if((q(indx+KNH4).ge.tinys .and. tmp.le.0.d0) .or.
     &       (tmp.ge.0.d0 .and. q(ng+inh3).ge.tiny2)) then
           hi(i,inh3)=tmp
          else
           hi(i,inh3)=0.d0
           if(q(indx+KNO3).ge.tinys) then
            hi(i,ihno3)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihno3)
           else
            hi(i,ihno3)=0.d0
           endif
          endif
          hi(i,ihcl)=0.d0
        endif ! KCL
        ! now, "hi < 0" means "q(aer) >= tinys" for HNO3 & HCl
        if(q(indx+KNH4).lt.tinys .and. hi(i,inh3) .lt.0.d0) then ! KNH4
          if(hi(i,ihno3).lt.0.d0 .and. hi(i,ihcl).lt.0.d0) then
           wf=hi(i,ihcl)/gmw(ihcl) /
     &        (hi(i,ihcl)/gmw(ihcl)+hi(i,ihno3)/gmw(ihno3))
           hi(i,ihcl)=hi(i,ihcl)-hi(i,inh3)/gmw(inh3)*gmw(ihcl)*wf
           hi(i,ihno3)=hi(i,ihno3)-hi(i,inh3)/gmw(inh3)*gmw(ihno3)
     &                 *(1-wf)
           if(hi(i,ihcl).gt.0.d0 .and. q(ng+ihcl).lt.tiny2) then
            hi(i,ihcl)=0.d0
            hi(i,ihno3)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihno3)
           endif
           if(hi(i,ihno3).gt.0.d0 .and. q(ng+ihno3).lt.tiny2) then
            hi(i,ihno3)=0.d0
            hi(i,ihcl)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihcl)
           endif
          elseif(hi(i,ihno3).lt.0.d0) then
           hi(i,ihno3)=hi(i,ihno3)-hi(i,inh3)/gmw(inh3)*gmw(ihno3)
           if(hi(i,ihno3).gt.0.d0 .and. q(ng+ihno3).lt.tiny2) then
            hi(i,ihno3)=0.d0
            hi(i,ihcl)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihcl)
           endif
          elseif(hi(i,ihcl).lt.0.d0) then
           hi(i,ihcl)=hi(i,ihcl)-hi(i,inh3)/gmw(inh3)*gmw(ihcl)
           if(hi(i,ihcl).gt.0.d0 .and. q(ng+ihcl).lt.tiny2) then
            hi(i,ihcl)=0.d0
            hi(i,ihno3)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihno3)
           endif
          endif
          hi(i,inh3)=0.d0
        endif ! KNH4
c   Special case for if a dry particle with NaNO3 and case 1.1.3 runs out of 
c   NH4NO3 to evaporate. In this case NO3 from NaNO3 will be lost from the aerosol
c   and NaCl will form which will drive towards a different equilibrium.
c   To prevent the inappropriate formation of NaCl when the NH4NO3 runs out
c   we will set the NO3 flux to zero
        if(icase(i).eq.113) then
          if(hi(i,ihno3).lt.0.d0 .and. c(i,13).le.tinys) then
           hi(i,inh3)=hi(i,inh3)-hi(i,ihno3)/gmw(ihno3)*gmw(inh3)
           hi(i,ihno3)=0.d0
           if(hi(i,inh3).gt.0.d0 .and. q(ng+inh3).lt.tiny2) then
            hi(i,inh3)=0.d0
            if(q(indx+KCL).ge.tinys) then
             hi(i,ihcl)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihcl)
            else
             hi(i,ihcl)=0.d0
            endif
           endif
          endif
        endif
c   Special case: A dry particle contains NH4Cl and uses case 1.1.2 but runs out of 
c   NaCl to evaporate. In this case NH4NO3 will form from NaNO3 regardless of the
c   potential of NH4NO3 to evaporate  To prevent the 
c   inappropriate formation of NH4NO3 when the NaCl runs out we will set the Cl flux
c   to equal the ammonia flux
        if(icase(i).eq.112) then
          if(hi(i,ihcl).lt.0.d0 .and. c(i,9).le.tinys) then
           hi(i,ihcl)=hi(i,inh3)/gmw(inh3)*gmw(ihcl)
           if(q(indx+KNO3).ge.tinys) then
            hi(i,ihno3)=-2.d0*hi(i,ih2so4)/gmw(ih2so4)*gmw(ihno3)
           else
            hi(i,ihno3)=0.d0
           endif
           if(hi(i,ihcl).gt.0.d0 .and. q(ng+ihcl).lt.tiny2) then
            hi(i,ihcl)=0.d0
            hi(i,inh3)=0.d0
           endif
          endif
        endif
c
        do k=1,norg-1                ! do organics
         kk=IHCL+k
         KnD=2.0d0*lamda(kk)/dp
         hi(i,kk)=12.0d0*diffus(kk)*(pgas(kk)-ps(i,kk))*
     &            f(KnD, delta(kk))*(gmw(kk)*1.0d-3) /
     &            (1.0d3*dp**2*rgas*temp)
c
c     if very small quantities are evaporating - organics (02/14/02)
c
         ii=nexti+k
         if(q(indx+ii).lt.tinys .and. hi(i,kk).lt.0.d0) hi(i,kk)=0.d0
        enddo
       endif                        ! no aerosol in this section
      enddo                         ! do nsecx


c       **************************************************************
c         MODIFIED this section to calculate as in Chock and Winkler
c		tmg(07/02/01)

c         Take the derivative of the total flux for each particle
c     Calculate the derivatives using Equation 11 of Pilinis et al. (1987)
c
c

      do i=1,nsecx
        ht(i)=0.0d0
        do k=1,ngas
          ht(i)=ht(i)+hi(i,k)
        enddo
      enddo

c     makes sure previous step wasn't too large (tmg, 07/02/01)
      if (tpres.gt.t) then
        do i=1,nsecx
           htavg(i)=(ht(i)+htlast(i))/2.
           merr(i) = qt(i)*(htavg(i)-htlast(i))*tmincr ! tmg (11/20/02)
           masserr = max(masserr,abs(merr(i)))
        enddo

c      if error small, increase tmincr
        if (masserr.lt.ltol/100) then
          tmincr= tmincr*2.0d0
        else if (masserr.lt.errtol/10.) then
          tmincr = tmincr*1.1d0
        endif
	  
c      if error large cut timestep in half and redo last iteration
        if (masserr.gt.errtol.and.tmincr.gt.1.0d-4) then
 30       tmincr = tmincr/2.0d0
          do k=1,ntotal
            q(k) = qlast(k)
          enddo
          do i=1,nsecx
            ht(i)=htlast(i)
            qt(i)=qtlast(i)
            do k=1,ngas
              hi(i,k)=hilast(i,k)
            enddo
          enddo
          tpres=tlast
        else if (tpres.eq.(tout)) then
 40       t=tpres
          return 
        endif
      endif 

c      make sure tout will not be exceeded
        if ((tpres+tmincr).gt.(tout)) then
          tmincr= tout-tpres
        endif

      masserr = 0.0d0

c       **************************************************************
c         MODIFIED this section to calculate as in Chock and Winkler
c		tmg(07/02/01)

c         Take the derivative of the total flux for each particle
c     Calculate the derivatives using Equation 11 of Pilinis et al. (1987)
c
      do i=1,nsecx
         indx=(i-1)*nsp
c         dqdt(indx+KH2O)= 0.0d0	       
c         dqdt(indx+KNa) = 0.0d0	       
         dqdt(indx+KSO4)= hi(i,IH2SO4)*qt(i)
         dqdt(indx+KNO3)= hi(i,IHNO3) *qt(i)
         dqdt(indx+KNH4)= hi(i,INH3)  *qt(i)
         dqdt(indx+KCL) = hi(i,IHCL)  *qt(i)
         do k=1,norg-1
            ii=nexti+k
            kk=IHCL+k
            dqdt(indx+ii) = hi(i,kk) *qt(i)
        enddo
      enddo

c
c     Done with aerosol rates. Start gases.
c
c     changed according to T-G model (tmg, 07/02/01)

c     initialize
      do i = 1,ngas
         dq(ng+i) = 0.0d0
      enddo

      do i=1,nsecx
        dqdt(ng+IH2SO4)=-hi(i,IH2SO4)*qt(i)
        dqdt(ng+IHNO3) =-hi(i,IHNO3) *qt(i)
        dqdt(ng+INH3)  =-hi(i,INH3)  *qt(i)
        dqdt(ng+IHCL)  =-hi(i,IHCL)  *qt(i)
        if (ht(i) .ne. 0.0d0) then
          dq(ng+IH2SO4)= dq(ng+IH2SO4) + 
     &      dqdt(ng+IH2SO4)/ht(i)*(exp(ht(i)*tmincr)-1.0d0)
          dq(ng+IHNO3) = dq(ng+IHNO3) 
     &      +dqdt(ng+IHNO3)/ht(i)*(exp(ht(i)*tmincr)-1.0d0)
          dq(ng+INH3) = dq(ng+INH3)
     &      +dqdt(ng+INH3)/ht(i)*(exp(ht(i)*tmincr)-1.0d0)
          dq(ng+IHCL) = dq(ng+IHCL)
     &      +dqdt(ng+IHCL)/ht(i)*(exp(ht(i)*tmincr)-1.0d0)
        else
          dq(ng+IH2SO4) = dq(ng+IH2SO4)+dqdt(ng+IH2SO4)*tmincr
          dq(ng+IHNO3) = dq(ng+IHNO3)+dqdt(ng+IHNO3)*tmincr
          dq(ng+INH3) = dq(ng+INH3)+dqdt(ng+INH3)*tmincr
          dq(ng+IHCL) = dq(ng+IHCL)+dqdt(ng+IHCL)*tmincr
        endif
        do k=1,norg-1
          ii=ng+IHCL+k
          km=IHCL+k
          dqdt(ii) =-hi(i,km)*qt(i)
          if (ht(i) .ne. 0.0d0) then
            dq(ii) = dq(ii)
     &        +dqdt(ii)/ht(i)*(exp(ht(i)*tmincr)-1.0d0)
          else
            dq(ii) = dq(ii) + dqdt(ii)*tmincr
          endif
        enddo
      enddo

c
c     Check derivative balance to ensure no change in mass will occur
c	Ignore for T-G method (tmg, 07/02/01)
c      To change the order of the summation later... (bkoo)
      dqdtsum=0.0d0

c      do i=1,neq
c	 dqdtsum=dqdtsum+dqdt(i)
c      enddo
c
c     -- an error of 1e-5 ug/m3/s will produce 1ug/m3 in 24 hours --
      if(abs(dqdtsum).ge.1e-5) then   
cgy          write(6,*)'mass error of ',dqdtsum,' ug/m3/s occured at ',t
cgy          write(90,*)'mass error of ',dqdtsum,' ug/m3/s occured at ',t
         do ispec=1,nsp
          dqsum(ispec)=0.0d0
          do isec=1,nsecx
             dqsum(ispec)=dqsum(ispec)+dqdt((isec-1)*nsp+ispec)
          enddo
          if(ispec.eq.kso4) dqsum(ispec)=dqsum(ispec)+dqdt(ng+ih2so4)
          if(ispec.eq.kno3) dqsum(ispec)=dqsum(ispec)+dqdt(ng+ihno3)
          if(ispec.eq.knh4) dqsum(ispec)=dqsum(ispec)+dqdt(ng+inh3)
          if(ispec.eq.kcl)  dqsum(ispec)=dqsum(ispec)+dqdt(ng+ihcl)
          if(ispec.gt.kcl) then
             dqsum(ispec)=dqsum(ispec)+dqdt(ng+ihcl+ispec-kcl)
          endif
ctmg          do isec=1,nsecx
cgy             write(90,21)t,ispec,isec+neqsec,dqsum(ispec)/emw(ispec),
cgy     &                dqdt((isec-1)*nsp+ispec)/emw(ispec)
ctmg           enddo
          enddo
      endif
 21   format('dqdtsum error: ',e14.5,2i5,2e14.5)
c
c     Add gas phase source term in ug/m3/sec
c	Modified according to T-G method (tmg, 07/02/01)

      do igas=1,ngas
       dqdt(ng+igas)=gsource(igas)
      enddo
c
c     Change the units of the gases from ug/m3/sec to ppm/sec
c     add pres (bkoo, 06/09/00)
      do i=1,ngas
       dqdt(ng+i)=dqdt(ng+i)*rgas*temp/(1.01325d5*pres*gmw(i))
       dq(ng+i)=dq(ng+i)*rgas*temp/(1.01325d5*pres*gmw(i))
      enddo

c	additional code for T-G method (tmg, 07/02/01)

      do k=1,ntotal
        qlast(k) = q(k)
      enddo

c       update concentrations according to new method
      do i = 1,ngas
        q(ng+i) = q(ng+i) + dq(ng+i)+ dqdt(ng+i)*tmincr
      enddo

      do i=1,nsecx
        htlast(i) = ht(i)
        qtlast(i) = qt(i)
        do k=1,ngas
          hilast(i,k)=hi(i,k)
        enddo
      enddo

      do i=1,nsecx
        do k=2,nsp
          if (ht(i).ne.0.0d0) then
            q((i-1)*nsp+k)=q((i-1)*nsp+k)+dqdt((i-1)*nsp+k)/ht(i)*
     &                     (exp(ht(i)*tmincr)-1.0d0)
          else
            q((i-1)*nsp+k)=q((i-1)*nsp+k)+dqdt((i-1)*nsp+k)*tmincr
          endif
        enddo
      enddo

cbk      if(aerm.eq.'HYBR') then
cbk       call negchk(tpres,q,nsecd)
cbk      else
cbk       call negchk(tpres,q,nsec)
cbk      endif
      call negchk(tpres,q,nsecx) ! bkoo (03/07/03)

      icount = icount + 1  ! tmg (12/17/02)

      tlast = tpres
      tpres = tpres+tmincr

      go to 20
      end


      subroutine dryin(i,dp,ht,pgas)

c     This subroutine calculates the intrasectional
c     coefficients for the 4 inorganic condensible gases
c     in the subcase of dry aerosols, for the new
c     multicomponent aerosol dynamics model (MADM)
c *************************************************************************
c               WRITTEN BY DR. CHRISTODOULOS PILINIS
c                          January 1998
c *************************************************************************
c           MODIFIED TO INCLUDE SUPERSATURATION FLAGS BY
C                       KEVIN P. CAPALDO
C                           JUNE 1998
c *************************************************************************
c  For the Ammonia rich case 9 subcases are used which utilize at most 2
c  of the following three reversible reactions
c   (1)  NH4NO3(s) <=> NH3(g) + HNO3(g)
c   (2)  NH4Cl(s)  <=> NH3(g) + HCl(g)
c   (3)  NaCl(s) + HNO3(g) <=>  NaNO3(s) + HCl(g)
c  Thermodynamically all three equations cannot be in equilibrium simultaneously
c  The choice of which reactions to use is determined by the compounds present
c  in the aerosol phase and the gas phase.  Because all 4 of the salts NaCl, NaNO3,
c  NH4Cl, and NH4NO3 cannot exist at the same time in the aerosol phase,
c  the assumption that NH4NO3 and NaCl will never coexist (these compounds will
c  react to form NaNO3 and NH4Cl) is used. Because of this assumption,
c  supersaturations in the gas phase do not always invoke the use of
c  the corresponding reaction.
c
      include 'dynamic.inc'
      include 'equaer.inc'
      dimension pgas(ngas)              ! partial pressures of condensible gase (Pa)
      dimension ht(nsecx)               ! total rate of change 
      integer icase(nsec)
      common /drycase/ icase
c     Equation for transition regime (Dahneke 1983)
      f(yKnD, alfa)=(1.0d0+yKnD)/(1.0d0+2.0d0*yKnD*(1.0d0+yKnD)/alfa)

c      tinys =1.0d-9          ! minimum mass in ug/m3

      if(qt(i).le.tinys) then       ! no aerosol in this section---->
         hi(i,IH2SO4)=0.0d0         ! no condensation/evaporation
         hi(i,INH3)  =0.0d0
         hi(i,IHNO3) =0.0d0
         hi(i,IHCL)  =0.0d0
         goto 200
      endif
c     change the units for R to use it for equilibrium constants
      RGAS1=RGAS/(1.01325d5*pres)   ! add pres (bkoo, 06/09/00)

      hno3kn=2.0d0*lamda(ihno3)/dp
      fhno3=f(hno3kn, delta(ihno3))
      dhno3=diffus(ihno3)*fhno3

      ammokn=2.0d0*lamda(inh3)/dp
      fammo=f(ammokn, delta(inh3))
      dammo=diffus(inh3)*fammo

      hclkn=2.0d0*lamda(ihcl)/dp
      fhcl=f(hclkn, delta(ihcl))
      dhcl=diffus(ihcl)*fhcl

      sulkn=2.0d0*lamda(IH2SO4)/dp
      fsul=f(sulkn, delta(IH2SO4))
      dsul=diffus(IH2SO4)*fsul
c
c     Calculate the rate constants and flags for condensation
c     Make the reaction consistent with the reactions used
c     by ISORROPIA. The equilibrium constants used by MADM
c     are caclulated as a combination of the equilibrium
c     constants used by ISORROPIA.
c
      RKNH4NO3=xk10/(rgas1*temp)**2
      RKNA=XK4*XK8/(XK3*XK9)
      RKNH4CL=xk6/(rgas1*temp)**2
      satnh4no3=pgas(ihno3)*pgas(inh3)/(rgas*temp)**2
      satnh4cl=pgas(ihcl)*pgas(inh3)/(rgas*temp)**2   
      if(satnh4no3.gt.rknh4no3) then
       ifgnh4no3 = 1                ! NH4NO3 forms
      else
       ifgnh4no3 = 0
      endif
      if(satnh4cl.gt.rknh4cl) then
       ifgnh4cl = 1                 ! NH4CL forms
      else
       ifgnh4cl = 0
      endif
      if(pgas(ihcl).gt.rkna*pgas(ihno3)) then
       ifgna = 1                    ! NaCl forms if NaNO3 exists
      else
       ifgna = 0                    ! NaNO3 forms if NaCl exists
      endif
c
      if(ifgnh4no3.eq.1.and.ifgnh4cl.eq.0.and.ifgna.eq.1) then
cgy       write(90,*)'error in gas phase conc. (or solid eq rates)'
       write(6,*)'error in gas phase conc. (or solid eq rates)'
       stop
      elseif(c(i,9).gt.tinys.and.c(i,13).gt.tinys) then
cgy       write(90,*)'error: NaCl and NH4NO3 should not coexist'
       write(6,*)'error: NaCl and NH4NO3 should not coexist'
       stop
      elseif(c(i,12).gt.tinys.and.
     &      (c(i,9).gt.tinys.or.c(i,11).gt.tinys)) then
cgy       write(90,*)'(NH4)2SO4 and NaCl or NaNO3 should not coexist'
       write(6,*)'(NH4)2SO4 and NaCl or NaNO3 should not coexist'
       stop
      endif
c
      k=IH2SO4                      ! do H2SO4
      hi(i,k)=12.0d0*dsul*pgas(k)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)
c
c     Now find the subcases to calculate the fluxes
c     for HNO3, NH3 and HCl
c
c     Ammonia poor case
      if(c(i,16).gt.tinys.or.c(i,17).gt.tinys.or.c(i,18).gt.tinys)
     &   goto 120                   ! Acidic aerosol
c
c     Ammonia Rich cases
c
c     No Sodium (except Na2SO4)
      if(c(i,9).le.tinys.and.c(i,11).le.tinys) then    ! ~NaCl and ~NaNO3
       if(c(i,13).gt.tinys.or.ifgnh4no3.eq.1) then     ! NH4NO3
        if(c(i,14).gt.tinys.or.ifgnh4cl.eq.1) goto 113 ! NH4Cl
        goto 116                                       ! ~NH4Cl
       elseif(c(i,14).gt.tinys.or.ifgnh4cl.eq.1) then
        goto 117                                       ! ~NH4NO3 and NH4Cl
       else
        goto 114                                       ! ~NH4NO3 and ~NH4Cl
       endif
c     NaCl present
      elseif(c(i,9).gt.tinys) then                     ! NaCl
       if(c(i,11).le.tinys.and.ifgna.eq.1) goto 117    ! ~NaNO3
       goto 112                                        ! NaNO3
c     NaNO3 present 
      else                                             ! NaNO3
       if(c(i,13).gt.tinys) then    ! NH4NO3 and NaCl wont form
        if(c(i,14).le.tinys.and.ifgnh4cl.eq.0.and.ifgna.eq.0)
     &     goto 116                                    ! ~NH4Cl
        goto 113                                       ! NH4Cl
       elseif(ifgna.eq.1) then                         ! NaCl may form
        if(ifgnh4no3.eq.1) goto 111                    ! NH4NO3 may form
        goto 112                                       ! NH4NO3 may form
       elseif(c(i,14).gt.tinys.or.ifgnh4cl.eq.1) then  ! NH4Cl
        goto 113                               
       else                                            ! ~NH4Cl
        goto 116
       endif
      endif
c
c     NH3 poor case, i.e. acidic aerosol (case 1.2)
c
 120  icase(i)=120
cgy   if(c(i,9).gt.tinys.or.c(i,11).gt.tinys.or.c(i,13).gt.tinys.or.
cgy  #   c(i,14).gt.tinys)write(90,*)'ERROR: SALTS AND ACID PRESENT ',
cgy  #   c(i,9),c(i,11),c(i,13),c(i,14)
c
      k=INH3                                           ! do NH3
      hi(i,k)=12.0d0*dammo*pgas(k)*(gmw(k)*1.0d-3)/
     #        (1.0d3*dp**2*rgas*temp)
      hi(i,IHCL)=0.0d0
      hi(i,IHNO3)=0.0d0

      goto 200
c
c     case 1.1.1 ***************************************************
c     Reactions 1 and 3
c     No NH4Cl but NaNO3, NaCl, and NH4NO3
 111  icase(i)=111
c     this case is not applied because NH4NO3 and NaCl cannot coexist
c
c     However, for the specific cases where NaCl and NH4NO3 are absent from
c     the aerosol phase but are suggested to form by supersaterations
c     of the gas phase components of reactions 1 and 3 (ifgnh4no3=1 and 
c     ifgna=1), we must determine which compound (NaCl or NH4NO3) will form.
c     This is done by determining whether NO3 will condense or evaporate
c     if both reactions 1 and 3 are employed.  If case 1.1.1 indicates net
c     condensation of NO3 then any NaCl that forms will be imeadiately 
c     converted to NaNO3 and NH4Cl, if net evaporation of NO3 is indicated
c     then reation 3 dominates and the NaCl formed is sufficient to react
c     with any NH4NO3 that condenses and keep the NH4NO3 concentration zero.
c
      que=(dammo*pgas(inh3)-dhno3*pgas(ihno3)-dhcl*pgas(ihcl)
     #    -2.0d0*dsul*pgas(ih2so4))/(rgas*temp)
c
      CHNO3=(-que+sqrt(que**2+4.0*dammo*RKNH4NO3*(dhno3+RKNA*dhcl
     #      )))/(2.0d0*(dhno3+RKNA*dhcl))
      CNH3=RKNH4NO3/CHNO3
      CHCL=RKNA*CHNO3

      k=IHNO3                                          ! do HNO3
      hi(i,k)=12.0d0*dhno3*(pgas(k)-CHNO3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=INH3                                           ! do NH3
      hi(i,k)=12.0d0*dammo*(pgas(k)-CNH3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=IHCL                                           ! do HCL
      hi(i,k)=12.0d0*dhcl*(pgas(k)-CHCL*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

c     if NO3 evaporates then NaCl forms and NH4NO3 does not form
      if(hi(i,ihno3).lt.0.0d0) goto 112
c     else NH4NO3 forms and NaCl does not form
      goto 113

c     case 1.1.2  ***************************************************
c     Reactions 2 and 3
c     No NH4NO3 but NaCl, NaNO3, and NH4Cl can all occur
 112  icase(i)=112
c
      que=(dammo*pgas(inh3)-dhno3*pgas(ihno3)-dhcl*pgas(ihcl)
     #    -2.0d0*dsul*pgas(ih2so4))/(rgas*temp)

      CHNO3=(-que+sqrt(que**2+4.0*dammo*RKNH4CL*(dhno3+RKNA*dhcl
     #      )/RKNA))/(2.0d0*(dhno3+RKNA*dhcl))
      CNH3=RKNH4CL/RKNA/CHNO3
      CHCL=RKNA*CHNO3

      k=IHNO3                                          ! do HNO3
      hi(i,k)=12.0d0*dhno3*(pgas(k)-CHNO3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=INH3                                           ! do NH3
      hi(i,k)=12.0d0*dammo*(pgas(k)-CNH3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=IHCL                                           ! do HCL
      hi(i,k)=12.0d0*dhcl*(pgas(k)-CHCL*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)
c
c     if NaCl doesn't exist and case 1.1.2 has been employed then
c     both NaNO3 must exist and reaction 3 indicates formation of
c     NaCl (ifgna=1).  If H2SO4 condenses then we must determine whether
c     there is sufficient NaCl formed to react with the condensing
C     H2SO4 (to form Na2SO4 and NH4Cl) or if NaNO3 reacts with H2SO4 
c     (and NH4NO3 is formed).  If the later is the case then 1.1.3 
c     is required.  To determine this the HCl flux is compared to the 
c     NH3 flux.  If JNH3 > JHCl (in molar units) then NaCl does not form.
      if(c(i,9).le.tinys.and.
     &   (hi(i,inh3)/gmw(inh3)).gt.(hi(i,ihcl)/gmw(ihcl))) goto 113
c
c     if there is no NH4 to evaporate then
      if(c(i,14).le.tinys.and.hi(i,inh3).lt.0.0d0) goto 115
c     if there is no NO3 to evaporate then 
      if(c(i,11).le.tinys.and.hi(i,ihno3).lt.0.0d0) goto 117
c     if there is no Cl to evaporate then 
      if(c(i,9).le.tinys.and.c(i,14).le.tinys.and.
     &   hi(i,ihcl).lt.0.0d0) goto 116
      goto 200
c
c     case 1.1.3  *********************************************
c     Reactios 1 and 2
c     No NaCl but NH4Cl, NH4NO3, and NaNO3 can occur
 113  icase(i)=113
      que=(dammo*pgas(inh3)-dhno3*pgas(ihno3)-dhcl*pgas(ihcl)
     #    -2.0d0*dsul*pgas(ih2so4))/(rgas*temp)

      CHNO3=(-que+sqrt(que**2+4.0*dammo*RKNH4NO3*(dhno3+dhcl
     #      *RKNH4CL/RKNH4NO3)))/(2.0d0*(dhno3+dhcl*RKNH4CL/RKNH4NO3))
      CNH3=RKNH4NO3/CHNO3
      CHCL=CHNO3*RKNH4CL/RKNH4NO3

      k=IHNO3                                          ! do HNO3
      hi(i,k)=12.0d0*dhno3*(pgas(k)-CHNO3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=INH3                                           ! do NH3
      hi(i,k)=12.0d0*dammo*(pgas(k)-CNH3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=IHCL                                           ! do HCL
      hi(i,k)=12.0d0*dhcl*(pgas(k)-CHCL*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)
c
c     if there is no NH4 to evaporate then
      if(c(i,13).le.tinys.and.c(i,14).le.tinys.and.
     &   hi(i,inh3).lt.0.0d0) then
       if(c(i,11).le.tinys) goto 114
       goto 119
      endif
c     if there is no Cl to evaporate then
      if(c(i,14).le.tinys.and.hi(i,ihcl).lt.0.0d0) goto 116
c     if there is no NO3 to evaporate then
      if(c(i,13).le.tinys.and.c(i,11).le.tinys.and.
     &   hi(i,ihno3).lt.0.0d0) goto 117

      goto 200
c
c     case 1.1.4  *********************************************
c     None of reactions 1, 2, or 3
c     Sulfates only (no nitrates or chlorides)
 114  icase(i)=114
      k=INH3                                           ! do NH3
      hi(i,k)=2.0d0*hi(i,IH2SO4)*(gmw(k)/gmw(IH2SO4))
      hi(i,k)=min(hi(i,k),12.0d0*dammo*pgas(k)*(gmw(k)*1.0d-3)/
     #        (1.0d3*dp**2*rgas*temp))
      hi(i,IHCL)=0.0d0
      hi(i,IHNO3)=0.0d0
      goto 200
c
c     case 1.1.5  *********************************************
c     Reaction 3 only with the flux of NH3 set to zero
c     Sodium system (No ammonium salts)
 115  icase(i)=115

      ss=(dhno3*pgas(ihno3)+dhcl*pgas(ihcl)
     #   +2.0d0*dsul*pgas(ih2so4))/(rgas*temp)

      CHNO3=ss/(dhno3+dhcl*RKNA)
      CHCL=CHNO3*RKNA

      k=IHNO3                                          ! do HNO3
      hi(i,k)=12.0d0*dhno3*(pgas(k)-CHNO3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=INH3                                           ! do NH3
      hi(i,k)=0.0d0

      k=IHCL                                           ! do HCL
      hi(i,k)=12.0d0*dhcl*(pgas(k)-CHCL*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

c     We assume here that case 1.1.5 is only called through case 1.1.2
c     if there is no Cl to evaporate then
      if(c(i,9).le.tinys.and.hi(i,ihcl).lt.0.0d0) goto 119
c     if there is no NO3 to evaporate then
      if(c(i,11).le.tinys.and.hi(i,ihno3).lt.0.0d0) goto 118

      goto 200
c
c     case 1.1.6  *********************************************
c     Only reaction 1
c     Nitrate system (no Chloride)
 116  icase(i)=116
      p=(dammo*pgas(inh3)-dhno3*pgas(ihno3)
     #  -2.0d0*dsul*pgas(ih2so4))/(rgas*temp)

      CHNO3=(-p+sqrt(p**2+4.0*dammo*RKNH4NO3*dhno3))
     #      /(2.0d0*dhno3)
      CNH3=RKNH4NO3/CHNO3

      k=IHNO3                                          ! do HNO3
      hi(i,k)=12.0d0*dhno3*(pgas(k)-CHNO3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=INH3                                           ! do NH3
      hi(i,k)=12.0d0*dammo*(pgas(k)-CNH3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=IHCL                                           ! do HCL
      hi(i,k)=0.0d0
c
c     if there is no NH4 to evaporate then 
      if(c(i,13).le.tinys.and.hi(i,inh3).lt.0.0d0) goto 119
c     if there is no NO3 to evaporate then 
      if(c(i,13).le.tinys.and.c(i,11).le.tinys.and.
     &   hi(i,ihno3).lt.0.0d0) goto 114

      goto 200
c
c     case 1.1.7  *********************************************
c     Only reaction 2
c     Chloride system (no Nitrates)
 117  icase(i)=117
      p=(-dammo*pgas(inh3)+dhcl*pgas(ihcl)
     #  +2.0d0*dsul*pgas(ih2so4))/(rgas*temp)

      RKNH4CL=xk6/(rgas1*temp)**2

      CHNO3=0.0d0
      CNH3=(-p+sqrt(p**2+4.0d0*dammo*dhcl*RKNH4CL))/(2.0*dammo)
      CHCL=RKNH4CL/CNH3

      k=IHNO3                                          ! do HNO3
      hi(i,k)=0.0d0

      k=INH3                                           ! do NH3
      hi(i,k)=12.0d0*dammo*(pgas(k)-CNH3*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

      k=IHCL                                           ! do HCL
      hi(i,k)=12.0d0*dhcl*(pgas(k)-CHCL*rgas*temp)*
     &        (gmw(k)*1.0d-3)/(1.0d3*dp**2*rgas*temp)

c     if there is no NH4 to evaporate then 
      if(c(i,14).le.tinys.and.hi(i,inh3).lt.0.0d0) goto 118
c     if there is no Cl to evaporate then 
      if(c(i,14).le.tinys.and.c(i,9).le.tinys.and.
     &   hi(i,ihcl).lt.0.0d0) goto 114
      goto 200
c
c     case 1.1.8  *********************************************
c     none of reactions 1, 2, or 3 but JNH3 is set to zero
c     NaCl system
 118  icase(i)=118
      hi(i,INH3)=0.0d0                                 ! do NH3
      hi(i,IHCL)=-2.0d0*hi(i,IH2SO4)/gmw(iH2SO4)*gmw(iHCL) ! HCL
      hi(i,IHNO3)=0.0d0                                ! do HNO3
      goto 200
c
c     case 1.1.9  *********************************************
c     none of reactions 1, 2, or 3 but JNH3 is set to zero
c     NaNO3 system
 119  icase(i)=119
      hi(i,INH3)=0.0D0                                 ! do NH3
      hi(i,IHNO3)=-2.0d0*hi(i,IH2SO4)/gmw(iH2SO4)*gmw(iHNO3) ! HNO3
      hi(i,IHCL)=0.0d0                                 ! do HCL

c 200  ht(i)=0.0d0
 200  ht(i)=ht(i)+hi(i,IH2SO4)+hi(i,IHCL)+hi(i,IHNO3)+hi(i,INH3)
c      write(6,*)'section ',i,'case ',icase(i)
      return
      end
