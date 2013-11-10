c=========================================================================
c  03/09/03: bkoo
c            - modified to eliminate organics (eqparto now deals with organics)
c            - moved call to wdiameter
c            - added call to eqparto
c  03/07/03: bkoo
c            - commented out last call to step (redundant)
c            - brought initial call to step out of IF statement
c              so that step can be called when HYBR is selected
c  12/05/02: tmg
c            - eliminated call to subroutine eqneut
c  10/22/02: tmg
c            - added parameter to subroutine step because not all sections
c              are always in use
c  01/25/02: tmg
c            - use dry and wet basis for diameter as appropriate
c  09/20/00: bkoo
c            - combine eqpart & eqparth
c  06/08/00: bkoo
c            - add coagulation & nucleation
c  Apr 2000: bkoo
c            - put a counter to avoid infinity loop
c            - fix negative concentrations
c            - need codes for organics evaporation
c=========================================================================
c
cgy supress error messages to units 39 and 90 (6/12/00)
c
c  EQILIBRIUM PARTITIONING OF CONDENSING OR EVAPORATING SPECIES
C  ACROSS A SIZE AND COMPOSITION RESOLVED AEROSOL DISTRIBUTION
C
C  THIS SUBROUTINE ASSUMES THAT GAS-AEROSOL EQUILIBRIUM IS ACHIEVED 
C  IN FOR THE SUM OF THE SIZE SECTIONAL COMPOSITIONS.
C
C  WRITTEN BY KEVIN P. CAPALDO
c  DEPARTMENT OF CHEMICAL ENGINEERING
C  CARNEGIE MELLON UNIVERSITY
C  JUNE 1998
C
c   SUBROUTINES REQUIRED
C     ISRPIA        EQUILIBRIUM THERMODYNAMIC SUBROUTINE ISORROPIA
C     ORGANICS      GIVES THE EQUILIBRIUM VAPOR PRESSURE OF ORGANIC SPECIES
C     STEP          CALLS ISRPIA IN REVERSE MODE FOR EACH SECTION
C     
C   VARIABLES
c     T             TIME IN SECONDS
C     Q(NTOTALx2)   ARRAY OF AEROSOL SPECIES FOR EACH SECTION (UG/M3) 
C                   THEN GASES (PPM)
C     DQDT(NTOTAL)  THE AEROSOL AND GAS DERIVATIVES (PER SECOND)
C     DQ(NSP)       THE TOTAL MASS OF EACH SPECIES TRANSFERED 
C                   FROM THE GAS TO THE AEROSOL PHASE (uMOLES/m3)
C          THE FOLLOWING TWO DIMENSIONAL ARRAYS ARE NECESARY FOR THE 
C          LINKING OF THE WEIGHTING FACTORS FOR NH4, NO3 AND CL.
C     FRQ(NSEC,nsp)   FRACTION OF DQ TRANSFERED TO EACH AEROSOL SECTION
C
      SUBROUTINE EQPART(t,q)

      include 'dynamic.inc'
      INCLUDE 'equaer.inc'                  ! ISORROPIA declarations

cbk      REAL*8 DQ(NSP),FRQ(nsec,nsp),WI(5),q0(3),q(ntotal)
cbk      real*8 qsav(ntotal), DQsav(nsp), accom(nsp)
      REAL*8 DQ(nexti),FRQ(nsec,nexti),WI(5),q0(3),q(ntotal) ! bkoo (03/09/03)
      real*8 qsav(nexti*nsec), DQsav(nexti), accom(nexti)    ! bkoo (03/09/03)
      real*8 qq, frq0(nsec,nexti) ! bkoo (10/07/03)
      logical done

      if(aerm.eq.'EQUI') then	
c
c     STEP 1/3: CALCULATE NUCLEATION RATE FOR THE WHOLE STEP
c
       call nucl(q)
c
c     STEP 2/3: CALCULATE COAGULATION RATE FOR THE WHOLE STEP
c
       call coagul(q)

cbk       call step(nsec,q) ! tmg (10/22/02)
      endif
      call step(nsecx2,q) ! bkoo (03/07/03)
      call wdiameter(q) ! bkoo (03/09/03)
      call eqparto(t,q) ! bkoo (03/09/03)
c
C     STEP 1:  DETERMINE BULK EQUILIBRIUM
c
C  INPUT:
c     IPROB = 1=Reverse prob, 0=Foreward prob
C  1. [WI] 
C     DOUBLE PRECISION array of length [NCOMP=5]. [NCOMP] is defined in 
C     include file 'isrpia.inc'.
C     Total aerosol concentrations, expressed in moles/m3.  
C     WI(1) - sodium,   expressed as [Na]
C     WI(2) - sulfate,  expressed as [H2SO4]
C     WI(3) - ammonium, expressed as [NH3]
C     WI(4) - nitrate,  expressed as [HNO3]
C     WI(5) - chloride, expressed as [HCl]
C
C  2. [RHI] 
C     DOUBLE PRECISION variable.  
C     Ambient relative humidity expressed on a (0,1) scale.
C
C  3. [TEMPI]
C     DOUBLE PRECISION variable. 
C     Ambient temperature expressed in Kelvins. 
C
C *** CONVERT INPUT CONCENTRATIONS TO moles/m3 **************************
C
      ng = nsp*nsecx2                       ! gases
      prod=rgas*temp/(1.01325D5*pres)       ! conversion from ppm to umoles/m3
                                            ! pres (bkoo, 06/09/00)
      WI(1) = 0.0D0
      WI(2) = Q(ng+ih2SO4) / PROD*1.0d-6
      WI(3) = q(ng+iNH3)   / PROD*1.0d-6
      WI(4) = q(ng+ihNO3)  / PROD*1.0d-6
      WI(5) = q(ng+ihCL)   / PROD*1.0d-6
      do k=1,nsecx2                         ! aerosols
         nn=(k-1)*nsp
         WI(1) = wi(1)+ q(nn+KNa) /emw(KNa) *1.D-6
         WI(2) = wi(2)+ q(nn+KSO4)/emw(KSO4)*1.D-6
         WI(3) = wi(3)+ q(nn+KNH4)/emw(KNH4)*1.D-6
         WI(4) = wi(4)+ q(nn+KNO3)/emw(KNO3)*1.D-6
         WI(5) = wi(5)+ q(nn+KCL) /emw(KCL) *1.D-6
      enddo
      RHI=rh
      TEMPI=temp
      IPROB = 0
C
C *** CALL ISORROPIA+ ***************************************************
C
      CALL ISRPIA ( WI, RHI, TEMPI, IPROB )

c     initialize DQ & ACCOM arrays - bkoo (05/24/01)
cbk      do i=1,nsp
      do i=1,nexti ! bkoo (03/09/03)
         dq(i)    = 0.0d0
         accom(i) = 1.0d0
      enddo
C
C   STEP 2:  DETERMINE THE TOTAL TRANSFER BETWEEN GAS AND AEROSOL
C
c      DQ(KH2O)=0.0D0                        ! WATER DONE LATER
c      DQ(KNA) =0.0D0                        ! SODIUM IS IN AEROSOL PHASE ONLY
      DQ(KSO4)=Q(ng+ih2so4)/prod            ! all H2SO4 condenses
      DQ(KNH4)=Q(NG+INH3)  /PROD - dmax1(GNH3 *1.0d6, 0.d0) ! bkoo (10/07/03)
      DQ(KNO3)=Q(NG+IHNO3) /PROD - dmax1(GHNO3*1.0d6, 0.d0) ! bkoo (10/07/03)
      DQ(KCL) =Q(NG+IHCL)  /PROD - dmax1(GHCL *1.0d6, 0.d0) ! bkoo (10/07/03)
C    DO ORGANICS  **ASSUME PS HAS NO SIZE OR COMPOSITION DEPENDANCE **
cbk  removed - bkoo (03/09/03)
cbk      DO IOG = 1,NORG-1
cbk         DQ(NEXTI+IOG)=(Q(NG+IHCL+IOG)-PS(1,IHCL+IOG)/
cbk     &                 (1.01325d-1*pres))/prod	! pres (bkoo, 06/09/00)
cbk      ENDDO
c
c   STEP 2B: PARTITION DQ's NH4, NO3, AND CL INTO NH4NO3, NH4CL, AND
C            EXCESS NH4 (which is associated with SO4)
C            DQ(KNO3) => now represents the transfer of NH4NO3
C            DQ(KCl) => now represents the transfer of NH4Cl
C            DQ(KNH4) => now represents the transfer of NH4 with SO4
c
      DQ(KNH4)=DQ(KNH4)-DQ(KNO3)-DQ(KCL)
c
c     Save initial aerosol concentrations
cbk      do i=1,ntotalx2
cbk        qsav(i)=q(i)
      do i=1,nsecx2                          ! bkoo (03/09/03)
       do ii=1,nexti                         ! bkoo (03/09/03)
        qsav((i-1)*nexti+ii)=q((i-1)*nsp+ii) ! bkoo (03/09/03)
       enddo                                 ! bkoo (03/09/03)
      enddo
c
c     Its possible to not be able to reach equilibrium for NH4Cl and NH4NO3
c     due to diferences between aerosol sectional vs. bulk composition 
      do i=1,3
        q0(i)=0.0d0
      enddo
      do isec=1,nsecx2
        q0(1)=q0(1)+q((isec-1)*nsp+kno3)    ! initial aerosol NO3
        q0(2)=q0(2)+q((isec-1)*nsp+knh4)    ! initial aerosol NH4
        q0(3)=q0(3)+q((isec-1)*nsp+kcl)     ! initial aerosol Cl
      enddo
C
c   STEP 3: DETERMINE THE RELATIVE RATES OF MASS TRANSFER FOR EACH SECTION
c       if we assume composition changes between size sections are not 
c       important, then the rate of transfer is proportional to the mass
c       mass transfer rate dependance on particle size

cbk   removed - bkoo (03/09/03)
cbk      accom(KH2O)=1.0
cbk      accom(KNa)=1.0
      accom(KSO4)=delta(ih2so4)
      accom(kno3)=delta(ihno3)
      accom(knh4)=delta(inh3)
      accom(kcl) =delta(ihcl)
cbk   removed - bkoo (03/09/03)
cbk      do isp=1,norg
cbk        accom(kcl+isp)=delta(ihcl+isp)
cbk      enddo
cbk      do isp=1,ninert
cbk        accom(kcl+norg+isp)=0.1
cbk      enddo

cbk      call wdiameter(q) ! tmg (01/25/02)

      rlambda=0.065d0
c
c calculate factors
c
cbk      do isp=1,nsp
      do isp=1,nexti ! bkoo (03/09/03)
       frqtot = 0.0
       do isec = 1,nsecx2
        frq(isec,isp) = qn(isec)*
     &      dsec(isec)/(1.0+rlambda/(accom(isp)*dsec(isec)))
        frqtot = frqtot + frq(isec,isp)
       enddo
c
c normalize
c
       do isec = 1,nsecx2
        frq(isec,isp) = frq(isec,isp)/frqtot
        frq0(isec,isp) = frq(isec,isp) ! save frq - bkoo (10/07/03)
       enddo
      enddo

      iter=0 ! counter for escaping out of an infinity loop (bkoo: Apr, 2000)
 80   iter=iter+1 ! moved - bkoo (10/07/03)
c     save dq's
cbk      do i=1,nsp
      do i=1,nexti ! bkoo (03/09/03)
        dqsav(i)=dq(i)
      enddo
c
c   FIRST condense all condensing species
c
cbk      do isp = 1,nsp
      do isp = 1,nexti ! bkoo (03/09/03)
       if(dq(isp).gt.0.0d0) then
        do isec=1,nsecx2
         INDX=(ISEC-1)*NSP
         q(indx+isp)=q(indx+isp)+frq(isec,isp)*dq(isp)*emw(isp)
c     special cases for remaping of dq
         if(isp.eq.kNO3.or.isp.eq.kCl) then
          q(indx+kNH4)= q(indx+knh4)+frq(isec,isp)*dq(isp)*emw(knh4)
         endif
        enddo
        dq(isp)=0.0d0
       endif
      enddo
c
c   SECOND evaporate all evaporating species
c   only NH4NO3 and NH4Cl can evaporate
c
C     CHECK FOR COMPLETE EVAPORATION
 100  frt=1.0d0
      frtcl = 0.0 ! bkoo (10/07/03)
      isp=kCl ! Cl first b/c NH4Cl forms before NH4NO3 when Na exists
      do m=1,2
       if(dq(isp).lt.0.0d0) then            ! evaporating
        do isec=1,nsecx2
         INDX=(ISEC-1)*NSP
         dqfx=DQ(ISP)*FRQ(ISEC,isp)
         IF(-dqfx.gt.tinys) then            ! evaporating significantly
          if(Q(INDX+ISP).LT.-dqfx*emw(isp)) then ! not enough NO3 or Cl
           frtq=-q(indx+isp)/dqfx/emw(isp)
           if(frtq.lt.frt) then
            frt=frtq
            ispsav=isp 	                    ! species that evaporates first
            isecsav=isec                    ! section that evaporates first
            if(q(indx+isp).lt.tinys) goto 150
           endif
          endif
          qq = q(indx+knh4) ! bkoo (10/07/03)
          if(isp.eq.KNO3)   ! bkoo (10/07/03)
     &              qq = qq + frtcl*frq(isec,KCL)*dq(KCL)*emw(knh4)
cbk          if(q(indx+knh4).lt.-dqfx*emw(kNH4)) then ! not enough NH4+
cbk           frtq=-q(indx+knh4)/dqfx/emw(kNH4)
          if(qq.lt.-dqfx*emw(kNH4)) then    ! not enough NH4+
           frtq=-qq/dqfx/emw(kNH4)
           if(frtq.lt.frt) then
            frt=frtq
            ispsav=isp                      ! species that evaporates first
            isecsav=isec                    ! section that evaporates first
            if(q(indx+knh4).lt.tinys) goto 150
           endif
          endif
         endif
        enddo
        if(isp.eq.KCL) frtcl = frt ! bkoo (10/07/03)
       endif
       isp=kNO3
      enddo
c partition species up to evaporation point
      done=.true.
      do m=1,2
       do isec=1,nsecx2
        INDX=(ISEC-1)*NSP
        q(indx+isp)=q(indx+isp)+frt*frq(isec,isp)*dq(isp)*emw(isp)
        q(indx+kNH4)=q(indx+knh4)+frt*frq(isec,isp)*dq(isp)*emw(knh4)
       enddo
       dq(isp)=(1-frt)*dq(isp)
       if(abs(dq(isp)).gt.1e-5) done=.false. ! tolerance for solution
       isp=kCl
      enddo
c check for complete solution
      if(done) then
       goto 200
      endif
c
c adjust frq to account for the totally evaporated species
c
 150  frqtot=0.0 ! bkoo (10/07/03) ...
      do isec=1,nsecx2
        if(isec.eq.isecsav) frq(isec,ispsav)=0.0d0
        frqtot = frqtot + frq(isec,ispsav)
      enddo
      if(frqtot.gt.tinys) then ! renormalize
        do isec=1,nsecx2
          frq(isec,ispsav)=frq(isec,ispsav)/frqtot
        enddo
      else
        if(iter.gt.itmaxeq) then ! moved - bkoo (10/07/03)
         if(min(dq(KNO3),dq(KCL)).lt.-0.3) then                               ! bkoo_dbg
         call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)                    ! bkoo_dbg
         write(*,'(A8,2E15.5,4I4)')'EQUI-F: ',dq(KNO3),dq(KCL),ichm,jchm,kchm ! bkoo_dbg
         endif                                                                ! bkoo_dbg
         goto 200
        endif
        goto 180 ! can't evaporate given species from any section
      endif

      goto 100
c
c   step 4: when sectional compositions prevent achievement of equilibrium
c           gas phase concentrations we need to adjust the non limiting
c           species so that that species does achieve equilibrium.
c           For example:  if NH4Cl cannot evaporate because a section 
c           runs out of NH4 and this is the only particle with any Cl
c           then KNH4Cl will not equal [NH3]final [HCL]final becuase 
c           this equilibrium could not be achieved.  However the 
c           equilibrium transfer flux of NH4NO3 was calculated assuming 
c           that all of the NH4Cl would evaporate.  Since it didn't 
c           KNH4NO3 will not equal [NH3]final [HNO3]final.  To correct
c           in this case we will evaporate enough NH4NO3 to establish
c           equilibrium. 
c
c           only evaporating species can encounter this problem
c           only if NH4 becomes zero in a section can this problem occur
c
 180  if(ispsav.eq.kno3) then
         isp2=kcl                     ! correcting species
         g2=ghcl*1.0d6                ! original equilibrium for HCL (umol/m3)
      else
         isp2=kno3                    ! correcting species
         g2=ghno3*1.0d6               ! original equilibrium for HNO3 (umol/m3)
      endif
      bb=gnh3*1.0d6+dq(ispsav)+g2
      dq(isp2)=(-bb+sqrt(bb**2-4*dq(ispsav)*g2))/2.0d0

      dq(isp2) = dmax1(dq(isp2), dq(knh4)+dqsav(ispsav)-dq(ispsav) ! bkoo (10/07/03)
     &                          +dqsav(isp2)-q(ng+inh3)/prod)

      if(iter.gt.itmaxeq+1) dq(isp2)=0.0d0 ! bkoo (11/14/01)
c reset aerosol concentrations and dq's for NO3, NH4 and Cl
      do isec=1,nsecx2
         indx=(isec-1)*nsp
cbk         q(indx+kno3)=qsav(indx+kno3) ! initial aerosol no3
cbk         q(indx+knh4)=qsav(indx+knh4) ! initial aerosol nh4
cbk         q(indx+kcl)=qsav(indx+kcl)   ! initial aerosol Cl
         indx2=(isec-1)*nexti          ! bkoo (03/09/03)
         q(indx+kno3)=qsav(indx2+kno3) ! bkoo (03/09/03)
         q(indx+knh4)=qsav(indx2+knh4) ! bkoo (03/09/03)
         q(indx+kcl)=qsav(indx2+kcl)   ! bkoo (03/09/03)
      enddo
c restore frq - bkoo (10/07/03)
      do isp=1,nexti
        do isec=1,nsecx2
          frq(isec,isp) = frq0(isec,isp)
        enddo
      enddo
c adjust DQ
cbk      do i=1,nsp
      do i=1,nexti ! bkoo (03/09/03)
         if(i.eq.ispsav) then
            dq(i)=dqsav(i)-dq(i)
         elseif(i.eq.isp2) then
            dq(i)=dqsav(i)-dq(i)
         elseif(i.eq.knh4) then
            dq(i)=dqsav(i)
         else
            dq(i)=0.0d0
         endif
      enddo
      goto 80
c      endif
C      
C   STEP 4: ASSIGN EQUILIBRIUM VALUES TO GASES
c   Its possible to not be able to reach equilibrium for NH4Cl and NH4NO3
c   do to diferences between aerosol sectional vs. bulk composition 
C   So we treat NH3, HNO3, and HCl conservatively
C
 200  Q(NG+IH2SO4)=0.0D0              ! all sulfuric acid condenses
      Q(NG+IHNO3)= Q(NG+IHNO3)+q0(1)*PROD/GMW(IHNO3)
      Q(NG+INH3) = Q(NG+INH3) +q0(2)*PROD/GMW(INH3)
      Q(NG+IHCL) = Q(NG+IHCL) +q0(3)*PROD/GMW(IHCL)
      do isec=1,nsecx2
       Q(NG+IHNO3)= Q(NG+IHNO3)-Q((ISEC-1)*NSP+KNO3)*PROD/GMW(IHNO3)
       Q(NG+INH3) = Q(NG+INH3) -Q((ISEC-1)*NSP+KNH4)*PROD/GMW(INH3)
       Q(NG+IHCL) = Q(NG+IHCL) -Q((ISEC-1)*NSP+KCL) *PROD/GMW(IHCL)
      enddo
c     correct negative NH3 - bkoo (10/07/03)
      dq(KNH4) = q(ng+inh3) / prod ! umol/m3
      if(dq(KNH4).lt.-tinys) then
        iter = 0
 300    frt = 1.d0
        do isec = 1,nsecx2
          indx=(isec-1)*nsp       
          if(frq0(isec,KNH4).gt.0.d0) frt = dmax1( dmin1(q(indx+KNH4)
     &             /(-dq(KNH4)*frq0(isec,KNH4)*emw(KNH4)),frt), 0.d0)
        enddo
        frqtot = 0.d0
        do isec = 1,nsecx2
          indx=(isec-1)*nsp
          q(indx+KNH4) = dmax1(q(indx+KNH4)
     &                  +frt*dq(KNH4)*frq0(isec,KNH4)*emw(KNH4),0.d0)
          if(q(indx+KNH4).lt.tinys) frq0(isec,KNH4) = 0.d0
          frqtot = frqtot + frq0(isec,KNH4)
        enddo
        q(ng+inh3) = q(ng+inh3) - frt*dq(KNH4)*prod
        ! check if we should evaporate more
        dq(KNH4) = (1.d0 - frt) * dq(KNH4)
        if(dq(KNH4).lt.-tinys) then
          if(frqtot.gt.tinys) then ! we have sections which are not empty    
            if(iter.le.itmaxeq) then ! check infinite loop
              iter = iter + 1
              do isec = 1,nsecx2
                frq0(isec,KNH4) = frq0(isec,KNH4) / frqtot
              enddo
              goto 300
            endif
          endif
          ! we need to evaporate more to achieve equilibrium
          ! but we completely evaporate the species in all sections
          ! or exceed itermax
          write(*,*)'EQUI-F2: dq(KNH4)=',dq(KNH4),' iter=',iter ! bkoo_dbg
        endif
      endif
C   DO ORGANICS  **ASSUME PS HAS NO COMPOSITION DEPENDANCE **
cbk removed - bkoo (03/09/03)
cbk      DO IOG = 1,NORG-1
cbk       Q(NG+IHCL+IOG) =PS(1,IHCL+NORG)/(1.01325d-1*pres) ! pres (bkoo, 06/09/00)
cbk      ENDDO
  
      call negchk(t,q,nsecx2)

      call eqneut(0, q) ! tmg (12/05/02)

      if(aerm.eq.'EQUI') then ! (tmg,01/31/02)
C
C   STEP 5: CALL EQUILIBRIUM CODE FOR EACH SECTION TO DETERMINE WATER
C
cbk        CALL STEP(nsec,Q)  ! tmg (10/22/02) bkoo (03/07/03)
C
C   STEP 6: Determine new diameter
C
        call ddiameter(q)
      endif
C
      RETURN
      END     

c ----------------------------------------------------------------------------
c neutralize the acidic or basicness of the aerosol in each section
c
c is   : 0 (false) for equilibrium section
c        1 (true)  for dynamic     section
c q    : aerosol [ug/m3] and gas [ppm] species array
c        dimension - ntotal  for equilibrium section of EQUI &
c                                dynamic     section of MADM
c                    ntotale for equilibrium section of HYBR
c                    ntotald for dynamic     section of HYBR
c
c CHTOT: Total charge balance [umole/m3]
c
      subroutine eqneut(is, q)

      include 'dynamic.inc'

      real*8      q(ntotal),chtot,chmove,chno3,chcl,chnh4,prod
      integer     nsx,ng,na,chspc(nexti)
      integer     jani,jcat ! bkoo (10/07/03)
      integer     kerr ! bkoo (01/06/04)

c charges of H2O, Na, SO4, NO3, NH4, Cl
      data chspc / 0, 1, -2, -1, 1, -1 /

      nsx  = nsecx * is + nsecx2 * (1 - is)
      ng   = nsp * nsx
      prod = rgas*temp/(1.01325d5*pres) ! PPM / prod = UMOLE/M3

c      if(is) call negchk(0.d0,q,nsecx)

      do i = 1, nsx
        na = (i-1)*nsp
        chtot = 0.0d0
        jani = 0 ! bkoo (10/07/03)
        jcat = 0 ! bkoo (10/07/03)
        do j = 2, nexti ! the first element is H2O
          chtot = chtot + DBLE(chspc(j)) * q(na+j) / emw(j)
          if(chspc(j).lt.0 .and. q(na+j).gt.tinys) jani = 1 ! bkoo (10/07/03)
          if(chspc(j).gt.0 .and. q(na+j).gt.tinys) jcat = 1 ! bkoo (10/07/03)
        enddo
c BASIC:
        if(chtot .gt. 0.0d0 .and. jani .eq. 0) then ! bkoo (10/07/03)
          chno3 = (q(ng+IHNO3)/prod)
          chcl  = (q(ng+IHCL) /prod)
          chnh4 = (q(na+KNH4) /emw(KNH4))
c first move HCl in
          chmove     = max(min(chtot, chcl),0.0d0)
          q(na+KCL)  = q(na+KCL)  + chmove * emw(KCL)
          q(ng+IHCL) = max(q(ng+IHCL) - chmove * prod,0.0d0)
          chtot      = chtot - chmove
c then move HNO3 in
          chmove     = max(min(chtot, chno3),0.0d0)
          q(na+KNO3) = q(na+KNO3) + chmove * emw(KNO3)
          q(ng+IHNO3)= max(q(ng+IHNO3)- chmove * prod,0.0d0)
          chtot      = chtot - chmove
c then move NH4 out
          chmove     = max(min(chtot, chnh4),0.0d0)
          q(na+KNH4) = max(q(na+KNH4) - chmove * emw(KNH4),0.0d0)
          q(ng+INH3) = q(ng+INH3) + chmove * prod
          chtot      = chtot - chmove
        endif                                       ! bkoo (10/07/03)
c ACIDIC:
        if(chtot .lt. 0.0d0 .and. jcat .eq. 0) then ! bkoo (10/07/03)
          chtot = -1.0d0 * chtot
          chno3 = (q(na+KNO3)/emw(KNO3))
          chcl  = (q(na+KCL) /emw(KCL))
          chnh4 = (q(ng+INH3)/prod)
c first move NH3 in
          chmove     = max(min(chtot, chnh4),0.0d0)
          q(na+KNH4) = q(na+KNH4) + chmove * emw(KNH4)
          q(ng+INH3) = max(q(ng+INH3) - chmove * prod,0.0d0)
          chtot      = chtot - chmove
c then move HCl out
          chmove     = max(min(chtot, chcl),0.0d0)
          q(na+KCL)  = max(q(na+KCL)  - chmove * emw(KCL),0.0d0)
          q(ng+IHCL) = q(ng+IHCL) + chmove * prod
          chtot      = chtot - chmove
c then move HNO3 out
          chmove     = max(min(chtot, chno3),0.0d0)
          q(na+KNO3) = max(q(na+KNO3) - chmove * emw(KNO3),0.0d0)
          q(ng+IHNO3)= q(ng+IHNO3)+ chmove * prod
          chtot      = chtot - chmove
        endif
      enddo
c
c     check the charge balance in each section - bkoo (01/06/04)
c
      do i = 1, nsx
        na = (i-1)*nsp
        kerr = 0
        if (q(na+KNA).gt.tinys .or. q(na+KNH4).gt.tinys) then
          if (q(na+KSO4).le.tinys .and. q(na+KNO3).le.tinys .and.
     &        q(na+KCL).le.tinys) kerr = 1
cbk          if ( (2.*q(na+KSO4)/emw(KSO4)+q(na+KNO3)/emw(KNO3)+q(na+KCL)
cbk     &            /emw(KCL))/(q(na+KNA)/emw(KNA)+q(na+KNH4)/emw(KNH4))
cbk     &            .gt.0.2e10) kerr = 1
        else
          if (q(na+KNO3).gt.tinys .or. q(na+KCL).gt.tinys) kerr = 1
        endif
        if (kerr.ne.0) then
          call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)
ctmg          write(iout,'(//,A)')'CHECK in EQNEUT'
ctmg          write(iout,*)' q(',i,'): ',(q(na+j),j=2,nexti)
ctmg          write(iout,*)' igrd,i,j,k: ',igrdchm,ichm,jchm,kchm
        endif
      enddo

      return
      end

