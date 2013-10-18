C=======================================================================
C
C *** MADM CODE
C *** SUBROUTINE EQUAER 
C *** THIS SUBROUTINE CALLS ISORROPIA+ FOR CALCULATING THE AEROSOL 
C     EQUILIBRIUM CONCENTRATIONS 
C
C ======================== ARGUMENTS / USAGE ===========================
C
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
C  OUTPUT:
C  1. [AEROSI]
C     DOUBLE PRECISION array of length [NINTI=21]. 
C     Aerosol species concentrations, expressed in ug/m3. 
C     c(K,01) - H2O             c(K,12) - (NH4)2SO4(s)
C     c(K,02) - H+(aq)          c(K,13) - NH4NO3(s)
C     c(K,03) - Na+(aq)         c(K,14) - NH4Cl(s)
C     c(K,04) - NH4+(aq)        c(K,15) - H2SO4(aq) (Undefined)
C     c(K,05) - Cl-(aq)         c(K,16) - NH4HSO4(s)
C     c(K,06) - SO4--(aq)       c(K,17) - NaHSO4(s)
C     c(K,07) - HSO4-(aq)       c(K,18) - (NH4)4H(SO4)2(s)
C     c(K,08) - NO3-(aq)        c(K,19( - OH-(aq)
C     c(K,09) - NaCl(s)         
C     c(K,10) - Na2SO4(s)       
C     c(K,11) - NaNO3(s)
C
C  2. [ps]
C     DOUBLE PRECISION array of length [NGAS]. 
C     Gaseous species partial pressures, expressed in Pascals.
C     ps(k,1) - NH3
C     ps(k,2) - HNO3
C     ps(k,3) - H2SO4
C     ps(k,4) - HCl 
C
C  3. [DRYF]
C     LOGICAL variable.
C     Contains information about the physical state of the system.
C     .TRUE. - There is no aqueous phase present
C     .FALSE.- There is an aqueous phase present
C 
C=======================================================================
C
      SUBROUTINE EQUAER (k, q, iprob)
C
      INCLUDE 'dynamic.inc'      ! MADM declarations
      INCLUDE 'equaer.inc'       ! ISORROPIA declarations
C
      real*8 q(ntotal)           ! local conc. array
      DIMENSION WI(NCOMP)
      dimension pgas(ngas)
      integer k                  ! master array section indices
      integer nsec0              ! last section of master array NOT in q
      integer nsec1              ! last section of master array occuring in q
C
C *** CONVERT INPUT CONCENTRATIONS TO moles/m3 **************************
C
      nn = (k-1)*nsp             ! aerosol
      ng = nsp*nsec              ! gases
      if(q(nn+KNa).le.tinys.and.q(nn+KSO4).le.tinys.and.
     #   q(nn+KNH4).le.tinys.and.q(nn+KNO3).le.tinys.and.
     #   q(nn+KCL).le.tinys) go to 999
 100  WI(1) = q(nn+KNa) /emw(KNa) *1.D-6
      WI(2) = q(nn+KSO4)/emw(KSO4)*1.D-6
      WI(3) = q(nn+KNH4)/emw(KNH4)*1.D-6
      WI(4) = q(nn+KNO3)/emw(KNO3)*1.D-6
      WI(5) = q(nn+KCL) /emw(KCL) *1.D-6
c iprob is always 1; the codes below are never used - bkoo
      if(iprob.eq.0) then        ! add gas phase
       prod=rgas*temp/(1.01325D5*pres)    ! conversion from ppm to umoles/m3
                                          ! pres (bkoo, 06/09/00)
       WI(1) = WI(1) + 0.0D0
       WI(2) = WI(2) + Q(ng+ih2SO4)/PROD*1.0d-6
       WI(3) = WI(3) + q(ng+iNH3)  /PROD*1.0d-6
       WI(4) = WI(4) + q(ng+ihNO3) /PROD*1.0d-6
       WI(5) = WI(5) + q(ng+ihCL)  /PROD*1.0d-6
      endif
C
C *** CALL ISORROPIA+ ***************************************************
C
      RHI=rh
      TEMPI=temp
c      IPROB = 1                  ! 1=Reverse prob, 0=Foreward prob

      metstbl=ims(k) ! bkoo (03/05/02)
      CALL ISRPIA ( WI, RHI, TEMPI, IPROB )
      metstbl=0      ! bkoo (03/05/02)
C
C *** CONVERT OUTPUT CONCENTRATIONS TO ug/m3 ****************************
C
C AEROSOL SPECIES
C
      q(nn+KH2O)= WATER *1.D9                  !  H2O
      C(k,01)   = q(nn+KH2O)                   !  H2O
      C(k,02)   = MOLAL(01)          *1.D6     !  H+(aq)          
      C(k,03)   = MOLAL(02)*intmw(03)*1.D6     !  Na+(aq)         
      C(k,04)   = MOLAL(03)*intmw(04)*1.D6     !  NH4+(aq)        
      C(k,05)   = MOLAL(04)*intmw(05)*1.D6     !  Cl-(aq)         
      C(k,06)   = MOLAL(05)*intmw(06)*1.D6     !  SO4--(aq)
      C(k,07)   = MOLAL(06)*intmw(07)*1.D6     !  HSO4-(aq)
      C(k,08)   = MOLAL(07)*intmw(08)*1.D6     !  NO3-(aq) 
      C(k,09)   = CNACL    *intmw(09)*1.D6     !  NaCl(s)  
      C(k,10)   = CNA2SO4  *intmw(10)*1.D6     !  Na2SO4(s)
      C(k,11)   = CNANO3   *intmw(11)*1.D6     !  NaNO3(s)
      C(k,12)   = CNH42S4  *intmw(12)*1.D6     !  (NH4)2SO4(s)
      C(k,13)   = CNH4NO3  *intmw(13)*1.D6     !  NH4NO3(s)
      C(k,14)   = CNH4CL   *intmw(14)*1.D6     !  NH4Cl(s)
      C(k,15)   = 0.0D0                    !  H2SO4(aq) (not used by isorropia)
      C(k,16)   = CNH4HS4  *intmw(16)*1.D6     !  NH4HSO4(s)      
      C(k,17)   = CNAHSO4  *intmw(17)*1.D6     !  NaHSO4(s)
      C(k,18)   = CLC      *intmw(18)*1.D6     !  (NH4)4H(SO4)2(s)
      C(k,19)   = COH      *intmw(19)*1.D6 !  bkoo (02/14/02)
c
c    Since issorropia ignores no3 and cl concentrations if the particle is 
c    found to be both acidic and dry, we need to evaporate these species
c    instantaneously.  Since NH4Cl forms from H2SO4 condensation before
c    NH4NO3, we start transfer cloride before nitrate.  
c    We then call issorropia again to calculate the new ci's.
c
c    Since statements except [if] and [write] have been blocked, however,
c    nothing actually happens. Therefore the whole part was made blocked.
c    (if unblocked, it should be modified so that pressure terms are included)
c    - bkoo
c 
c      if(dryf) then
c       if(c(k,18).gt.tinys.or.c(k,17).gt.tinys.or.c(k,16).gt.tinys) then
c        irecalc=0
c        acid=clc*1.d6+cnahso4*1.d6+cnh4hs4*1.d6
c        if(q(nn+kcl).gt.tinys) then
c         if(q(nn+kcl).lt.acid*emw(kcl)) then
c          write(6,*) tcom,'SOLID ACID ',k,q(nn+kcl),
c     &                                       'ug/m3 cl- moved to gas1'
c          write(90,*) tcom,'SOLID ACID ',k,q(nn+kcl),
c     &                                       'ug/m3 cl- moved to gas1'
c          q(ng+ihcl)=q(ng+ihcl)+
c     &                         q(nn+kcl)*rgas*temp/(1.015d5*gmw(ihcl))
c          q(nn+kcl)=0.0d0
c          acid=acid-q(nn+kcl)/float(emw(kcl))
c         else
c          write(6,*) tcom,'SOLID ACID ',k,acid*emw(kcl),
c     &                                  'ug/m3 cl- moved to gas2'
c          write(90,*) tcom,'SOLID ACID ',k,acid*emw(kcl),
c     &                                  'ug/m3 cl- moved to gas2'
c          q(ng+ihcl)=q(ng+ihcl)+
c     &                 acid*emw(kcl)*rgas*temp/(1.015d5*gmw(ihcl))
c          q(nn+kcl)=q(nn+kcl)-acid*emw(kcl)
c          acid=0.0d0
c         endif
c         irecalc=1
c        endif
c        if(q(nn+kno3).gt.tinys.and.acid.gt.0.0d0) then
c         if(q(nn+kno3).lt.acid*emw(kno3)) then
c          write(6,*) tcom,'SOLID ACID ',k,q(nn+kno3),
c     &                                      'ug/m3 NO3- moved to gas1'
c          write(90,*) tcom,'SOLID ACID ',k,q(nn+kno3),
c     &                                      'ug/m3 NO3- moved to gas1'
c          q(ng+ihno3)=q(ng+ihno3)+
c     &                         q(nn+kno3)*rgas*temp/(1.015d5*gmw(ihno3))
c          q(nn+kno3)=0.0d0
c         else
c          write(6,*) tcom,'SOLID ACID ',k,acid*emw(kno3),
c     &                                  'ug/m3 NO3- moved to gas2'
c          write(90,*) tcom,'SOLID ACID ',k,acid*emw(kno3),
c     &                                  'ug/m3 NO3- moved to gas2'
c          q(ng+ihno3)=q(ng+ihno3)+
c     &                 acid*emw(kno3)*rgas*temp/(1.015d5*gmw(ihno3))
c          q(nn+kno3)=q(nn+kno3)-acid*emw(kno3)
c         endif
c         irecalc=1
c        endif
c        if(irecalc.eq.1) goto 100
c       endif
c      endif

C
C GASEOUS SPECIES
C
c Originally Pilinis tried to slow things down by capping the "ps" values
c at 1.0. I had to add the acidity constraint to achieve control in a more
c robust way. Because the time scales of concern for atmospheric modeling
c are large enough, artificially slowing things down to the degree in HYBR
c seemed appropriate to improve stability. There is nothing special however
c about a "ps" of 1.0 though. It just seemed large enough not to seriously
c affect accuracy and small enough to allow the code to run stably in a
c reasonable amount of time. For MADM I expect I took this out because
c I was more concerned about accuracy and less concerned about runtime.
c Also, MADM is somewhat more stable since there is no division between
c "equilibrium" and "dynamic" particles. - Kevin, 5/30/00
c (bkoo, 06/06/00)
c      if(aerm.eq.'MADM') then ! bkoo (01/08/01)
         ps(k,INH3  ) = GNH3 *rgas*temp	!  NH3(g)	(Pa)
         ps(k,IHNO3 ) = GHNO3*rgas*temp	!  HNO3(g)	(Pa)
         ps(k,IH2SO4) = 0.0D0		!  H2SO4(g)	(Pa)
         ps(k,IHCL  ) = GHCL *rgas*temp	!  HCL(g)	(Pa)
c      else
c         ps(k,INH3  ) = min(1.0d0,GNH3 *rgas*temp)
c         ps(k,IHNO3 ) = min(1.0d0,GHNO3*rgas*temp)
c         ps(k,IH2SO4) = 0.0D0
c         ps(k,IHCL  ) = min(1.0d0,GHCL *rgas*temp)
c      endif
C
C 'DRY' FLAG
C
      dry(k) = DRYF
c
      RETURN


999   CONTINUE
      q(nn+KH2O)= 0.0                   !  H2O
      C(k,01)   = 0.0                   !  H2O
      C(k,02)   = 0.0D0                 !  H+(aq)          
      C(k,03)   = 0.0D0                 !  Na+(aq)         
      C(k,04)   = 0.0D0                 !  NH4+(aq)        
      C(k,05)   = 0.0D0                 !  Cl-(aq)         
      C(k,06)   = 0.0D0                 !  SO4--(aq)
      C(k,07)   = 0.0D0                 !  HSO4-(aq)
      C(k,08)   = 0.0D0                 !  NO3-(aq) 
      C(k,09)   = 0.0D0                 !  NaCl(s)  
      C(k,10)   = 0.0D0                 !  Na2SO4(s)
      C(k,11)   = 0.0D0                 !  NaNO3(s)
      C(k,12)   = 0.0D0                 !  (NH4)2SO4(s)
      C(k,13)   = 0.0D0                 !  NH4NO3(s)
      C(k,14)   = 0.0D0                 !  NH4Cl(s)
      C(k,15)   = 0.0D0             !  H2SO4(aq) (not used by isorropia)
      C(k,16)   = 0.0D0                 !  NH4HSO4(s)      
      C(k,17)   = 0.0D0                 !  NaHSO4(s)
      C(k,18)   = 0.0D0                 !  (NH4)4H(SO4)2(s)
      C(k,19)      = 0.0D0          !  bkoo (02/14/02)
      ps(k,INH3  ) = 0.0D0              !  NH3(g)    (Pa)
      ps(k,IHNO3 ) = 0.0D0              !  HNO3(g)   (Pa)        
      ps(k,IH2SO4) = 0.0D0              !  H2SO4(g)  (Pa)        
      ps(k,IHCL  ) = 0.0D0              !  HCL(g)    (Pa)

      RETURN
C
C *** END OF SUBROUTINE EQUAER ******************************************
C
      END
