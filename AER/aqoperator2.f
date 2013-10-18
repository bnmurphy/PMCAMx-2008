      subroutine aqoperator2(tbeg,deltat,gas,aerosol,lwc,
     &  rh,temp,p,iaq)

c
c       TOTAL                        GAS             AQUEOUS   
c-----------------------------------------------------------------------
c     1. FORMALDEHYDE TOTAL        1. SO2       1.  S(IV)
c     2. FORMIC ACID TOTAL         2. H2O2      2.  H2O2(aq)
c                                  3. NH3       3.  NITRATE
c                                               4.  CHLORIDE
c                                               5.  AMMONIUM
c                                               6.  SULFATE
c                                               7.  HSO5-
c                                               8.  HMSA
c                                              
c
c THE AQUEOUS-PHASE EQUATIONS ARE SOLVED FOR EACH SPECIES
c
c
c.....Y MATRIX STRUCTURE  EVERYTHING IS (ug/m3)
c     ONLY 21 ODEs ARE SOLVED FOR THE WHOLE DISTRIBUTION
c
c  YAQ(1) = TOTAL FORMALDEHYDE
c 
c  YAQ(2) = SO2 (g)
c  YAQ(3) = H2O2 (g)
c  YAQ(4) = NH3 (g)
c  YAQ(5) = HCOOH (total)
c
c  YAQ(6) = SO2 (aq) SECTION 1
c  YAQ(7) = H2O2 (aq) SECTION 1
c  YAQ(8) = NITRATE (aq) SECTION 1
c  YAQ(9) = CHLORIDE (aq) SECTION 1        
c  YAQ(10) = AMMONIUM (aq) SECTION 1
c  YAQ(11) = SULFATE (aq) SECTION 1
c  YAQ(12) = HSO5- (aq) SECTION 1
c  YAQ(13) = HMSA (aq) SECTION 1
c
c  YAQ(14) = SO2 (aq) SECTION 2
c  YAQ(15) = H2O2 (aq) SECTION 2
c  YAQ(16) = NITRATE (aq) SECTION 2
c  YAQ(17) = CHLORIDE (aq) SECTION 2        
c  YAQ(18) = AMMONIUM (aq) SECTION 2
c  YAQ(19) = SULFATE (aq) SECTION 2
c  YAQ(20) = HSO5- (aq) SECTION 2
c  YAQ(21) = HMSA (aq) SECTION 2
c
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
      real*4 gas(ngas_aq), aerosol(nsect,naers)    
      real*4 sodium(nsect)
      real*4 deltat, rh, temp
      real*4 gascon(ngas_aq)
      real*4 yaq(meqn2)
      real*4 watcont(2)
      real*4 salt(2), crustal(2), lwc(2)
      real*4 check(3,naers)
c
      common / gascon / gascon
      common / epcomy / ymin, hmax
      common / integration / wat1,relh1,therm,dt,dtemp
      common / aqrate / sodium, pres
      common / count/ icount, iprint
      common/aqchem2/watcont,crustal,salt,watvap
c
c     INITIALIZATION (ONLY ON FIRST CALL)
c     (CALCULATION OF SECTIONAL DIAMETERS, LOADING OF MOLEC. WEIGHTS)
c
      if (iaq .eq. 1) then
      call dropinit
      call aqdist(fdist,fdist2)
      iaq = 0
      endif
      pres = p                                   ! pressure in atm
      call constants(temp)
c
c     ZERO THE YAQ MATRIX
c
      do i=1, meqn2
      yaq(i)=0.0
      enddo
c      
c     TRANSFER VIA COMMON VARIABLES TO AQRATE
c
      watcont(1) = lwc(1)
      watcont(2) = lwc(2)
      relh1 = rh
      therm = temp
      dt = deltat
      icount=0
      iprint=20
c
c     TRANSFER OF ALL GAS-PHASE CONCENTRATIONS TO GASCON AND THEN
c     THROUGH COMMON BLOCK TO RATES
c
      do i=1, ngas_aq
      gascon(i) = gas(i)
      enddo
c
c     UNIT CONVERSION FACTORS     
c      
      fso2=1000.*p*wso2/(8.314e-2*temp)     !SO2   conver. factor from ppm to ug/m3
      fh2o2=1000.*p*wh2o2/(8.314e-2*temp)   !H2O2  conver. factor from ppm to ug/m3
      fhcho=1000.*p*whcho/(8.314e-2*temp)   !HCHO  conver. factor from ppm to ug/m3
      fhcooh=1000.*p*whcooh/(8.314e-2*temp) !HCOOH conver. factor from ppm to ug/m3
      fammon=1000.*p*wnh3/(8.314e-2*temp)   !NH3   conver. factor from ppm to ug/m3
c    
c      
c   *** BEGINNING OF AQUEOUS-PHASE CALCULATION
c   LOADING OF THE CONCENTRATIONS IN THE YAQ VECTOR
c
c   THE TOTAL FORMALDEHYDE/FORMIC ACID IS TRANSFERRED AS GAS-PHASE SPECIES
c    (IT IS NOT INCLUDED IN THE AQUEOUS OR AEROSOL MATRICES)
      yaq(1) = gas(nghcho)*fhcho           ! total HCHO (ug/m3)
c
c   GAS-PHASE SPECIES
c
      yaq(2) = gas(ngso2)*fso2             ! total SO2 (ug/m3)
      yaq(3) = gas(ngh2o2)*fh2o2           ! total H2O2 (ug/m3)                    
      yaq(4) = gas(nga) *fammon            ! NH3(g) in ug/m3
      yaq(5) = gas(nghcooh)*fhcooh         ! HCOOH (total) (ug/m3)
c
c   *** CALCULATION FOR FIRST AQUEOUS SECTION ***
c
c
c A FIXED SECTIONAL DIAMETER IS USED FOR THE ACTIVATION
c CALCULATION. THIS IS COMPATIBLE WITH THE ASSUMPTION THAT
c A CLOUD OR FOG EXISTS IF THE RH EXCEEDS A THRESHOLD VALUE
c (FOR EXAMPLE RH>85%).
c
c     ** THE DROPLETS ARE ASSUMED TO BY WITHOUT SO2 OR
c        H2O2 IN THE BEGINNING OF A TIMESTEP **
c      (THIS IS DONE TO AVOID CARRYING THE S(IV)/H2O2 CONCENTRATIONS
c      IN THE 3D SIMULATION **

      yaq(6) = 1.e-10     ! SO2 (aq) in ug/m3
      yaq(7) = 1.e-10     ! H2O2 (aq) in ug/m3
c
      do isect=1, nsect
      if (daer(isect) .ge. dactiv .AND. daer(isect) .lt. dsep) then
      yaq(8) = yaq(8) + aerosol(isect,nan)       ! NITRATE (aq) in ug/m3 
      yaq(9) = yaq(9) + aerosol(isect,nac)       ! CHLORIDE (aq) in ug/m3
      yaq(10) = yaq(10) + aerosol(isect,naa)     ! AMMONIUM (aq) in ug/m3   
      yaq(11) = yaq(11) + aerosol(isect,na4)     ! SULFATE (aq) in ug/m
c
      yaq(12) = yaq(12) + aerosol(isect,nahso5)       ! HSO5- in ug/m3
      yaq(13) = yaq(13) + aerosol(isect,nahmsa)       ! HMSA in ug/m3
      endif
      enddo

c   *** CALCULATION FOR SECOND AQUEOUS SECTION ***
c
      yaq(14) = 1.e-10     ! SO2 (aq) in ug/m3
      yaq(15) = 1.e-10     ! H2O2 (aq) in ug/m3
c
      do isect=1, nsect
      if (daer(isect) .ge. dsep) then
      yaq(16) = yaq(16) + aerosol(isect,nan)       ! NITRATE (aq) in ug/m3 
      yaq(17) = yaq(17) + aerosol(isect,nac)       ! CHLORIDE (aq) in ug/m3
      yaq(18) = yaq(18) + aerosol(isect,naa)     ! AMMONIUM (aq) in ug/m3   
      yaq(19) = yaq(19) + aerosol(isect,na4)     ! SULFATE (aq) in ug/m3
      yaq(20) = yaq(20) + aerosol(isect,nahso5)       ! HSO5- in ug/m3
      yaq(21) = yaq(21) + aerosol(isect,nahmsa)       ! HMSA in ug/m3      
      endif
      enddo
c
c     CALCULATION OF THE DISSOLUTION OF THE AVAILABLE HNO3 AND HCl
c     TO THE DROPLETS/PARTICLES. WE ASSUME THAT THE TIMESCALE OF DISSOLUTION
c     IS SMALL AND THAT ALL HNO3 AND HCl ARE DISSOLVED (ZERO VAPOR PRESSURE)
c
      fhno3=1000.*whno3*pres/(8.314e-2*temp) !HNO3 conver. factor from ppm to ug/m3
      fhcl=1000.*whcl*pres/(8.314e-2*temp)   !HCl  conver. factor from ppm to ug/m3
c
      hno3 = gas(ngcn)*fhno3   ! HNO3(g) (in ug/m3)
      hcl =  gas(ngcc)*fhcl    ! HCl(g) (in ug/m3)
c      
c     SAVE THE OLD AQUEOUS-PHASE CONCENTRATIONS BEFORE THE INTEGRATION
c
      tnitold1= yaq(8)
      chlorold1= yaq(9)
      ammonold1= yaq(10)
      sulfateold1= yaq(11)
      hso5old1 = yaq(12)
      hmsaold1 = yaq(13)
      tnitold2= yaq(16)
      chlorold2= yaq(17)
      ammonold2= yaq(18)
      sulfateold2= yaq(19)
      hso5old2 = yaq(20)
      hmsaold2 = yaq(21)      
c
c      ** ADD THE GAS-PHASE CONCENTRATIONS TO THE TOTAL FOR HCl AND HNO3 **
c
c      * REMOVE IT FROM THE GAS-PHASE *
      gas(ngcn) = 0.0          ! all HNO3 is dissolved
      gas(ngcc) = 0.0          ! all HCl is dissolved
c      * CALCULATE PARTITIONING FACTORS *
      z1 = lwc(1) / diamet1
      z2 = lwc(2) / diamet2
      z = z1+z2
c
      fr1 = z1/z               ! fraction that goes to section 1
      fr2 = z2/z               ! fraction that goes to section 2
c
      yaq(8) = yaq(8)+hno3*fr1     ! HNO3 increase Section 1
      yaq(9) = yaq(9)+hcl*fr1      ! HCl increase Section 1
      yaq(16) = yaq(16)+hno3*fr2   ! HNO3 increase Section 2
      yaq(17) = yaq(17)+hcl*fr2    ! HCl increase Section 2
c
c
c.....THE YAQ MATRIX HAS BEEN INITIALIZED
c     CALCULATION OF SODIUM VECTOR (NEEDED FOR THE INTEGRATION ROUTINE)
c       (TRANSFERRED VIA COMMON STATEMENT AQRATES1)
c
      do isect=1, nsect
      sodium(isect) = aerosol(isect,nas)
      enddo
c
c     VARIABLES FOR INTEGRATION
c
      stime = tbeg*60.                        ! in seconds
      stout = (tbeg+deltat) * 60.              ! in seconds
c
c
c     CALCULATE CRUSTAL SPECIES CONCENTRATION INSIDE DROPLETS
c     (IT IS USED IN THE AQUEOUS-PHASE CHEMISTRY CALCULATION TO
c      ESTIMATE Fe and Mn CONCENTRATIONS
c
      crustal(1) = aerosol(5,kcru)+aerosol(6,kcru)+aerosol(7,kcru)
      crustal(2) = aerosol(8,kcru)+aerosol(9,kcru)+aerosol(10,kcru)
      salt(1) = aerosol(5,ksod)+aerosol(6,ksod)+aerosol(7,ksod)
      salt(2) = aerosol(8,ksod)+aerosol(9,ksod)+aerosol(10,ksod)
c
c     INTEGRATE
c 
      call aqintegr2(meqn2, yaq, stime, stout, temp, p)
c
c     CALCULATE THE DIFFERENCES
c
      dnit1 = yaq(8) - tnitold1
      dchlor1 = yaq(9) - chlorold1
      dammon1 = yaq(10) - ammonold1
      dsulf1 = yaq(11) - sulfateold1
      dhso51 = yaq(12) - hso5old1
      dhmsa1 = yaq(13) - hmsaold1
c
      dnit2 = yaq(16) - tnitold2
      dchlor2 = yaq(17) - chlorold2
      dammon2 = yaq(18) - ammonold2
      dsulf2 = yaq(19) - sulfateold2
      dhso52 = yaq(20) - hso5old2
      dhmsa2 = yaq(21) - hmsaold2      
c
c     ADD THE MASS CHANGE TO THE CORRESPONDING SECTION 
c
c      do isect=5, 7  ! for sections 0.1 to 10 um
      do isect=4, 6  ! for sections 0.04 to 40 um
      aerosol(isect,nan)=aerosol(isect,nan)+dnit1*fdist2(isect)    ! NITRATE (aq) in ug/m3
      aerosol(isect,nac)=aerosol(isect,nac)+dchlor1*fdist2(isect)  ! CHLORIDE (aq) in ug/m3
      aerosol(isect,naa)=aerosol(isect,naa)+dammon1*fdist2(isect)  ! AMMONIUM (aq) in ug/m3   
      aerosol(isect,na4)=aerosol(isect,na4)+dsulf1*fdist2(isect)   ! SULFATE (aq) in ug/m3
      aerosol(isect,nahso5)=aerosol(isect,nahso5)+dhso51*fdist2(isect)      
      aerosol(isect,nahmsa)=aerosol(isect,nahmsa)+dhmsa1*fdist2(isect)
      enddo
c
c      do isect=8, 10 ! for sections 0.1 to 10 um
      do isect=7, 8  ! for sections 0.04 to 40 um
      aerosol(isect,nan)=aerosol(isect,nan)+dnit2*fdist2(isect)    ! NITRATE (aq) in ug/m3
      aerosol(isect,nac)=aerosol(isect,nac)+dchlor2*fdist2(isect)  ! CHLORIDE (aq) in ug/m3
      aerosol(isect,naa)=aerosol(isect,naa)+dammon2*fdist2(isect)  ! AMMONIUM (aq) in ug/m3   
      aerosol(isect,na4)=aerosol(isect,na4)+dsulf2*fdist2(isect)   ! SULFATE (aq) in ug/m3
      aerosol(isect,nahso5)=aerosol(isect,nahso5)+dhso52*fdist2(isect)      
      aerosol(isect,nahmsa)=aerosol(isect,nahmsa)+dhmsa2*fdist2(isect)      
      enddo

c     CHECK IF THE ALGORITHM RESULTED IN NEGATIVE CONCENTRATIONS
c     REPORT THE CORRECTIONS AT THE WARNING FILE  
c
c  check(1,isp) = total aerosol concentration for species isp     
c  check(2,isp) = concentration below zero for species isp
c  check(3,isp) = (-) total mass added to the system to get a positive conc.

      do j = 1,3
      do isp = 1,naers
      check(j,isp) = 0.0
      enddo
      enddo
      
      do isp = 1,naers
      do isect=1,nsect
      check(1,isp) = check(1,isp) + aerosol(isect,isp)
      enddo
      enddo
      
c If negative in a section, change to a small positive value
c Account for this mass change by subtracting from other 
c non-negative sections.      
      
      do isp = 1,naers
      do isect = 1,nsect
      if (aerosol(isect,isp) .lt. 0.0) then
      check(2,isp) = check(2,isp)+aerosol(isect,isp)
cbk      check(3,isp) = check(3,isp)+check(2,isp)-1.e-12
      check(3,isp) = check(3,isp)+aerosol(isect,isp)-1.e-12 ! bkoo (10/07/03)
      aerosol(isect,isp) = 1.e-12
      endif
      enddo
      enddo
      
      do isect = 1,nsect
      
      if (aerosol(isect,nan) .gt. 1.e-12) then
      aerosol(isect,nan)=aerosol(isect,nan)+(aerosol(isect,nan))*
     & check(3,nan)/(check(1,nan)-check(2,nan))
      endif
      if (aerosol(isect,nac) .gt. 1.e-12) then
      aerosol(isect,nac)=aerosol(isect,nac)+(aerosol(isect,nac))*
     & check(3,nac)/(check(1,nac)-check(2,nac))
      endif
      if (aerosol(isect,naa) .gt. 1.e-12) then
      aerosol(isect,naa)=aerosol(isect,naa)+(aerosol(isect,naa))*
     & check(3,naa)/(check(1,naa)-check(2,naa))
      endif                  

c For HSO5 and HMSA, if negative and there is not enough mass in the
c remaining aerosol sections to account for the added mass, subtract
c from sulfate.
     
      if (check(1,nahso5)-check(2,nahso5) .lt. 0.-check(3,nahso5)) then
      check(3,na4) = check(3,na4)+check(3,nahso5)*96./113.
      else if (aerosol(isect,nahso5) .gt. 1.e-12) then
      aerosol(isect,nahso5)=aerosol(isect,nahso5)+
     & (aerosol(isect,nahso5))*check(3,nahso5)/
     & (check(1,nahso5)-check(2,nahso5))
      endif
      
      if (check(1,nahmsa)-check(2,nahmsa) .lt. 0.-check(3,nahmsa)) then
      check(3,na4) = check(3,na4)+check(3,nahmsa)*96./111.
      else if (aerosol(isect,nahmsa) .gt. 1.e-12) then
      aerosol(isect,nahmsa)=aerosol(isect,nahmsa)+
     & (aerosol(isect,nahmsa))*check(3,nahmsa)/
     & (check(1,nahmsa)-check(2,nahmsa))
      endif  
                  
      if (aerosol(isect,na4) .gt. 1.e-12) then
      aerosol(isect,na4)=aerosol(isect,na4)+(aerosol(isect,na4))*
     & check(3,na4)/(check(1,na4)-check(2,na4))
      endif
     
      enddo 

 

c
c     RETURN REST OF YAQ RESULTS TO THEIR ORIGINAL MATRICES (GAS,AEROSOL)
c
c     ** GAS-PHASE SPECIES **     
c     IMPORTANT : THE DISSOLVED S(IV) and H2O2 ARE TRANSFERRED
c                 BACK TO THE GAS PHASE
c     
      gas(nghcho) = yaq(1)/fhcho              ! total HCHO (ppm)
      gas(ngso2)= (yaq(2)+yaq(6)*0.7901+yaq(14)*0.7901)/fso2  ! SO2(g) (ppm)
      gas(ngh2o2)=(yaq(3)+yaq(7)+yaq(15))/fh2o2             ! H2O2(g) (ppm)                    
      gas(nga)=yaq(4)/fammon                  ! NH3(g) (ppm)
      gas(nghcooh)= yaq(5)/fhcooh             ! total HCOOH (ppm)
c
c THE CODE NEGLECTS THE CHANGE IN THE SIZE DISTRIBUTION SHAPE
c OF THE NONVOLATILE AEROSOL COMPONENTS BECAUSE OF THE SULFATE
c PRODUCTION. THE CHANGE IN THE SULFATE DISTRIBUTION IS CALCULATED
c USING THE DISTRIBUTION FUNCTIONS WHILE THE CHANGE IN THE
c DISTRIBUTION OF THE VOLATILE INORGANIC AEROSOL COMPONENTS WILL
c CALCULATED BY THE AEROSOL MODULE AFTER CLOUD EVAPORATION
c
c
c.....END OF CODE
      return
      end
