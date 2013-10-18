      subroutine aqoperator1(tbeg,deltat,gas,aerosol,lwc,
     &  rh,temp,p,iaq)
c
c
c
c.....INPUTS
c     TBEG : BEGINNING OF THE INTEGRATION (in min)
c     DELTAT : INTEGRATION TIME (in min)
c     GAS(NGAS) : GAS-PHASE CONCENTRATION VECTOR (in ppm)
c     AEROSOL(NSECT,NAERS) : AEROSOL SPECIES CONCENTRATION (in ug/m3) 
c     LWC : LIQUID WATER CONTENT (in g/m3)
c     RH : RELATIVE HUMIDITY  (in 0-1 scale)
c     TEMP : TEMPERATURE (in K)
c     P : PRESSURE (in atm)
c     IAQ : FLAG: 1 AT FIRST CALL, 0 THERE AFTER
c
c.....OUTPUTS
c     GAS(NGAS) : UPDATED GAS-PHASE CONCENTRATIONS (in ppm)
c     AEROSOL(NSECT,NAERS) : UPDATED AEROSOL SPECIES CONCENTRATION (in ug/m3) 
c     IAQ : FLAG SET TO ZERO
c
c.....VARIABLES
c     GCON(NGAS) : GAS-PHASE CONCENTRATIONS (in ppm) (THROUGH COMMON TO RATES)
c
c
c.....DIFFERENTIAL EQUATIONS ARE SOLVED FOR THE FOLLOWING SPECIES      
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
c.....Y MATRIX STRUCTURE  EVERYTHING IS (ug/m3)
c     ONLY 11 ODEs ARE SOLVED FOR THE WHOLE DISTRIBUTION
c
c  YAQ(1) = TOTAL FORMALDEHYDE
c  YAQ(2) = TOTAL FORMIC ACID
c  YAQ(3) = SO2 (g)
c  YAQ(4) = H2O2 (g)
c  YAQ(5) = NH3 (g)
c  YAQ(6) = SO2 (aq)
c  YAQ(7) = H2O2 (aq)
c  YAQ(8) = NITRATE (aq)
c  YAQ(9) = CHLORIDE (aq)          
c  YAQ(10) = AMMONIUM (aq)
c  YAQ(11) = SULFATE (aq)
c  YAQ(12) = HSO5- (aq)
c  YAQ(13) = HMSA (aq)
c     
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
      real*4 gas(ngas_aq), aerosol(nsect,naers)
      real*4 sodium(nsect)
      real*4 deltat, lwc, rh, temp
      real*4 gascon(ngas_aq)
      real*4 yaq(meqn1)
      real*4 watcont, watvap
      real*4 check(3,naers)
c
      common / gascon / gascon
      common / epcomy / ymin, hmax
      common / integration / wat1,relh1,therm,dt,dtemp
      common / aqrate / sodium, pres
      common / count/ icount, iprint
      common/aqchem1/watcont,watvap,crustal,salt
c
c     INITIALIZATION (ONLY ON FIRST CALL)
c     (CALCULATION OF SECTIONAL DIAMETERS, LOADING OF MOLEC. WEIGHTS)

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
      do i=1, meqn1
      yaq(i)=0.0
      enddo
c      
c     TRANSFER VIA COMMON VARIABLES TO AQRATE
c
      wat1 = lwc
      watcont = lwc
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
      fso2=1000.*pres*wso2/(8.314e-2*temp)     !SO2   conver. factor from ppm to ug/m3
      fh2o2=1000.*pres*wh2o2/(8.314e-2*temp)   !H2O2  conver. factor from ppm to ug/m3
      fhcho=1000.*pres*whcho/(8.314e-2*temp)   !HCHO  conver. factor from ppm to ug/m3
      fhcooh=1000.*pres*whcooh/(8.314e-2*temp) !HCOOH conver. factor from ppm to ug/m3
      fammon=1000.*pres*wnh3/(8.314e-2*temp)   !NH3   conver. factor from ppm to ug/m3
c    
c      
c   *** BEGINNING OF AQUEOUS-PHASE CALCULATION
c   LOADING OF THE CONCENTRATIONS IN THE YAQ VECTOR
c
c   THE TOTAL FORMALDEHYDE/FORMIC ACID ARE TRANSFERRED AS GAS-PHASE SPECIES
c    (THEY ARE NOT INCLUDED IN THE AQUEOUS OR AEROSOL MATRICES)
      yaq(1) = gas(nghcho)*fhcho           ! total HCHO (ug/m3)
      yaq(2) = gas(nghcooh)*fhcooh         ! total HCOOH (ug/m3)
c
c   GAS-PHASE SPECIES
c
      yaq(3) = gas(ngso2)*fso2             ! total SO2 (ug/m3)
      yaq(4) = gas(ngh2o2)*fh2o2           ! total H2O2 (ug/m3)                    
      yaq(5) = gas(nga) *fammon            ! NH3(g) in ug/m3
c
c   AQUEOUS-PHASE SPECIES
c
      do isect=1, nsect
c
c A FIXED SECTIONAL DIAMETER IS USED FOR THE ACTIVATION
c CALCULATION. THIS IS COMPATIBLE WITH THE ASSUMPTION THAT
c A CLOUD OR FOG EXISTS IF THE RH EXCEEDS A THRESHOLD VALUE
c (FOR EXAMPLE RH>85%).
c
c
c     ** THE DROPLETS ARE ASSUMED TO BY WITHOUT SO2 OR
c        H2O2 IN THE BEGINNING OF A TIMESTEP **
c      (THIS IS DONE TO AVOID CARRYING THE S(IV)/H2O2 CONCENTRATIONS
c      IN THE 3D SIMULATION **
      yaq(6) = 1.e-10     ! SO2 (aq) in ug/m3
      yaq(7) = 1.e-10     ! H2O2 (aq) in ug/m3
c
      if (daer(isect) .ge. dactiv) then
      yaq(8) = yaq(8) + aerosol(isect,nan)       ! NITRATE (aq) in ug/m3 
      yaq(9) = yaq(9) + aerosol(isect,nac)       ! CHLORIDE (aq) in ug/m3
      yaq(10) = yaq(10) + aerosol(isect,naa)     ! AMMONIUM (aq) in ug/m3   
      yaq(11) = yaq(11) + aerosol(isect,na4)     ! SULFATE (aq) in ug/m3
      yaq(12) = yaq(12) + aerosol(isect,nahso5)  ! HSO5- in ug/m3
      yaq(13) = yaq(13) + aerosol(isect,nahmsa)  ! HMSA in ug/m3
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
      tnitold=yaq(8)
      chlorold=yaq(9)
      ammonold=yaq(10)
      sulfateold=yaq(11)
      hso5old=yaq(12)
      hmsaold=yaq(13)
c
c      ** ADD THE GAS-PHASE CONCENTRATIONS TO THE TOTAL FOR HCl AND HNO3 **
      gas(ngcn) = 0.0          ! all HNO3 is dissolved
      gas(ngcc) = 0.0          ! all HCl is dissolved
      yaq(8) = yaq(8)+hno3     ! HNO3 increase
      yaq(9) = yaq(9)+hcl      ! HCl increase
c
c
c.....THE YAQ MATRIX HAS BEEN INITIALIZED
c
c     CALCULATION OF SODIUM VECTOR (NEEDED FOR THE INTEGRATION ROUTINE)
c       (TRANSFERRED VIA COMMON STATEMENT AQRATESA)
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
c     CALCULATE CRUSTAL SPECIES CONCENTRATION INSIDE DROPLETS
c     (IT IS USED IN THE AQUEOUS-PHASE CHEMISTRY CALCULATION TO
c      ESTIMATE Fe and Mn CONCENTRATIONS
c
      crustal = 0.
      salt =0.
      do isect=1,nsect
c
      if (dd(isect) .ge. dactiv) then
      crustal=crustal+aerosol(isect,nar)
      salt=salt+aerosol(isect,nas)
      endif
c
      enddo   
c
c     INTEGRATE
c 
      sbef=0.
      do i=1,nsect
       sbef=sbef+aerosol(i,na4)*32./96.+aerosol(i,nahso5)*32./113.+
     &           aerosol(i,nahmsa)*32./111.
      enddo
      sbef=sbef+yaq(3)*32./64.

      call aqintegr1(meqn1, yaq, stime, stout, temp, p)

c
c     CALCULATE THE DIFFERENCES
c
      dnit = yaq(8) - tnitold
      dchlor = yaq(9) - chlorold
      dammon = yaq(10) - ammonold
      dsulf = yaq(11) - sulfateold
      dhso5 = yaq(12) - hso5old
      dhmsa = yaq(13) - hmsaold
c
c     ADD THE MASS CHANGE TO THE CORRESPONDING SECTION 
c
      do isect=1, nsect
      aerosol(isect,nan)=aerosol(isect,nan)+dnit*fdist(isect)    ! NITRATE (aq) in ug/m3
      aerosol(isect,nac)=aerosol(isect,nac)+dchlor*fdist(isect)  ! CHLORIDE (aq) in ug/m3
      aerosol(isect,naa)=aerosol(isect,naa)+dammon*fdist(isect)  ! AMMONIUM (aq) in ug/m3   
      aerosol(isect,na4)=aerosol(isect,na4)+dsulf*fdist(isect)   ! SULFATE (aq) in ug/m3
      aerosol(isect,nahso5)=aerosol(isect,nahso5)
     & +dhso5*fdist(isect)   					! hso5 (aq) in ug/m3
      aerosol(isect,nahmsa)=aerosol(isect,nahmsa)
     & +dhmsa*fdist(isect)   					! hmsa (aq) in ug/m3
      enddo
      
c
      saft=0.
      do i=1,nsect
       saft=saft+aerosol(i,na4)*32./96.+aerosol(i,nahso5)*32./113.+
     &           aerosol(i,nahmsa)*32./111.
      enddo
      saft=saft+(yaq(3)+0.7901*yaq(6))*32./64.
c     CHECK IF THE ALGORITHM RESULTED IN NEGATIVE CONCENTRATIONS
c     REPORT THE CORRECTIONS AT THE WARNING FILE  
c
c      do isect=1, nsect
c      do isp=1, naers     
c      if (aerosol(isect,isp) .lt. 0.0) then
c      write(2,*)'NEGATIVE CONCENTRATION AT', stime
c      write(2,*)'SECTION=',isect, 'SPECIES=',isp, aerosol(isect,isp)
c     
c      write(6,*)'NEGATIVE CONCENTRATION AT', stime	!added by TMR
c      write(6,*)'SECTION=',isect, ' SPECIES=',isp, aerosol(isect,isp)
c
cckf      aerosol(isect,isp)=1.e-12 
c      endif 
c      enddo 
c      enddo 

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

      saft1=0.
      do i=1,nsect
       saft1=saft1+aerosol(i,na4)*32./96.+aerosol(i,nahso5)*32./113.+
     &           aerosol(i,nahmsa)*32./111.
      enddo
      saft1=saft1+(yaq(3)+0.7901*yaq(6))*32./64.

c
c     RETURN YAQ RESULTS TO THEIR ORIGINAL MATRICES (GAS,AEROSOL)
c
c     ** GAS-PHASE SPECIES **     
c     IMPORTANT : THE DISSOLVED S(IV) and H2O2 ARE TRANSFERRED
c                 BACK TO THE GAS PHASE
c     
      gas(nghcho) = yaq(1)/fhcho              ! total HCHO (ppm)
      gas(nghcooh)= yaq(2)/fhcooh             ! total HCOOH (ppm)
      gas(ngso2)= (yaq(3)+yaq(6)*0.7901)/fso2       ! SO2(g) (ppm)
      gas(ngh2o2)=(yaq(4)+yaq(7))/fh2o2             ! H2O2(g) (ppm)                    
      gas(nga)=yaq(5)/fammon               ! NH3(g) (ppm)
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
