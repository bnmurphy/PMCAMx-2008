      subroutine RAQCHEM ( TEMP, PRES_PA, TAUCLD, 
     &                     WCAVG, GAS, AEROSOL,
     &                     idiag, iout, igrd, iaq, jaq, kaq)
camx
c  This is a stand-alone version of the aqueous phase chemistry from
c  CMAQ version 4 (November 2002). Deposition calculations are
c  commented out. Aerosol concentrations by size mode have been replaced
c  by total concentrations.
camx
c      SUBROUTINE AQCHEM ( JDATE, JTIME, TEMP, PRES_PA, TAUCLD, PRCRATE,
c     &                    WCAVG, WTAVG, AIRM, ALFA0, ALFA2, ALFA3, GAS,
c     &                    AEROSOL, GASWDEP, AERWDEP, HPWDEP )
C-----------------------------------------------------------------------
C
C  DESCRIPTION:
C    Compute concentration changes in cloud due to aqueous chemistry,
C    scavenging and wet deposition amounts.
C
C  Revision History:
C      No   Date   Who	What
C      -- -------- ---  -----------------------------------------
C      0  / /86    CW   BEGIN PROGRAM - Walceks's Original Code
C      1  / /86    RB   INCORPORATE INTO RADM
C      2  03/23/87 DH   REFORMAT
C      3  04/11/88 SJR  STREAMLINED CODE - ADDED COMMENTS
C      4  08/27/88 SJR  COMMENTS, MODIFIED FOR RPM
C      4a 03/15/96 FSB  Scanned hard copy to develop Models3
C                       Version.
C      5  04/24/96 FSB  Made into Models3 Format
C      6  02/18/97 SJR  Revisions to link with Models3
C      7  08/12/97 SJR  Revised for new concentration units (moles/mole)
C                       and new treatment of nitrate and nitric acid
C      8  01/15/98 sjr  revised to add new aitken mode scavenging
C                       and aerosol number scavenging
C      9  12/15/98 David Wong at LM:
C             -- change division of XL, TEMP to multiplication of XL, TEMP
C                reciprocal, respectively
C             -- change / TOTOX / TSIV to / ( TOTOX * TSIV )
C     10  03/18/99 David Wong at LM:
C             -- removed "* 1.0" redundant calculation at TEMP1 calculation
C     11  04/27/00 sjr  Added aerosol surface area as modeled species
C
C  Reference:
C     Walcek & Taylor, 1986, A theoretical Method for computing
C      vertical distributions of acidity and sulfate within cumulus
C      clouds, J. Atmos Sci.,  Vol. 43, no. 4 pp 339 - 355
C
C  Called by:  AQMAP
C
C  Calls the following subroutines:  none
C
C  Calls the following functions:  HLCONST
C
C  ARGUMENTS     TYPE      I/O       DESCRIPTION
C  ---------     ----  ------------  --------------------------------
C  GAS(ngas)     real  input&output  Concentration for species i=1,11
C  GASWDEP(ngas) real     output     wet deposition for species
C                                    (1) = SO2   conc (mol/mol of S02)
C                                    (2) = HNO3  conc (mol/mol of HNO3)
C                                    (3) = N2O5  conc (mol/mol of N2O5)
C                                    (4) = CO2   conc (mol/mol of CO2)
C                                    (5) = NH3   conc (mol/mol of NH3)
C                                    (6) = H2O2  conc (mol/mol of H2O2)
C                                    (7) = O3    conc (mol/mol of O3)
C                                    (8) = FOA   conc (mol/mol of FOA)
C                                    (9) = MHP   conc (mol/mol of MHP)
C                                    (10)= PAA   conc (mol/mol of PAA)
C                                    (11)= H2SO4 conc (mol/mol of H2SO4)
C
camx
c
c  AEROSOL concentrations for species i=1,9 in camx version
c                                    (1) = SO4    conc (mol/mol)
c                                    (2) = NH4    conc (mol/mol)
c                                    (3) = NO3    conc (mol/mol)
c                                    (4) = CAC    conc (mol/mol)
c                                    (5) = MGC    conc (mol/mol)
c                                    (6) = NACL   conc (mol/mol)
c                                    (7) = A3FE   conc (mol/mol)
c                                    (8) = B2MN   conc (mol/mol)
c                                    (9) = KCL    conc (mol/mol)
c
camx
C  AEROSOL(naer) real input&output   Concentration for species i=1,21
C  AERWDEP(naer) real     output     wet deposition for species
C                                    (1) = SO4AKN conc (mol/mol)
C                                    (2) = SO4ACC conc (mol/mol)
C                                    (3) = NH4AKN conc (mol/mol)
C                                    (4) = NH4ACC conc (mol/mol)
C                                    (5) = NO3AKN conc (mol/mol)
C                                    (6) = NO3ACC conc (mol/mol)
C                                    (7) = NO3COR conc (mol/mol)
C                                    (8) = ORGAKN conc (mol/mol)
C                                    (9) = ORGACC conc (mol/mol)
C                                    (10)= PRIAKN conc (mol/mol)
C                                    (11)= PRIACC conc (mol/mol)
C                                    (12)= PRICOR conc (mol/mol)
C                                    (13)= CACO3  conc (mol/mol)
C                                    (14)= MGCO3  conc (mol/mol)
C                                    (15)= NACL   conc (mol/mol)
C                                    (16)= A3FE   conc (mol/mol)
C                                    (17)= B2MN   conc (mol/mol)
C                                    (18)= KCL    conc (mol/mol)
C                                    (19)= NUMAKN conc (#/mol)
C                                    (20)= NUMACC conc (#/mol)
C                                    (21)= NUMCOR conc (#/mol)
C                                    (22)= SRFAKN conc (m2/mol)
C                                    (23)= SRFACC conc (m2/mol)
C
C-----------------------------------------------------------------------

      IMPLICIT NONE

c      INCLUDE SUBST_CONST          ! constants
c      INCLUDE SUBST_XSTAT          ! M3EXIT status codes
c      INCLUDE 'AQ_PARAMS.EXT'      ! aqueous chemistry shared parameters
camx
c Essential parameters from the include files
c
       integer     NGAS     ! local number of gasses for aqchem
       integer     NAER     ! local number of aerosols for aqchem
       REAL        MOLVOL   ! Molar volume at STP [ L/mol ] Non MKS units 
       REAL        STDATMPA ! standard atmosphere  [ Pa ]
       REAL        STDTEMP  ! Standard Temperature [ K ]
       parameter ( NGAS = 11 )
       parameter ( NAER = 9  )
       PARAMETER ( MOLVOL = 22.41410 )
       PARAMETER ( STDATMPA = 101325.0 )
       PARAMETER ( STDTEMP = 273.15 )
camx
c      CHARACTER*120 XMSG           ! Exit status message
c      DATA          XMSG / ' ' /

C...........PARAMETERS and their descriptions:

      INTEGER      NUMOX           ! number of oxidizing reactions
      PARAMETER  ( NUMOX =  5 )

      REAL         H2ODENS         ! density of water at 20 C and 1 ATM
      PARAMETER  ( H2ODENS = 1000.0 )  ! (kg/m3)

      INTEGER      NLIQS           ! number of liquid phase species
      PARAMETER  ( NLIQS = 33 )

      REAL         ONETHIRD       ! 1/3
      PARAMETER  ( ONETHIRD = 1.0 / 3.0 )

      REAL         TWOTHIRDS       ! 2/3
      PARAMETER  ( TWOTHIRDS = 2.0 / 3.0 )
camx
c........... Gas species pointers
      integer      lso2
      integer      lhno3
      integer      ln2o5
      integer      lco2
      integer      lnh3
      integer      lh2o2
      integer      lo3
      integer      lfoa
      integer      lmhp
      integer      lpaa
      integer      lh2so4
      parameter    (lso2   = 1 )
      parameter    (lhno3  = 2 )
      parameter    (ln2o5  = 3 )
      parameter    (lco2   = 4 )
      parameter    (lnh3   = 5 )
      parameter    (lh2o2  = 6 )
      parameter    (lo3    = 7 )
      parameter    (lfoa   = 8 )
      parameter    (lmhp   = 9 )
      parameter    (lpaa   = 10)
      parameter    (lh2so4 = 11)
c........... Aerosol species pointers
      integer      lso4
      integer      lnh4
      integer      lno3
      integer      lcaco3
      integer      lmgco3
      integer      lnacl
      integer      la3fe
      integer      lb2mn
      integer      lkcl
      parameter    (lso4   = 1 )
      parameter    (lnh4   = 2 )
      parameter    (lno3   = 3 )
      parameter    (lcaco3 = 4 )
      parameter    (lmgco3 = 5 )
      parameter    (lnacl  = 6 )
      parameter    (la3fe  = 7 )
      parameter    (lb2mn  = 8 )
      parameter    (lkcl   = 9 )
camx
C...........ARGUMENTS and their descriptions

c      INTEGER      JDATE           ! current model date, coded YYYYDDD
c      INTEGER      JTIME           ! current model time, coded HHMMSS

c      REAL         AIRM            ! total air mass in cloudy layers (mol/m2)
c      REAL         ALFA0           ! scav coef for aitken aerosol number
c      REAL         ALFA2           ! scav coef for aitken aerosol sfc area
c      REAL         ALFA3           ! scav coef for aitken aerosol mass
c      REAL         HPWDEP          ! hydrogen wet deposition (mm mol/liter)
c      REAL         PRCRATE         ! precip rate (mm/hr)
      REAL         PRES_PA         ! pressure (Pa)
      REAL         TAUCLD          ! timestep for cloud (s)
      REAL         TEMP            ! temperature (K)
      REAL         WCAVG           ! liquid water content (kg/m3)
c      REAL         WTAVG           ! total water content (kg/m3)
      REAL         GAS    ( NGAS ) ! gas phase concentrations (mol/molV)
      REAL         AEROSOL( NAER ) ! aerosol concentrations (mol/molV)
c      REAL         GASWDEP( NGAS ) ! gas phase wet deposition array (mm mol/liter)
c      REAL         AERWDEP( NAER ) ! aerosol wet deposition array (mm mol/liter)
      integer      idiag           ! unit number for diagnostic output
      integer      iout            ! unit number for message output
      integer      igrd            ! grid number
      integer      iaq             ! i grid cell index (column)
      integer      jaq             ! j grid cell index (row)
      integer      kaq             ! k grid cell index (layer)

C...........LOCAL VARIABLES (scalars) and their descriptions:

      CHARACTER*16 PNAME          ! driver program name
      SAVE         PNAME
      DATA         PNAME / 'AQCHEM' /

      INTEGER      I20C            ! loop counter for do loop 20
      INTEGER      I30C            ! loop counter for do loop 30
      INTEGER      ITERAT          ! # iterations of aqueaous chemistry solver
      INTEGER      I7777C          ! aqueous chem iteration counter
      INTEGER      ICNTAQ          ! aqueous chem iteration counter
c      INTEGER      LIQ             ! loop counter for liquid species
      INTEGER      IOX             ! index over oxidation reactions

c      REAL         DEPSUM
c      REAL         BETASO4
      REAL         A               ! iron's anion concentration
      REAL         AC              ! H+ concentration in cloudwater (mol/liter)
      REAL         ACT1            ! activity corretion factor!single ions
      REAL         ACT2            ! activity factor correction!double ions
      REAL         ACTB            !
      REAL         AE              ! guess for H+ conc in cloudwater (mol/liter)
      REAL         B               ! manganese's anion concentration
      REAL         PRES_ATM        ! pressure (Atm)
      REAL         BB              ! lower limit guess of cloudwater pH
      REAL         CA              ! Calcium conc in cloudwater (mol/liter)
      REAL         CAA             ! inital Calcium in cloudwater (mol/liter)
c      REAL         NO3CORA         ! initial NO3COR in cloudwater (mol/liter)
      REAL         CL              ! Cl-  conc in cloudwater (mol/liter)
      REAL         CLA             ! initial Cl in cloudwater (mol/liter)
      REAL         CO2H            ! Henry's Law constant for CO2
      REAL         CO21            ! First dissociation constant for CO2
      REAL         CO22            ! Second dissociation constant for CO2
      REAL         CO212           ! CO21*CO22
      REAL         CO212H          ! CO2H*CO21*CO22
      REAL         CO21H           ! CO2H*CO21
      REAL         CO2L            ! CO2 conc in cloudwater (mol/liter)
      REAL         CO3             ! CO3= conc in cloudwater (mol/liter)
      REAL         CO3A            ! initial CO3 in cloudwater (mol/liter)
c      REAL         CTHK1           ! cloud thickness (m)
      REAL         DTRMV           !
      REAL         DTS6            !
      REAL         FA              ! functional value ??
      REAL         FB              ! functional value ??
      REAL         FE              ! Fe+++ conc in cloudwater (mol/liter)
      REAL         FEA             ! initial Fe in cloudwater (mol/liter)
      REAL         FNH3            ! frac weight of NH3 to total ammonia
c      REAL         FNH4ACC         ! frac weight of NH4 acc to total ammonia
      REAL         FHNO3           ! frac weight of HNO3 to total NO3
c      REAL         FNO3ACC         ! frac weight of NO3 acc to total NO3
c      REAL         FNO3COR         ! frac weight of NO3 cor to total NO3
c      REAL         FRACACC         ! frac ACC that was from accum mode
c      REAL         FRACCOR         ! frac NO3 that was from coarse mode
      REAL         FOA1            ! First dissociation constant for FOA
      REAL         FOAH            ! Henry's Law constant for FOA
      REAL         FOA1H           ! FOAH*FOA1
      REAL         FOAL            ! FOA conc in cloudwater (mol/liter)
      REAL         FTST            !
      REAL         GM              !
      REAL         GM1             !
      REAL         GM1LOG          !
      REAL         GM2             ! activity correction factor
      REAL         GM2LOG          !
      REAL         HA              !
      REAL         HB              !
      REAL         H2OW            !
      REAL         H2O2H           ! Henry's Law Constant for H2O2
      REAL         H2O2L           ! H2O2 conc in cloudwater (mol/liter)
      REAL         HCO2            ! HCO2 conc in cloudwater (mol/liter)
      REAL         HCO3            ! HCO3 conc in cloudwater (mol/liter)
      REAL         HNO3H           ! Henry's Law Constant for HNO3
      REAL         HNO31           ! First dissociation constant for HNO3
      REAL         HNO31H          !
      REAL         HNO3L           ! HNO3 conc in cloudwater (mol/liter)
      REAL         HSO3            ! HSO3 conc in cloudwater (mol/liter)
      REAL         HSO4            ! HSO4 concn in cloudwater (mol/liter)
      REAL         HTST            !
      REAL         K               ! K conc in cloudwater (mol/liter)
      REAL         KA              ! initial K in cloudwater (mol/liter)
      REAL         LGTEMP          ! log of TEMP
c      REAL         M3NEW           ! accumulation mode mass at time t
c      REAL         M3OLD           ! accumulation mode mass at time 0
      REAL         MG              !
      REAL         MGA             ! inital Mg in cloudwater (mol/liter)
      REAL         MHPH            ! Henry's Law Constant for MHP
      REAL         MHPL            ! MHP conc in cloudwater (mol/liter)
      REAL         MN              ! Mn++ conc in cloudwater (mol/liter)
      REAL         MNA             ! initial Mn in cloudwater (mol/liter)
      REAL         NA              ! Na conc in cloudwater (mol/liter)
      REAL         NAA             ! initial Na in cloudwater (mol/liter)
      REAL         NH31            ! First dissociation constant for NH3
      REAL         NH3H            ! Henry's Law Constant for NH3
      REAL         NH3DH20         !
      REAL         NH31HDH         !
      REAL         NH3L            ! NH3 conc in cloudwater (mol/liter)
      REAL         NH4             ! NH4+ conc in cloudwater (mol/liter)
c      REAL         NH4AKNA         ! init NH4 akn conc in cloudwater (mol/liter)
      REAL         NH4ACCA         ! init NH4 acc conc in cloudwater (mol/liter)
c      REAL         NITAER          ! total aerosol nitrate 
      REAL         NO3             ! NO3 conc in cloudwater (mol/liter)
      REAL         NO3ACCA         ! init NO3 acc conc in cloudwater (mol/liter)
c      REAL         NO3AKNA         ! init NO3 akn conc in cloudwater (mol/liter)
      REAL         O3H             ! Henry's Law Constant for O3
      REAL         O3L             ! O3 conc in cloudwater (mol/liter)
      REAL         OH              ! OH conc in cloudwater (mol/liter)
c      REAL         ORGN            ! ORGANIC aerosol in cloudwater (mol/liter)
c      REAL         ORGACCA         ! init ORG ACC aerosol in cloudwater (mol/liter)
c      REAL         ORGAKNA         ! init ORG AKN aerosol in cloudwater (mol/liter)
      REAL         PAAH            ! Henry's Law Constant for PAA
      REAL         PAAL            ! PAA conc in cloudwater (mol/liter)
      REAL         PCO20           ! total CO2 partial pressure (atm)
      REAL         PCO2F           ! gas only CO2 partial pressure (atm)
      REAL         PFOA0           ! total ORGANIC acid partial pressure (atm)
      REAL         PFOAF           ! gas only ORGANIC ACID partial press (atm)
      REAL         PH2O20          ! total H2O2 partial pressure (atm)
      REAL         PH2O2F          ! gas only H2O2 partial pressure (atm)
      REAL         PHNO30          ! total HNO3 partial pressure (atm)
      REAL         PHNO3F          ! gas only HNO3 partial pressure (atm)
      REAL         PMHP0           ! total MHP partial pressure (atm)
      REAL         PMHPF           ! gas only MHP partial pressure (atm)
      REAL         PNH30           ! total NH3 partial pressure (atm)
      REAL         PNH3F           ! gas only NH3 partial pressure (atm)
      REAL         PO30            ! total O3 partial pressure (atm)
      REAL         PO3F            ! gas only O3 partial pressure (atm)
      REAL         PPAA0           ! total PAA partial pressure (atm)
      REAL         PPAAF           ! gas only PAA partial pressure (atm)
c      REAL         PRIM            ! PRIMARY acc+akn aerosol in cloudwater (mol/liter)
c      REAL         PRIMCOR         ! PRIMARY coarse aerosol in cloudwater (mol/liter)
c      REAL         PRIACCA         ! init PRI ACC aerosol in cloudwater (mol/liter)
c      REAL         PRIAKNA         ! init PRI AKN aerosol in cloudwater (mol/liter)
c      REAL         PRICORA         ! init PRI COR aerosol in cloudwater (mol/liter)
      REAL         PSO20           ! total SO2 partial pressure (atm)
      REAL         PSO2F           ! gas only SO2 partial pressure (atm)
      REAL         RATE            !
      REAL         RECIPA1         !
      REAL         RECIPA2         !
      REAL         RECIPAP1        ! one over pressure (/atm)
      REAL         RH2O2           !
      REAL         RMHP            !
      REAL         RPAA            !
      REAL         RT              ! gas const * temperature (liter atm/mol)
      REAL         SCVEFF          ! Scavenging efficiency (%)
      SAVE         SCVEFF
      DATA         SCVEFF / 100.0 / ! currently set to 100%
      REAL         SIV             ! dissolved so2 in cloudwater (mol/liter)
      REAL         SK6             !
      REAL         SK6TS6          !
      REAL         SO21            ! First dissociation constant for SO2
      REAL         SO22            ! Second dissociation constant for SO2
      REAL         SO2H            ! Henry's Law Constant for SO2
      REAL         SO212           ! SO21*SO22
      REAL         SO212H          ! SO21*SO22*SO2H
      REAL         SO21H           ! SO21*SO2H
      REAL         SO2L            ! SO2 conc in cloudwater (mol/liter)
      REAL         SO3             ! SO3= conc in cloudwater (mol/liter)
      REAL         SO4             ! SO4= conc in cloudwater (mol/liter)
      REAL         STION           ! ionic strength
      REAL         TAC             !
      REAL         TEMP1           !
      REAL         TIMEW           ! cloud chemistry clock (sec)
      REAL         TOTOX           !
      REAL         TOTAMM          ! total ammonium
      REAL         TOTNIT          ! total nitrate
      REAL         TS6             ! SO4 conc in cloudwater (mol/liter)
c      REAL         TS6AKNA         ! init SO4 akn conc in cloudwater (mol/liter)
      REAL         TS6ACCA         ! init SO4 acc conc in cloudwater (mol/liter)
      REAL         TSIV            !
      REAL         TST             !
c      REAL         XC1             ! (/mm)
c      REAL         XC2             ! (liter-atm/mol/mm)
      REAL         XL              ! conversion factor (liter-atm/mol)
      REAL         ONE_OVER_XL     ! 1.0 / XL
      REAL         PRES_ATM_OVER_XL     ! PRES_ATM / XL
      REAL         XLCO2           !
      REAL         XLH2O2          !
      REAL         XLHNO3          !
      REAL         XLMHP           !
      REAL         XLNH3           !
      REAL         XLO3            !
      REAL         XLPAA           !
      REAL         XLSO2           !

C...........LOCAL VARIABLES (arrays) and their descriptions:

      REAL         LIQUID( NLIQS ) ! liquid concentration array (mol/liter)
c      REAL         WETDEP( NLIQS ) ! wet deposition array (mm mol/liter)
      REAL         DSIVDT( 0:NUMOX ) ! rate of so2 oxid incloud (mol/liter/sec)
      REAL         DS4   ( 0:NUMOX ) ! S(IV) oxidized over timestep DTW(0)
      REAL         DTW   ( 0:NUMOX ) ! cloud chemistry timestep (sec)

      REAL         ONE_OVER_TEMP     ! 1.0 / TEMP

C...........EXTERNAL FUNCTIONS and their descriptions:

      REAL          HLCONST
      EXTERNAL      HLCONST

C*********************************************************************
C     begin body of subroutine AQCHEM

      ONE_OVER_TEMP = 1.0 / TEMP

C...check for bad temperature, cloud air mass, or pressure

c      IF ( TEMP .LE. 0.0 ) THEN
c        IF ( AIRM .LE. 0.0 ) THEN
c          IF ( PRES_PA .LE. 0.0 ) THEN
      IF ( TEMP .LE. 0.0 .OR. PRES_PA .LE. 0.0 ) THEN
c            XMSG = 'MET DATA ERROR'
c            CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         write(iout,'(//,A)') 'ERROR in RAQCHEM:'
         write(iout,'(/,A)') 
     &            'Error in RADM aqueous chemistry (RAQCHEM)'
         write(iout,'(A,3E10.3)') 
     &            'Invalid met data: T, P',TEMP,PRES_PA
c     &            'Invalid met data: T, AIRM, P',TEMP,AIRM,PRES_PA
         write(iout,'(a,4i4)') 
     &         'Grid, icell, jcell, kcell:', igrd, iaq, jaq, kaq
         call camxerr()
c          END IF
c        END IF
      END IF

C...compute several conversion factors

      ICNTAQ = 0
      ITERAT = 0
      RT = ( MOLVOL / STDTEMP ) * TEMP             ! R * T (liter atm / mol)
      PRES_ATM = PRES_PA /  STDATMPA               ! pressure (atm)
c      CTHK1 = AIRM * RT / ( PRES_ATM * 1000.0 )    ! cloud thickness (m)
      XL   = WCAVG * RT / H2ODENS     ! conversion factor (l-atm/mol)
      ONE_OVER_XL = 1.0 / XL
      PRES_ATM_OVER_XL = PRES_ATM / XL
      TST  = 0.999
      GM   = SCVEFF / 100.0
      ACT1 = 1.0
      ACT2 = 1.0
      GM2  = 1.0
      TIMEW = 0.0
      RECIPAP1 = 1.0 / PRES_ATM
c      XC1  = 1.0 / ( WCAVG * CTHK1 )
c      XC2  = RT / ( 1000.0 * CTHK1 )

C...set equilibrium constants as a function of temperature
C...   Henry's law constants

      SO2H  = HLCONST( 'SO2             ', TEMP )
      CO2H  = HLCONST( 'CO2             ', TEMP )
      NH3H  = HLCONST( 'NH3             ', TEMP )
      H2O2H = HLCONST( 'H2O2            ', TEMP )
      O3H   = HLCONST( 'O3              ', TEMP )
      HNO3H = HLCONST( 'HNO3            ', TEMP )
      MHPH  = HLCONST( 'METHYLHYDROPEROX', TEMP )
      PAAH  = HLCONST( 'PEROXYACETIC_ACI', TEMP )
      FOAH  = HLCONST( 'FORMIC_ACID     ', TEMP )

C...dissociation constants

      FOA1 = 1.71E-4

C...From ERT book

      SK6 = 10.0**( 1180.0 * ONE_OVER_TEMP - 5.95 )

C...From Maahs <1982>

      SO21 = 10.0**( 853.0 * ONE_OVER_TEMP ) / 5.495E4
      SO22 = 10.0**( 621.9 * ONE_OVER_TEMP ) / 1.897E9

C...From Edwards et al. <1978>

      LGTEMP = ALOG( TEMP )

      CO21 = 
     &  10.0**( -5251.5 * ONE_OVER_TEMP - 15.94   * LGTEMP + 102.269 )
      CO22 =
     &   10.0**( -5401.4 * ONE_OVER_TEMP - 15.4096 * LGTEMP + 95.574 )
      H2OW =
     &   10.0**( -5839.5 * ONE_OVER_TEMP - 9.7618  * LGTEMP + 61.206 )

C...From Morgan and Maas <1931>

      NH31 = 10.0**( -189.1 * ONE_OVER_TEMP - 4.117 )

C...Schwartz and White <1981, 1983>

      HNO31 = 15.4

C...Kinetic oxidation rates
C...   From Chamedies (1982)

      TEMP1 = ONE_OVER_TEMP - 1.0 / 298.0

      RH2O2 = 80000.0 * EXP( -3650.0 * TEMP1 )

C...From Kok

      RMHP = 1.75E7 * EXP( -3801.0 * TEMP1 )
      RPAA = 3.64E7 * EXP( -3994.0 * TEMP1 )

C...make initializations

c      DO LIQ = 1, NLIQS
c        WETDEP( LIQ ) = 0.0
c      END DO

      DO IOX = 0, NUMOX
        DSIVDT( IOX ) = 0.0
        DTW   ( IOX ) = 0.0
        DS4   ( IOX ) = 0.0
      END DO

C...compute the initial accumulation aerosol 3rd moment
c
c      M3OLD = ( AEROSOL( LSO4ACC ) * SGRAERMW( LSO4ACC ) / 1.8e6
c     &      +   AEROSOL( LNH4ACC ) * SGRAERMW( LNH4ACC ) / 1.8e6
c     &      +   AEROSOL( LNO3ACC ) * SGRAERMW( LNO3ACC ) / 1.8e6
c     &      +   AEROSOL( LORGACC ) * SGRAERMW( LORGACC ) / 2.0e6
c     &      +   AEROSOL( LPRIACC ) * SGRAERMW( LPRIACC ) / 2.2e6 )
ccc     &      * 6.0 / PI    ! cancels out in division at end of subroutine
c
C...compute fractional weights for several species
camx
      FHNO3   = 1.0
      TOTNIT = GAS( LHNO3 ) + AEROSOL( LNO3 )
      IF ( TOTNIT .GT. 0.0 ) FHNO3   = GAS( LHNO3 ) / TOTNIT
c
      FNH3    = 1.0
      TOTAMM = GAS( LNH3 ) + AEROSOL( LNH4 )
      IF ( TOTAMM .GT. 0.0 ) FNH3    = GAS( LNH3 ) / TOTAMM
camx
c
c      NITAER = AEROSOL( LNO3ACC ) + AEROSOL( LNO3COR )
c      IF ( NITAER .GT. 0.0 ) THEN
c        FRACACC = AEROSOL( LNO3ACC ) / NITAER
c        FRACCOR = AEROSOL( LNO3COR ) / NITAER
c      ELSE
c        FRACACC = 1.0
c        FRACCOR = 0.0
c      END IF
c      
c      TOTNIT = GAS( LHNO3 ) + AEROSOL( LNO3ACC ) + AEROSOL( LNO3COR )
c      IF ( TOTNIT .GT. 0.0 ) THEN
c        FHNO3   = GAS( LHNO3 ) / TOTNIT
c        FNO3ACC = AEROSOL( LNO3ACC ) / TOTNIT
c        FNO3COR = AEROSOL( LNO3COR ) / TOTNIT
c      ELSE
c        FHNO3   = 1.0
c        FNO3ACC = 0.0
c        FNO3COR = 0.0
c      END IF
c
c      TOTAMM = GAS( LNH3 ) + AEROSOL( LNH4ACC )
c      IF ( TOTAMM .GT. 0.0 ) THEN
c        FNH3    = GAS( LNH3 ) / TOTAMM
c        FNH4ACC = AEROSOL( LNH4ACC ) / TOTAMM
c      ELSE
c        FNH3    = 1.0
c        FNH4ACC = 0.0
c      END IF
c
C...initial concentration from accumulation-mode aerosol loading (mol/liter)
C...  an assumption is made that all of the accumulation-mode
C...  aerosol mass in incorporated into the cloud droplets
camx
      TS6ACCA = ( AEROSOL( LSO4 )
     &        +   GAS    ( LH2SO4  ) ) * PRES_ATM_OVER_XL
      NO3ACCA =   AEROSOL( LNO3 )   * PRES_ATM_OVER_XL
      NH4ACCA =   AEROSOL( LNH4 )   * PRES_ATM_OVER_XL
camx
c      TS6ACCA = ( AEROSOL( LSO4ACC )
c     &        +   GAS    ( LH2SO4  ) ) * PRES_ATM_OVER_XL
c      NO3ACCA =   AEROSOL( LNO3ACC )   * PRES_ATM_OVER_XL
c      NH4ACCA =   AEROSOL( LNH4ACC )   * PRES_ATM_OVER_XL
c      ORGACCA =   AEROSOL( LORGACC )   * PRES_ATM_OVER_XL
c      PRIACCA =   AEROSOL( LPRIACC )   * PRES_ATM_OVER_XL

C...initial concentration from coarse-mode aerosol loading (mol/liter)
C...  an assumption is made that all of the coarse-mode
C...  aerosol mass in incorporated into the cloud droplets

      CLA     = ( AEROSOL( LNACL   )
     &         * PRES_ATM_OVER_XL )
c     &        +   AEROSOL( LKCL    ) ) * PRES_ATM_OVER_XL
c      NO3CORA =   AEROSOL( LNO3COR )   * PRES_ATM_OVER_XL
      CAA     =   AEROSOL( LCACO3  )   * PRES_ATM_OVER_XL
      MGA     =   AEROSOL( LMGCO3  )   * PRES_ATM_OVER_XL
      NAA     =   AEROSOL( LNACL   )   * PRES_ATM_OVER_XL
      KA      =   AEROSOL( LKCL    )   * PRES_ATM_OVER_XL
      FEA     =   AEROSOL( LA3FE   )   * PRES_ATM_OVER_XL
      MNA     =   AEROSOL( LB2MN   )   * PRES_ATM_OVER_XL
      CO3A    = ( AEROSOL( LCACO3  )
     &        +   AEROSOL( LMGCO3  ) ) * PRES_ATM_OVER_XL
c      PRICORA =   AEROSOL( LPRICOR )   * PRES_ATM_OVER_XL

C...set constant factors that will be used in later multiplications (moles/atm)

      XLH2O2  = H2O2H * XL
      XLO3    = O3H   * XL
      XLMHP   = MHPH  * XL
      XLPAA   = PAAH  * XL
      XLSO2   = SO2H  * XL
      XLNH3   = NH3H  * XL
      XLHNO3  = HNO3H * XL
      XLCO2   = CO2H  * XL

      SO212   = SO21  * SO22
      SO21H   = SO21  * SO2H
      SO212H  = SO212 * SO2H
      CO212   = CO21  * CO22
      CO21H   = CO21  * CO2H
      CO212H  = CO22  * CO21H
      NH3DH20 = NH31  / H2OW
      NH31HDH = NH3H  * NH3DH20
      FOA1H   = FOA1  * FOAH
      HNO31H  = HNO31 * HNO3H

C...If kinetic calculations are made, return to this point

      I20C = 0
20    CONTINUE

      I20C = I20C + 1
      IF ( I20C .GE. 1000 ) THEN
c        XMSG = 'EXCESSIVE LOOPING AT I20C'
c        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         write(iout,'(//,A)') 'ERROR in RAQCHEM:'
         write(iout,'(/,A)') 'Error in RADM aqueous chemistry (AQCHEM)'
         write(iout,'(A)') 'EXCESSIVE LOOPING AT I20C'
         write(iout,'(a,4i4)') 
     &         'Grid, icell, jcell, kcell:', igrd, iaq, jaq, kaq
         call camxerr()
      END IF

C...set aitken-mode aerosol loading (mol/liter)
c
c      NO3AKNA = AEROSOL( LNO3AKN ) * PRES_ATM_OVER_XL
c     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
c      NH4AKNA = AEROSOL( LNH4AKN ) * PRES_ATM_OVER_XL
c     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
c      TS6AKNA = AEROSOL( LSO4AKN ) * PRES_ATM_OVER_XL
c     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
c      ORGAKNA = AEROSOL( LORGAKN ) * PRES_ATM_OVER_XL
c     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
c      PRIAKNA = AEROSOL( LPRIAKN ) * PRES_ATM_OVER_XL
c     &        * ( 1.0 - EXP( -ALFA3 * TIMEW ) )
c
C...Initial gas phase partial pressures (atm)
C...   = initial partial pressure - amount deposited partial pressure

      PSO20  = GAS( LSO2  ) * PRES_ATM
     &       + DS4( 0 ) * XL
c     &       - ( WETDEP(  8 ) + WETDEP(  9 ) + WETDEP( 10 ) ) * XC2
      PNH30  = GAS( LNH3  ) * PRES_ATM
     &       + ( NH4ACCA ) * XL
c     &       + ( NH4ACCA + NH4AKNA ) * XL
c     &       - ( WETDEP(  2 ) + WETDEP( 15 ) ) * XC2
      PHNO30 = ( GAS( LHNO3 ) + 2.0 * GAS( LN2O5 ) ) * PRES_ATM
     &       + ( NO3ACCA ) * XL
c      PHNO30 = ( GAS( LHNO3 ) + 2.0 * GAS( LN2O5 ) ) * PRES_ATM
c     &       + ( NO3ACCA + NO3CORA + NO3AKNA ) * XL
c     &       - ( WETDEP( 14 ) + WETDEP( 32 ) ) * XC2
      PH2O20 = GAS( LH2O2 ) * PRES_ATM 
      PO30   = GAS( LO3   ) * PRES_ATM 
c      PH2O20 = GAS( LH2O2 ) * PRES_ATM - WETDEP( 17 ) * XC2
c      PO30   = GAS( LO3   ) * PRES_ATM - WETDEP( 18 ) * XC2
      PFOA0  = GAS( LFOA  ) * PRES_ATM
c     &       - ( WETDEP( 22 ) + WETDEP( 23 ) ) * XC2
      PMHP0  = GAS( LMHP  ) * PRES_ATM
      PPAA0  = GAS( LPAA  ) * PRES_ATM
c      PMHP0  = GAS( LMHP  ) * PRES_ATM - WETDEP( 24 ) * XC2
c      PPAA0  = GAS( LPAA  ) * PRES_ATM - WETDEP( 25 ) * XC2
      PCO20  = GAS( LCO2  ) * PRES_ATM
     &       + CO3A * XL
c     &       - ( WETDEP( 11 ) + WETDEP( 12 ) + WETDEP( 13 ) ) * XC2

C...don't allow gas concentrations to go below zero

      PSO20  = MAX( PSO20,  0.0 )
      PNH30  = MAX( PNH30,  0.0 )
      PH2O20 = MAX( PH2O20, 0.0 )
      PO30   = MAX( PO30,   0.0 )
      PFOA0  = MAX( PFOA0,  0.0 )
      PMHP0  = MAX( PMHP0,  0.0 )
      PPAA0  = MAX( PPAA0,  0.0 )
      PCO20  = MAX( PCO20,  0.0 )
      PHNO30 = MAX( PHNO30, 0.0 )

C...Molar concentrations of soluble aerosols
C...   = Initial amount - amount deposited  (mol/liter)

      TS6     = TS6ACCA - DS4( 0 )
      CL      = CLA
      CA      = CAA
      MG      = MGA
      NA      = NAA
      K       = KA
      FE      = FEA
      MN      = MNA
      A       = 3.0 * FE
      B       = 2.0 * MN
c      TS6     = TS6ACCA  + TS6AKNA
c     &        - ( WETDEP(  6 ) + WETDEP(  7 ) ) * XC1
c     &        - DS4( 0 )
c      CL      = CLA      -   WETDEP( 16 )  * XC1
c      CA      = CAA      -   WETDEP(  3 )  * XC1
c      MG      = MGA      -   WETDEP( 29 )  * XC1
c      NA      = NAA      -   WETDEP(  4 )  * XC1
c      K       = KA       -   WETDEP( 30 )  * XC1
c      FE      = FEA      -   WETDEP( 19 )  * XC1
c      MN      = MNA      -   WETDEP( 20 )  * XC1
c      ORGN    = ORGACCA + ORGAKNA - WETDEP( 27 )  * XC1
c      PRIM    = PRIACCA + PRIAKNA - WETDEP( 28 )  * XC1
c      PRIMCOR = PRICORA  -   WETDEP( 33 )  * XC1
c      A       = 3.0 * FE
c      B       = 2.0 * MN
c
C...don't allow aerosol concentrations to go below zero

      TS6     = MAX( TS6,     0.0 )
      CL      = MAX( CL,      0.0 )
      CA      = MAX( CA,      0.0 )
      MG      = MAX( MG,      0.0 )
      NA      = MAX( NA,      0.0 )
c      K       = MAX( K,       0.0 )
      FE      = MAX( FE,      0.0 )
      MN      = MAX( MN,      0.0 )
c      ORGN    = MAX( ORGN,    0.0 )
c      PRIM    = MAX( PRIM,    0.0 )
c      PRIMCOR = MAX( PRIMCOR, 0.0 )
      A       = MAX( A,       0.0 )
      B       = MAX( B,       0.0 )

      SK6TS6 = SK6 * TS6

C...find solution of the equation using a method of reiterative
C...  bisections Make initial guesses for pH:   between .01  to  10.

      HA =  0.01
      HB = 10.0

      I7777C = 0
7777  CONTINUE

      I7777C = I7777C + 1
      IF ( I7777C .GE. 1000 ) THEN
c        XMSG = 'EXCESSIVE LOOPING AT I7777C'
c        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         write(iout,'(//,A)') 'ERROR in RAQCHEM:'
         write(iout,'(/,A)') 'Error in RADM aqueous chemistry (AQCHEM)'
         write(iout,'(A)') 'EXCESSIVE LOOPING AT I7777C'
         write(iout,'(a,4i4)') 
     &         'Grid, icell, jcell, kcell:', igrd, iaq, jaq, kaq
         call camxerr()
      END IF

      HA = MAX( HA - 0.8, 0.1 )
      HB = MIN( HB + 0.8, 9.9 )
      AE = 10.0**( -HA )

      RECIPA1 = 1.0 / ( AE * ACT1 )
      RECIPA2 = 1.0 / ( AE * AE * ACT2 )

C...calculate final gas phase partial pressure of SO2, NH3, HNO3
C...  HCOOH, and CO2 (atm)

      PSO2F = PSO20 / ( 1.0 + XLSO2 * ( 1.0 + SO21 * RECIPA1
     &      + SO212 * RECIPA2 ) )

      PNH3F = PNH30 / ( 1.0 + XLNH3 * ( 1.0 + NH3DH20 * AE ) )

      PFOAF = PFOA0 / ( 1.0 + XL * ( FOAH + FOA1H * RECIPA1 ) )

      PHNO3F = PHNO30 / ( 1.0 + XLHNO3 * ( 1.0 + HNO31 * RECIPA1 ) )

      PCO2F = PCO20 / ( 1.0 + XLCO2 * ( 1.0 + CO21 * RECIPA1
     &      + CO212 * RECIPA2 ) )

C...calculate liquid phase concentrations (moles/liter)

      SO4  = SK6TS6 / ( AE * GM2 + SK6 )
      HSO4 = TS6 - SO4
      SO3  = SO212H  * PSO2F  * RECIPA2
      HSO3 = SO21H   * PSO2F  * RECIPA1
      CO3  = CO212H  * PCO2F  * RECIPA2
      HCO3 = CO21H   * PCO2F  * RECIPA1
      OH   = H2OW    * RECIPA1
      NH4  = NH31HDH * PNH3F  * AE
      HCO2 = FOA1H   * PFOAF  * RECIPA1
      NO3  = HNO31H  * PHNO3F * RECIPA1

C...compute functional value

      FA = AE + NH4 + 2.0 *  (CA + MG - CO3 - SO3 - SO4 ) - OH - HCO3
     &   - HSO3 - NO3 - HSO4 - HCO2

C...Start iteration and bisection ****************<<<<<<<

      I30C = 0
30    CONTINUE

      I30C = I30C + 1
      IF ( I30C .GE. 1000 ) THEN
c        XMSG = 'EXCESSIVE LOOPING AT I30C'
c        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         write(iout,'(//,A)') 'ERROR in RAQCHEM:'
         write(iout,'(/,A)') 'Error in RADM aqueous chemistry (AQCHEM)'
         write(iout,'(A)') 'EXCESSIVE LOOPING AT I30C'
         write(iout,'(a,4i4)') 
     &         'Grid, icell, jcell, kcell:', igrd, iaq, jaq, kaq
         call camxerr()
      END IF

      BB = ( HA + HB ) / 2.0
      AE = 10.0**( -BB )

      ICNTAQ = ICNTAQ + 1
      IF ( ICNTAQ .GE. 3000 ) THEN
c        XMSG = 'Maximum AQCHEM total iterations exceeded'
c        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         write(iout,'(//,A)') 'ERROR in RAQCHEM:'
         write(iout,'(/,A)') 'Error in RADM aqueous chemistry (AQCHEM)'
         write(iout,'(A)') 'EXCESSIVE LOOPING AT I3000C'
         write(iout,'(a,4i4)') 
     &         'Grid, icell, jcell, kcell:', igrd, iaq, jaq, kaq
         call camxerr()
      END IF

      RECIPA1 = 1.0 / ( AE * ACT1 )
      RECIPA2 = 1.0 / ( AE * AE * ACT2 )

C...calculate final gas phase partial pressure of SO2, NH3, HNO3
C...  HCOOH, and CO2 (atm)

      PSO2F = PSO20 / ( 1.0 + XLSO2
     &	    * ( 1.0 + SO21 * RECIPA1 + SO212 * RECIPA2 ) )

      PNH3F = PNH30 / ( 1.0 + XLNH3 * ( 1.0 + NH3DH20 * AE ) )

      PHNO3F = PHNO30 / ( 1.0 + XLHNO3 * ( 1.0 + HNO31 * RECIPA1 ) )

      PFOAF = PFOA0 / ( 1.0 + XL * ( FOAH + FOA1H * RECIPA1 ) )

      PCO2F = PCO20 / ( 1.0 + XLCO2 * ( 1.0 + CO21 * RECIPA1
     &      + CO212 * RECIPA2 ) )

C...calculate liquid phase concentrations (moles/liter)

      SO4  = SK6TS6 / ( AE * GM2 + SK6 )
      HSO4 = TS6 - SO4
      SO3  = SO212H  * PSO2F  * RECIPA2
      HSO3 = SO21H   * PSO2F  * RECIPA1
      CO3  = CO212H  * PCO2F  * RECIPA2
      HCO3 = CO21H   * PCO2F  * RECIPA1
      OH   = H2OW    * RECIPA1
      NH4  = NH31HDH * PNH3F  * AE
      HCO2 = FOA1H   * PFOAF  * RECIPA1
      NO3  = HNO31H  * PHNO3F * RECIPA1

C...compute functional value

      FB = AE + NH4 + 2.0 * ( CA + MG - CO3 - SO3 - SO4 ) - OH - HCO3
     &   - HSO3 - NO3 - HSO4 - HCO2

C...Calculate and check the sign of the product of the two functional values

      FTST = FA * FB
      IF ( FTST .LE. 0.0 ) THEN
        HB = BB
      ELSE
        HA = BB
        FA = FB
      END IF

C...Check convergence of solutions

      HTST = HA / HB
      IF ( HTST .LE. TST ) GO TO 30

C...end of zero-finding routine ****************<<<<<<<<<<<<

C...compute Ionic strength and activity coefficient by the Davies equation

      STION = 0.5 * (AE + NH4 + OH + HCO3 + HSO3
     &      + 4.0 * (SO4 + CO3 + SO3 + CA + MG + MN)
     &      + NO3 + HSO4 + 9.0 * FE + NA + CL + A + B + HCO2)
c     &      + NO3 + HSO4 + 9.0 * FE + NA + K + CL + A + B + HCO2)
      GM1LOG = -0.509 * ( SQRT( STION )
     &       / ( 1.0 + SQRT( STION ) ) - 0.2 * STION )
      GM2LOG = GM1LOG * 4.0
      GM1  = 10.0**GM1LOG
      GM2  = MAX( 10.0**GM2LOG, 1.0E-30 )
      ACTB = ACT1
      ACT1 = MAX( GM1 * GM1, 1.0E-30 )
      ACT2 = MAX( GM1 * GM1 * GM2, 1.0E-30 )

C...check for convergence and possibly go to 7777, to recompute
C...  Gas and liquid phase concentrations

      TAC = ABS( ACTB - ACT1 ) / ACTB
      IF ( TAC .GE. 1.0E-2 ) GO TO 7777

C...return an error if the pH is not in range

ccc      IF ( ( HA .LT. 0.02 ) .OR. ( HA .GT. 9.49 ) ) THEN
      IF ( ( HA .LT. 0.1 ) .OR. ( HA .GT. 9.9 ) ) THEN
c        print *, ha
c        XMSG = 'PH VALUE OUT OF RANGE'
c        CALL M3EXIT ( PNAME, JDATE, JTIME, XMSG, XSTAT2 )
         write(iout,'(//,A)') 'ERROR in RAQCHEM:'
         write(iout,'(/,A)') 'Error in RADM aqueous chemistry (AQCHEM)'
         write(iout,'(A,f5.1)') 'PH VALUE OUT OF RANGE: ', ha
         write(iout,'(a,4i4)') 
     &         'Grid, icell, jcell, kcell:', igrd, iaq, jaq, kaq
         call camxerr()
      END IF

C...Make those concentration calculations which can be made outside
C...  of the function.

      SO2L = SO2H * PSO2F
      AC = 10.0**( -BB )
      SIV = SO3 + HSO3 + SO2L

C...Calculate final gas phase concentrations of oxidants (atm)

      PH2O2F = ( PH2O20 + XL * DS4( 1 ) ) / ( 1.0 + XLH2O2 )
      PO3F   = ( PO30   + XL * DS4( 2 ) ) / ( 1.0 + XLO3   )
      PMHPF  = ( PMHP0  + XL * DS4( 4 ) ) / ( 1.0 + XLMHP  )
      PPAAF  = ( PPAA0  + XL * DS4( 5 ) ) / ( 1.0 + XLPAA  )

      PH2O2F = MAX( PH2O2F, 0.0 )
      PO3F   = MAX( PO3F,   0.0 )
      PMHPF  = MAX( PMHPF,  0.0 )
      PPAAF  = MAX( PPAAF,  0.0 )

C...Calculate liquid phase concentrations of oxidants (moles/liter)

      H2O2L = PH2O2F * H2O2H
      O3L   = PO3F   * O3H
      MHPL  = PMHPF  * MHPH
      PAAL  = PPAAF  * PAAH
      FOAL  = PFOAF  * FOAH
      NH3L  = PNH3F  * NH3H
      CO2L  = PCO2F  * CO2H
      HNO3L = PHNO3F * HNO3H

C...load the liquid concentration array with current values

      LIQUID(  1 ) = AC
      LIQUID(  2 ) = NH4
      LIQUID(  3 ) = CA
      LIQUID(  4 ) = NA
      LIQUID(  5 ) = OH
      LIQUID(  6 ) = SO4
      LIQUID(  7 ) = HSO4
      LIQUID(  8 ) = SO3
      LIQUID(  9 ) = HSO3
      LIQUID( 10 ) = SO2L
      LIQUID( 11 ) = CO3
      LIQUID( 12 ) = HCO3
      LIQUID( 13 ) = CO2L
      LIQUID( 14 ) = NO3
      LIQUID( 15 ) = NH3L
      LIQUID( 16 ) = CL
      LIQUID( 17 ) = H2O2L
      LIQUID( 18 ) = O3L
      LIQUID( 19 ) = FE
      LIQUID( 20 ) = MN
      LIQUID( 21 ) = A
      LIQUID( 22 ) = FOAL
      LIQUID( 23 ) = HCO2
      LIQUID( 24 ) = MHPL
      LIQUID( 25 ) = PAAL
      LIQUID( 26 ) = 0.0
c      LIQUID( 27 ) = ORGN
c      LIQUID( 28 ) = PRIM
      LIQUID( 29 ) = MG
      LIQUID( 30 ) = K
      LIQUID( 31 ) = B
      LIQUID( 32 ) = HNO3L
c      LIQUID( 33 ) = PRIMCOR

C...if the maximum cloud lifetime has not been reached, the compute
C...  the next timestep.

      IF ( TIMEW .LT. TAUCLD ) THEN

C...make kinetics calculations
C...  note: DS4(i) and DSIV(I) are negative numbers!

        DTRMV = 300.0
c        IF ( ( CTHK1 .GT. 1.0E-10 ) .AND. ( PRCRATE .GT. 1.0E-10 ) )
c     &     DTRMV = 3.6 * WTAVG * 1000.0 * CTHK1 / PRCRATE  ! <<<uma found bug, was .36
c        DTRMV = MIN( DTRMV, 300.0 )
        ITERAT = ITERAT + 1

C...Define the total S(iv) available for oxidation

        TSIV = PSO20 * ONE_OVER_XL

C...Calculate sulfur iv oxidation rate due to H2O2

        DSIVDT( 1 ) = -RH2O2 * H2O2L * SO2L / ( 0.1 + AC )
        IF ( ( DSIVDT( 1 ) .EQ. 0.0 ) .OR. ( TSIV .LE. 1.0E-30 ) ) THEN
          DTW(1) = DTRMV
        ELSE
          TOTOX = PH2O20 * ONE_OVER_XL
          RATE = -DSIVDT( 1 ) / ( TOTOX * TSIV )
          DTW( 1 ) = 0.05 / ( RATE * MAX( TOTOX, TSIV ) )
        END IF

C...Calculate sulfur iv oxidation rate due to O3

        IF ( BB .GE. 2.7 ) THEN
          DSIVDT( 2 ) = -4.19E5 * ( 1.0 + 2.39E-4 / AC ) * O3L * SIV
        ELSE
          DSIVDT( 2 ) = -1.9E4 * SIV * O3L / SQRT( AC )
        END IF
        IF ( ( DSIVDT( 2 ) .EQ. 0.0 ) .OR. ( TSIV .LE. 1.0E-30 ) ) THEN
          DTW( 2 ) = DTRMV
        ELSE
          TOTOX = PO30 * ONE_OVER_XL
          RATE = -DSIVDT( 2 ) / ( TOTOX * TSIV )
          DTW( 2 ) = 0.01 / ( RATE * MAX( TOTOX, TSIV ) )
        END IF

C...Calculate sulfur iv oxidation rate due to 02 catalyzed by Mn++
C...  and Fe+++  See Table IV Walcek & Taylor ( 1986)

        IF ( BB .GE. 4.0 )  THEN  ! 4.0  < pH

          IF ( SIV .LE. 1.0E-5 ) THEN
            DSIVDT( 3 ) = -5000.0 * MN * HSO3
          ELSE IF ( SIV .GT. 1.0E-5 ) THEN
            DSIVDT( 3 ) = -( 4.7 * MN * MN / AC
     &                  + 1.0E7 * FE * SIV * SIV )
          END IF  ! end of first pass through SIV conc.

        ELSE          ! pH , + 4.0

	  IF ( SIV .LE. 1.0E-5 ) THEN
            DSIVDT( 3 ) = -3.0 * ( 5000.0 * MN * HSO3
     &                  + 0.82 * FE * SIV / AC )
          ELSE
            DSIVDT( 3 ) = -( 4.7 * MN * MN / AC
     &                  + ( 0.82 * FE * SIV / AC )
     &                  * ( 1.0 + 1.7E3 * MN**1.5 / ( 6.3E-6 + FE ) ) )
          END IF ! end of second pass through SIV conc.

        END IF  ! end of pass through pH

        IF ( ( DSIVDT( 3 ) .EQ. 0.0 ) .OR. ( TSIV .LE. 1.0E-30 ) ) THEN
          DTW( 3 ) = DTRMV
        ELSE
          RATE = -DSIVDT( 3 ) / TSIV
          DTW( 3 ) = 0.1 / RATE
        END IF

C...Calculate sulfur oxidation rate due to MHP

        DSIVDT( 4 ) = -RMHP * AC * MHPL * HSO3
        IF ( ( DSIVDT( 4 ) .EQ. 0.0 ) .OR. ( TSIV .LE. 1.0E-30 ) ) THEN
          DTW( 4 ) = DTRMV
        ELSE
          TOTOX = PMHP0 * ONE_OVER_XL
          RATE = -DSIVDT( 4 ) / ( TOTOX * TSIV )
          DTW( 4 ) = 0.1 / ( RATE * MAX( TOTOX, TSIV ) )
        END IF

C...Calculate sulfur oxidation due to PAA

        DSIVDT( 5 ) = -RPAA * HSO3 * PAAL * ( AC + 1.65E-5 )
        IF ( ( DSIVDT( 5 ) .EQ. 0.0 ) .OR. ( TSIV .LE. 1.0E-30 ) ) THEN
          DTW( 5 ) = DTRMV
        ELSE
          TOTOX = PPAA0 * ONE_OVER_XL
          RATE = -DSIVDT( 5 ) / ( TOTOX * TSIV )
          DTW( 5 ) = 0.1 / ( RATE * MAX( TOTOX, TSIV ) )
        END IF

C...Calculate total sulfur iv oxidation rate

        DSIVDT( 0 ) = 0.0
        DO IOX = 1, NUMOX
          DSIVDT( 0 ) = DSIVDT( 0 ) + DSIVDT( IOX )
        END DO

C...Calculate a minimum time step required

        DTW( 0 ) = MIN( DTW( 1 ), DTW( 2 ), DTW( 3 ),
     &                  DTW( 4 ), DTW( 5 ) )

C...check for large time step

	IF ( DTW( 0 ) .GT. 8.0E+37 ) THEN
          write(iout,'(/,A)') 
     &          'Warning in RADM aqueous chemistry (AQCHEM)'
c          WRITE(idiag,1001) 
c     &          PRCRATE, DSIVDT(0), TS6, DTW(0), CTHK1, WTAVG
          WRITE(idiag,1001) 
     &          DSIVDT(0), TS6, DTW(0)
        ELSE

C...calculate the change in sulfur iv for this time step

c60        DTS6 = ABS( DTW( 0 ) * ( -DSIVDT( 0 ) - TS6 * PRCRATE
c     &         / ( 3.6 * CTHK1 * WTAVG * 1000.0 ) ) )
60        DTS6 = ABS( DTW( 0 ) * DSIVDT( 0 ) )

C...If DSIV(0), sulfur iv oxidized during this time step would be
C... less than 5% of sulfur oxidized since time 0, then double DT

          IF ( DTW( 0 ) .LE. TAUCLD ) THEN
            IF ( DTS6 .LT. 0.05 * TS6 ) THEN
              DTW( 0 ) = DTW( 0 ) * 2.0
	      GO TO 60
            END IF
          END IF
        END IF
        DTW( 0 ) = MIN( DTW( 0 ), DTRMV )

C...If the total time after this time increment will be greater than
C...  TAUCLD sec., then set DTW(0) so that total time will be TAUCLD

        IF ( TIMEW + DTW( 0 ) .GT. TAUCLD ) DTW( 0 ) = TAUCLD - TIMEW
        IF ( TS6 .LT. 1.0E-11 ) DTW( 0 ) = TAUCLD - TIMEW
c        IF ( ITERAT .GT. 100 ) DTW( 0 ) = TAUCLD - TIMEW
        IF ( ITERAT .GT. 100 ) then
          write(idiag,'(a,4i4,2f7.1)') 
     &    'AQCHEM: iterat > 100:', igrd, iaq, jaq, kaq, timew, taucld
          DTW( 0 ) = TAUCLD - TIMEW
        endif
C...Set DSIV(I), I = 0,NUMOX, the amount of S(IV) oxidized by each
C... individual oxidizing agent, as well as the total.

        DO IOX = 0, NUMOX
          DS4( IOX ) = DS4( IOX ) + DTW( 0 ) * DSIVDT( IOX )
        END DO

C...Compute depositions and concentrations for each species
c
c        DO LIQ = 1, NLIQS
c          WETDEP( LIQ ) = WETDEP( LIQ )
c     &                  + PRCRATE * LIQUID( LIQ ) * DTW( 0 ) * WCAVG
c     &                  * 1000.0 / ( 3600.0 * WTAVG * 1000.0 )
c        END DO

        TIMEW = TIMEW + DTW( 0 )

C...Return to make additional calculations

        GO TO 20
      END IF

C...At this point, TIMEW=TAUCLD
C...  compute the scavenging coefficient for SO4 which will be used for
C...  scavenging aerosol number in the accumulation and coarse mode
c
c      DEPSUM = ( WETDEP( 6 ) + WETDEP( 7 ) ) * XC1
c
c      IF ( ( TS6ACCA + TS6AKNA - DS4( 0 ) ) .NE. 0.0 ) THEN
c        BETASO4 = DEPSUM / ( ( TS6ACCA + TS6AKNA - DS4( 0 ) ) * TAUCLD )
c      ELSE
c        BETASO4 = 0.0
c      END IF
c
C...Compute the output concentrations and wetdeposition amounts

      TOTAMM = ( PNH3F  + ( NH4 + NH3L  ) * XL ) * RECIPAP1
      TOTNIT = ( PHNO3F + ( NO3 + HNO3L ) * XL ) * RECIPAP1

c
C...gas-phase species wet deposition (mm mol/lit)
c
c      GASWDEP( LSO2   ) = WETDEP(  8 ) + WETDEP(  9 ) + WETDEP( 10 )
c      GASWDEP( LNH3   ) = WETDEP( 15 )
c      GASWDEP( LH2O2  ) = WETDEP( 17 )
c      GASWDEP( LO3    ) = WETDEP( 18 )
c      GASWDEP( LCO2   ) = WETDEP( 11 ) + WETDEP( 12 ) + WETDEP( 13 )
c      GASWDEP( LFOA   ) = WETDEP( 22 ) + WETDEP( 23 )
c      GASWDEP( LMHP   ) = WETDEP( 24 )
c      GASWDEP( LPAA   ) = WETDEP( 25 )
c      GASWDEP( LHNO3  ) = WETDEP( 32 )
c      GASWDEP( LN2O5  ) = 0.0
c      GASWDEP( LH2SO4 ) = 0.0
c
C...gas concentrations (mol/molV)

      GAS( LSO2   ) = ( PSO2F   + XL *  SIV )   * RECIPAP1
      GAS( LNH3   ) = FNH3 * TOTAMM
      GAS( LH2O2  ) = ( PH2O2F  + XL *  H2O2L ) * RECIPAP1
      GAS( LO3    ) = ( PO3F    + XL *  O3L )   * RECIPAP1
      GAS( LCO2   ) = ( PCO2F   + XL *  CO2L )  * RECIPAP1
      GAS( LFOA   ) = ( PFOAF   + XL * ( FOAL
     &              +  HCO2 ) ) * RECIPAP1
      GAS( LMHP   ) = ( PMHPF   + XL *  MHPL )  * RECIPAP1
      GAS( LPAA   ) = ( PPAAF   + XL *  PAAL )  * RECIPAP1
      GAS( LHNO3  ) = FHNO3 * TOTNIT
      GAS( LN2O5  ) = 0.0 ! assume all into aerosol
      GAS( LH2SO4 ) = 0.0 ! assume all into aerosol

C...aerosol species wet deposition (mm mol/lit)
C...  there is no wet deposition of aitken particles, they attached
C...  to the accumulation mode particles
c
c      AERWDEP( LSO4AKN ) = 0.0
c      AERWDEP( LSO4ACC ) = WETDEP(  6 ) + WETDEP(  7 )
c      AERWDEP( LNH4AKN ) = 0.0
c      AERWDEP( LNH4ACC ) = WETDEP(  2 )
c      AERWDEP( LNO3AKN ) = 0.0
c      AERWDEP( LNO3ACC ) = WETDEP( 14 ) * FRACACC
c      AERWDEP( LNO3COR ) = WETDEP( 14 ) * FRACCOR
c      AERWDEP( LORGAKN ) = 0.0
c      AERWDEP( LORGACC ) = WETDEP( 27 )
c      AERWDEP( LPRIAKN ) = 0.0
c      AERWDEP( LPRIACC ) = WETDEP( 28 )
c      AERWDEP( LPRICOR ) = WETDEP( 33 )
c      AERWDEP( LNACL   ) = WETDEP(  4 )
c      AERWDEP( LA3FE   ) = WETDEP( 19 )
c      AERWDEP( LB2MN   ) = WETDEP( 20 )
c      AERWDEP( LCACO3  ) = WETDEP(  3 )
c      AERWDEP( LMGCO3  ) = WETDEP( 29 )
c      AERWDEP( LKCL    ) = WETDEP( 30 )
c      AERWDEP( LNUMAKN ) = 0.0
c      AERWDEP( LNUMACC ) = AEROSOL( LNUMACC ) * AIRM
c     &                   * ( 1.0 - EXP( -BETASO4 * TAUCLD ) )
c      AERWDEP( LNUMCOR ) = AEROSOL( LNUMCOR ) * AIRM
c     &                   * ( 1.0 - EXP( -BETASO4 * TAUCLD ) )
c      AERWDEP( LSRFAKN ) = 0.0
c      AERWDEP( LSRFACC ) = 0.0
c
cC...aerosol concentrations (mol/molV)
c
camx
      AEROSOL( LSO4    ) = TS6    * XL * RECIPAP1
      AEROSOL( LNH4    ) = TOTAMM * (1.0 - FNH3)
      AEROSOL( LNO3    ) = TOTNIT * (1.0 - FHNO3)
      AEROSOL( LNACL   ) = NA     * XL * RECIPAP1
      AEROSOL( LA3FE   ) = FE     * XL * RECIPAP1
      AEROSOL( LB2MN   ) = MN     * XL * RECIPAP1
      AEROSOL( LCACO3  ) = CA     * XL * RECIPAP1
      AEROSOL( LMGCO3  ) = MG     * XL * RECIPAP1
      AEROSOL( LKCL    ) = K      * XL * RECIPAP1
camx
c      AEROSOL( LSO4AKN ) = AEROSOL( LSO4AKN ) * EXP( -ALFA3 * TAUCLD )
c      AEROSOL( LSO4ACC ) = TS6    * XL * RECIPAP1
c      AEROSOL( LNH4AKN ) = AEROSOL( LNH4AKN ) * EXP( -ALFA3 * TAUCLD )
c      AEROSOL( LNH4ACC ) = FNH4ACC * TOTAMM
c      AEROSOL( LNO3AKN ) = AEROSOL( LNO3AKN ) * EXP( -ALFA3 * TAUCLD )
c      AEROSOL( LNO3ACC ) = FNO3ACC * TOTNIT
c      AEROSOL( LNO3COR ) = FNO3COR * TOTNIT
c      AEROSOL( LORGAKN ) = AEROSOL( LORGAKN ) * EXP( -ALFA3 * TAUCLD )
c      AEROSOL( LORGACC ) = ORGN   * XL * RECIPAP1
c      AEROSOL( LPRIAKN ) = AEROSOL( LPRIAKN ) * EXP( -ALFA3 * TAUCLD )
c      AEROSOL( LPRIACC ) = PRIM   * XL * RECIPAP1
c      AEROSOL( LPRICOR ) = PRIMCOR* XL * RECIPAP1
c      AEROSOL( LNACL   ) = NA     * XL * RECIPAP1
c      AEROSOL( LA3FE   ) = FE     * XL * RECIPAP1
c      AEROSOL( LB2MN   ) = MN     * XL * RECIPAP1
c      AEROSOL( LCACO3  ) = CA     * XL * RECIPAP1
c      AEROSOL( LMGCO3  ) = MG     * XL * RECIPAP1
c      AEROSOL( LKCL    ) = K      * XL * RECIPAP1
c      AEROSOL( LNUMAKN ) = AEROSOL( LNUMAKN ) * EXP( -ALFA0 * TAUCLD )
c      AEROSOL( LNUMACC ) = AEROSOL( LNUMACC ) * EXP( -BETASO4 * TAUCLD )
c      AEROSOL( LNUMCOR ) = AEROSOL( LNUMCOR ) * EXP( -BETASO4 * TAUCLD )
c
cC...compute the final accumulation aerosol 3rd moment
c
c      M3NEW = ( AEROSOL( LSO4ACC ) * SGRAERMW( LSO4ACC ) / 1.8e6
c     &      +   AEROSOL( LNH4ACC ) * SGRAERMW( LNH4ACC ) / 1.8e6
c     &      +   AEROSOL( LNO3ACC ) * SGRAERMW( LNO3ACC ) / 1.8e6
c     &      +   AEROSOL( LORGACC ) * SGRAERMW( LORGACC ) / 2.0e6
c     &      +   AEROSOL( LPRIACC ) * SGRAERMW( LPRIACC ) / 2.2e6 )
cCCC     &      * 6.0 / PI      ! cancels out in division below
c
c      AEROSOL( LSRFAKN ) = AEROSOL( LSRFAKN ) * EXP( -ALFA2 * TAUCLD )
c      AEROSOL( LSRFACC ) = AEROSOL( LSRFACC )
c     &                   * ( EXP( -BETASO4 * TAUCLD * ONETHIRD ) )
c     &                   * ( M3NEW / MAX( M3OLD, 1.0E-30) ) ** TWOTHIRDS
c
cC...store the amount of hydrogen deposition
c
c      HPWDEP = WETDEP( 1 )
c
      RETURN

C...formats

c1001  FORMAT (1X,'STORM RATE=', F6.3, 'DSIVDT(0) =', F10.5,
c     &       'TS6=', F10.5, 'DTW(0)=', F10.5, 'CTHK1=', F10.5,
c     &       'WTAVG=', F10.5)
1001  FORMAT (1X, 'DSIVDT(0) =', F10.5,
     &       'TS6=', F10.5, 'DTW(0)=', F10.5)

      END
