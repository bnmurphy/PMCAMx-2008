c     This file contains the necessary parameteres for the new generation
c     multicomponent aerosol dynamics model (MADM)
c *************************************************************************
c     WRITTEN BY DR. CHRISTODOULOS PILINIS
c              December, 1996
c
c     Modified for CAMx and secondary organic aerosols by Fred Lurmann May 2000
c
c     03/09/03: bkoo
c               - commented out parameter nsoap (nsoap is now assigned in soap.com)
c *************************************************************************
c
      include 'section.inc'
c  ** changes in nsec necesitate changes in formats for subroutine SPCPRINT **
      parameter(neqsec=4)             ! Number of equilibrium sections
cbk      parameter(neqsec=2)
      parameter(nsecd=nsec-neqsec)    ! Number of dynamic sections
      parameter(nexti=6)              ! Number of external inorganic aerosol species
c
      parameter(imovsec=1)         ! 0 is fixed sections 1 is moving sections
c      parameter(imovsec=0)         ! 0 is fixed sections 1 is moving sections
c
c  external inorganic species are those to be used by the 3D-gas-aerosol model
c  while internal inorganic species are the detailed species, handled by MADM

      parameter( KH2O=1 )       ! pointer for external aerosol water
      parameter(  KNa=2 )       ! pointer for external aerosol sodium
      parameter( KSO4=3 )       ! pointer for external aerosol sulfate
      parameter( KNO3=4 )       ! pointer for external aerosol nitrate
      parameter( KNH4=5 )       ! pointer for external aerosol ammonium
      parameter(  KCL=6 )       ! pointer for external aerosol chlorides
c  
c     BNM 07/21/09  CHANGED THESE POINTERS TO ADD NTSOA
      parameter (KAPO1 = 7)!pointer for external FRESH POA
      parameter (KAPO2 = 8)!pointer for external FRESH POA
      parameter (KAPO3 = 9)!pointer for external FRESH POA
      parameter (KAPO4 =10)!pointer for external FRESH POA
      parameter (KAPO5 =11)!pointer for external FRESH POA
      parameter (KAPO6 =12)!pointer for external FRESH POA
      parameter (KAPO7 =13)!pointer for external FRESH POA
      parameter (KAPO8 =14)!pointer for external FRESH POA
c
      parameter (KA0O1 =15)!pointer for external OXIDIZED POA
      parameter (KA0O2 =16)!pointer for external OXIDIZED POA
      parameter (KAOO3 =17)!pointer for external OXIDIZED POA
      parameter (KAOO4 =18)!pointer for external OXIDIZED POA
      parameter (KAOO5 =19)!pointer for external OXIDIZED POA
      parameter (KAOO6 =20)!pointer for external OXIDIZED POA
      parameter (KAOO7 =21)!pointer for external OXIDIZED POA
      parameter (KAOO8 =22)!pointer for external OXIDIZED POA
c
      parameter (KABS1 =23)!pointer for external BIOGENIC SOA
      parameter (KABS2 =24)!pointer for external BIOGENIC SOA
      parameter (KABS3 =25)!pointer for external BIOGENIC SOA
      parameter (KABS4 =26)!pointer for external BIOGENIC SOA
      parameter (KABS5 =27)!pointer for external BIOGENIC SOA
c
      parameter (KAAS1 =28)!pointer for external ANTHROPOGENIC SOA
      parameter (KAAS2 =29)!pointer for external ANTHROPOGENIC SOA
      parameter (KAAS3 =30)!pointer for external ANTHROPOGENIC SOA
      parameter (KAAS4 =31)!pointer for external ANTHROPOGENIC SOA
      parameter (KAAS5 =32)!pointer for external ANTHROPOGENIC SOA
c
      parameter (KANS1 =33)!pointer for external NON-TRADITIONAL SOA
      parameter (KANS2 =34)!pointer for external NON-TRADITIONAL SOA
      parameter (KANS3 =35)!pointer for external NON-TRADITIONAL SOA
      parameter (KANS4 =36)!pointer for external NON-TRADITIONAL SOA
      parameter (KANS5 =37)!pointer for external NON-TRADITIONAL SOA
      parameter (KANS6 =38)!pointer for external NON-TRADITIONAL SOA
      parameter (KANS7 =39)!pointer for external NON-TRADITIONAL SOA
      parameter (KANS8 =40)!pointer for external NON-TRADITIONAL SOA
c
      parameter (KPOM  =41)!pointer for external Primary Organic Aerosol
      parameter (KEC  = 42)!pointer for external Elemental Carbon Aerosol
      parameter (KCRUS =43)!pointer for external Crustal Aerosol 

c
c
c      parameter(ninti=18)       ! Number of internal inorganic aerosol species
      parameter(ninti=19)       ! bkoo (02/14/02)

c      1 = H2O
c      2 = H+
c      3 = Na+
c      4 = NH4+
c      5 = Cl-
c      6 = SO4--
c      7 = HSO4-
c      8 = NO3-
c      9 = NaCl
c     10 = Na2SO4
c     11 = NaNO3
c     12 = (NH4)2SO4
c     13 = NH4NO3
c     14 = NH4Cl
c     15 = H2SO4 (aq)
c     16 = NH4HSO4
c     17 = NaHSO4
c     18 = (NH4)3H(SO4)2
c     19 = OH-                  ! bkoo (02/14/02)
c
       parameter(ninert=2)    	! Number inert aerosol species (EC & Crustal)
       parameter(norg=35)       ! Number of organic aerosol species (Both primary and secondary) BNM 07/21/09
cbk       parameter(nsoap=4)    ! Number or secondary organic aerosol species
cbk       parameter(ngas=4+nsoap)        ! Number of gases of interest, 4 inorganic + nsoap
       parameter(ngas=4+34) ! bkoo (03/09/03)  [4+CONDENSIBLE ORGS]
       parameter(INH3  =1)            ! Pointer of NH3 in the gas phase
       parameter(IHNO3 =2)            ! Pointer of HNO3 in the gas phase
       parameter(IH2SO4=3)            ! Pointer of H2SO4 in the gas phase
       parameter(IHCL  =4)            ! Pointer of HCl in the gas phase
c
      parameter (ICPO1 = 5)!pointer for external FRESH POA
      parameter (ICPO2 = 6)!pointer for external FRESH POA
      parameter (ICPO3 = 7)!pointer for external FRESH POA
      parameter (ICPO4 = 8)!pointer for external FRESH POA
      parameter (ICPO5 = 9)!pointer for external FRESH POA
      parameter (ICPO6 =10)!pointer for external FRESH POA
      parameter (ICPO7 =11)!pointer for external FRESH POA
      parameter (ICPO8 =12)!pointer for external FRESH POA
c
      parameter (ICOO1 =13)!pointer for external OXIDIZED POA
      parameter (ICOO2 =14)!pointer for external OXIDIZED POA
      parameter (ICOO3 =15)!pointer for external OXIDIZED POA
      parameter (ICOO4 =16)!pointer for external OXIDIZED POA
      parameter (ICOO5 =17)!pointer for external OXIDIZED POA
      parameter (ICOO6 =18)!pointer for external OXIDIZED POA
      parameter (ICOO7 =19)!pointer for external OXIDIZED POA
      parameter (ICOO8 =20)!pointer for external OXIDIZED POA
c
      parameter (ICBS1 =21)!pointer for external BIOGENIC SOA
      parameter (ICBS2 =22)!pointer for external BIOGENIC SOA
      parameter (ICBS3 =23)!pointer for external BIOGENIC SOA
      parameter (ICBS4 =24)!pointer for external BIOGENIC SOA
      parameter (ICBS5 =25)!pointer for external BIOGENIC SOA
c
      parameter (ICAS1 =26)!pointer for external ANTHROPOGENIC SOA
      parameter (ICAS2 =27)!pointer for external ANTHROPOGENIC SOA
      parameter (ICAS3 =28)!pointer for external ANTHROPOGENIC SOA
      parameter (ICAS4 =29)!pointer for external ANTHROPOGENIC SOA
      parameter (ICAS5 =30)!pointer for external ANTHROPOGENIC SOA
c
      parameter (ICNS1 =31)!pointer for external NON-TRADITIONAL SOA
      parameter (ICNS2 =32)!pointer for external NON-TRADITIONAL SOA
      parameter (ICNS3 =33)!pointer for external NON-TRADITIONAL SOA
      parameter (ICNS4 =34)!pointer for external NON-TRADITIONAL SOA
      parameter (ICNS5 =35)!pointer for external NON-TRADITIONAL SOA
      parameter (ICNS6 =36)!pointer for external NON-TRADITIONAL SOA
      parameter (ICNS7 =37)!pointer for external NON-TRADITIONAL SOA
      parameter (ICNS8 =38)!pointer for external NON-TRADITIONAL SOA
c
       parameter(nsp=nexti+norg+ninert)
       parameter(ntotal=nsec*nsp+ngas)
       parameter(ntotald=(nsec-neqsec)*nsp+ngas)
       parameter(ntotale=neqsec*nsp+ngas)
       parameter(nsxgas=nsec*ngas)
       parameter(nsxint=nsec*ninti)
       parameter(naer=ntotal-ngas)
c
c      SWITCHES FOR TURNING PROCESSES ON OR OFF
c
       parameter (inucl = 1)           ! 1 is on, 0 is off
       parameter (icoag = 1)           ! 1 is on, 0 is off
c
c      ADDITIONAL RESOLUTION FOR COAGULATION CALCULATIONS
c
       parameter (nres = 3)
c
       implicit double precision (a-h,o-z)
       logical dry(nsec)              ! true, for dry aerosol       

       dimension diffus(nsp)          ! diffusion coefficients of condensible gases (SI units)
       dimension delta(nsp)           ! accomodation coefficients of condensible gases 
       dimension vel(nsp)             ! Mean speed of condensible gas molecules
       double precision lamda(nsp)    ! Mean free path of condensible gases
       dimension dsec(nsecp1)          ! Sectional diameters (1-->lower limit, nsec+1-->upper limit) 
c                                                        units=umeters
       dimension qt(nsec)              ! Total aerosol concentration at each section
       dimension qtt(nsp)              ! Total concentrations of aerosol species (ug/m3)
       dimension qn(nsec)              ! 0th moment (number) of total mass for seciton (qt/dsec**3)
       dimension ps(nsec,nsp)         ! Gas phase partial pressures of 4 inorganic and norg organic
c                                        species at the particle surface (Pa)
c
       dimension qsource(naer)              ! Aerosol module master array for primary aer. emissions of
                                   ! aerosol species in ugr/m3/sec.
        dimension gsource(ngas)         ! Gas phase source in ug/m3/sec.
c
       dimension c(nsec, ninti)       ! Concentrations of internal inorganic aerosol species (ugr/m3) 
       dimension emw(nsp)             ! Molecular weights of external inorganics + organics
       dimension gmw(ngas)            ! Molecular weights of condensible gases ( inorganics + organics )
       double precision intmw(ninti)  ! Molecular weights of internal inorganics
       dimension hi(nsec,ngas)        ! Hi coefficients of condensible gases ( inorganics + organics )
c                                       from Equation 8 of Pilinis et al. (1987)
       double precision qtot0(nexti)  ! initial aerosol plus gas concentrations of inorganic ions
c
c      parameter(lrw=22 + ntotal * ( 7 + ntotal + 9))  ! added 7 so that the min is 16 as required by lsoda
       parameter(lrw=22 +  9*ntotal + ntotal**2 )      ! for lsode
       parameter(liw=20+ntotal)
       dimension rwork(lrw)                ! work array for lsoda
       dimension iwork(liw)                ! integer work array for lsoda
       dimension dsecf(nsecp1)             ! sectional diameters of the fixed grid
       character*4 aerm			! aerosol module
					!  EQUI for equilibrium model
					!  MADM for dynamic model (MADM)
					!  HYBR for hybrid model
       character*20 ver	 		! version info.

        common /module/ pi, dsec, ps, qsource, c, emw, intmw, qtt, rh,
     & temp, dry, rgas, diffus, vel, delta, lamda,hi,qt,gmw,qn,tcom,
     & gsource, dsecf
c
      common / tinys / tinys  ! minimum non-zero concentration (ug/m3)
      common / massbal / qtot0,ifgprt,ifgbal
      common / aerm / aerm
      common / ver / ver
      common / dt / dt,pres
c
c     added by bkoo (09/19/00)
      common / arraysize / ntotalx, nsecx, ntotalx2, nsecx2
c
c     added by bkoo (02/26/02)
      parameter(atoli=1.0d-6) ! absolute tolerance for aerosols
      parameter(rtoli=1.0d-3) ! relative tolerance for all species
      parameter(govra=1.0d-2) ! ratio of PPM to UG/M3 for most species

c     added by bkoo (03/05/02)
      logical isfirst         ! first call to ISRPIA of the dynamic step?
      logical dold(nsec)      ! previous state of the aerosol (dry?)
      integer ims(nsec)       ! 0 - deliquescent, 1 - metastable
      integer icdry(nsec)     ! counter for WET/DRY oscillation
      common / ifms / isfirst, dold, ims, icdry
