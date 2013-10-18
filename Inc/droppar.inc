c UPDATED DROPPAR.INC FOR THE SIZE RESOLVED MODEL
c LAST UPDATE : 10 June 1998
c*************************************************************************
c                                DROPPAR.INC
c*************************************************************************
c
c                AQUEOUS-PHASE PARAMETERS AND VARIABLES
c
c NSP : NUMBER OF TOTAL PARTICLE SECTIONS
c
       parameter (nsp_aq = nsect)
c
c NBOUNDS : NUMBER OF SECTIONAL BOUNDARIES
c
c AQUEOUS-PHASE COMPONENTS
c
c   IMPORTANT : ALL COMPONENTS HAVE THE SAME POSITIONS IN      
c               BOTH AEROSOL AND AQUEOUS MATRICES
c               NEVER CHANGE THIS CONVENTION BECAUSE AQMAIN
c               DEPENDS ON IT
c
      parameter (ksod = nas)               ! Na(+)
      parameter (khyd = nah)               ! H(+)
      parameter (kamm = naa)               ! NH4(+)
      parameter (knit = nan)               ! NO3(-)
      parameter (kchl = nac)               ! Cl(-)
      parameter (ksvi = na4)               ! S(VI)
      parameter (kwat = naw)               ! H2O
      parameter (kec1 = nae)               ! EC
      parameter (koc  = nao)               ! OC
      parameter (kcru = nar)               ! Crustal
c      parameter (khso5 = 1)                ! HSO5-
c      parameter (khmsa = 2)                ! HMSA
c      parameter (kform = 3)                ! Formic Acid
c
c GASES IN LOCAL ARRAY
c
c     INCORRECT  ******* TO BE FIXED *************      
      parameter (ngso2 = 11)
      parameter (ngh2o2 = 12)
      parameter (nghcho = 13)
      parameter (nghcooh = 14)
      parameter (nghno2 = 15)
      parameter (ngno = 16)
      parameter (ngno2 = 17)
      parameter (ngo3 = 18)
      parameter (ngpan = 19)
      parameter (ngoh = 20)
      parameter (ngho2 = 21)
      parameter (ngno3 = 22)
      parameter (ngch3o2 = 23)
      parameter (ngch3o2h = 24)
      parameter (ngch3oh = 25)
      parameter (ngch3co3h = 26)
c                                         
c     NUMBER OF EQUATIONS FOR AQUEOUS-PHASE CHEMISTRY SOLUTION
c
       parameter (meqn1 = 13)
       parameter (meqn2 = 21)
c      
c     ACTIVATION DIAMETER (Dry)
c
      real*4 dactiv,dsep,diamet1,diamet2,avdiam
      parameter (dactiv = 0.7e-6)       ! in m
c
c     SEPARATION DIAMETER FOR TWO SECTIONS (DRY)
c 
      parameter (dsep = 2.5e-6)         ! in m
c
c     WET DIAMETERS OF THE TWO SECTIONS
c
      parameter (diamet1 = 19.e-6)      ! in m
      parameter (diamet2 = 29.e-6)      ! in m
c     If bulk, diameter = avdiam
      parameter (avdiam = 20.e-6)     
c
c     CHOICE OF EXPRESSION FOR IRON CHEMISTRY
c               = 0 (No iron/manganese chemistry)
c          KIRON = 1 (Phenomenological, Martin et al., 1991)
c                = 2 (Martin, 1984)
c
      parameter (kiron = 2)            ! was 1
c
c     CHOICE OF TURNING ON OR OFF RADICAL CHEMISTY
c     (IT IS BETTER TO TURN IT OFF DURING THE NIGHT)
c
      parameter (iradical = 0)       ! kmf (07/24/02)
c
c     CHOICE OF TURNING OFF CHLORINE CHEMISTRY
c
      real*4 chlorine
      parameter (chlorine = 0.0)
c      
c     CHOICE OF TURNING OFF NO2 RXN (104) ! kmf (08/27/02)
c
      real*4 nitrogen
      parameter (nitrogen = 0.0)      ! uncertain rxn rate - producing too much SO4 here
c
c     PARAMETER FOR SCALING OF PHOTOLYSIS RATES
c
      real*4 photo
      parameter (photo = 1.0)  
c      
c     CALCULATE THE DISTRIBUTION FACTORS BASED ON THE LIQUID
c     WATER CONTENT
c     FDIST(isect) : IS THE FRACTION OF THE MASS CHANGE THAT WILL
c                    GO TO THE CORRESPONDING SECTION
c 
      double precision fdist(nsect), fdist2(nsect)
      common /fd/ fdist, fdist2
c for sections from 0.01 to 10 um
c      data fdist /0.,0.,0.,0.0,0.13, 0.27, 0.1, 0.2, 0.2, 0.1/
c      data fdist2/0.,0.,0.,0.0,0.267,0.533,0.2, 0.4, 0.4, 0.2/      
c for sections from 0.04 to 40 um
cgem      data fdist /0.0, 0.0, 0.0, 0.1, 0.3, 0.1, 0.35, 0.15, 0.0, 0.0/
cgem      data fdist2/0.0, 0.0, 0.0, 0.2, 0.6, 0.2, 0.70, 0.30, 0.0, 0.0/
c
c     Fraction of crustal material that is alkaline 
c
      real*4 caratio                   ! kmf (08/01/02)
      parameter (caratio = 0.0)        ! was 0.1  
c
c     
c     Size resolution 1=bulk,
c                     2=2-section sr
c                     3=variable sr)
c     bulk treats all droplets as homogeneous
c     2-section size resolved divides droplets into 2 separate 
c     	homogeneous sections
c     variable size-resolved - executes either the bulk or size-resolved
c	calculations based on initial conditions (NH3, HNO3, SO2, H2O2, 
c 	lwc, and caratio)
c
      parameter (ires = 3)
c
c     fraction of liquid water content that goes to each s.r. section
c
      real*4 frac1,frac2
      parameter (frac1 = 0.8)               ! fraction of lwc in sect. 1
      parameter (frac2 = 0.2)               ! fraction of lwc in sect. 2
c      
c
c     ASSUMPTION : Fe(3+) and Mn(2+) = 0.003%, 0.001% of Crustal mass
c      
      real*4 firon,fman
      parameter (firon = 0.00003)
      parameter (fman = 0.00001)
CKF      
      real*4 minso2
      parameter (minso2 = 1.5e-4)
