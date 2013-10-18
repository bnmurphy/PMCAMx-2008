      subroutine cpamech3( r, nr, pa, npa, npa_init )
c
c-----CAMx v4.02 030709
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
C***********************************************************************
C  Called by CHEM, uses integrated reaction rates (R) to calculate
C  important chemical process rates (PA)
C
C    *** NOTE:  This code is hardwired for the CAMx CB4 vers. 3 ***
C
C  Arguments:  R should be the average rate during the most recent
C              chemistry time step, multipled by the time step, eg:
C              the chem solver calculates 0.5*dt*(R+Rold) and passes it
C              to this routine. It's OK to use simply R*dt, but this
C              is a slightly less accurate form of the integration.
C              Important reaction rates are summed below to calculate PA 
C              which is passed back to the Chem routine and accumulated for 
C              one hour and then output. 
C              NR is the number of reactions.
C              NPA is the number of chemical process rates calculated below.
C              RJNO2 is the photo rate constant multiplied by dt, dimensionless.
C              R is assumed to be ppm.
C              PA is converted to ppb.
C               
C** REVISION HISTORY:
C**   NO.    DATE   WHO  WHAT
C**   ---    -----  ---  --------------------------------------------
C**   01     10/00  GST  Created for use with the CB4 photochemical mechanism.
C**                      This is hardwired for the CAMx CB4 vers 3.
C**                      
C**                         
C***********************************************************************
C
      include 'camx.prm'
      include 'filunit.com'
      include 'tracer.com'
      include 'procan.com'
C
C... Local variable declarations
C
      INTEGER  NR, NPA, NPA_INIT
      REAL     R(NR), PA(NPA)
      REAL     cyc_HONO_p, cyc_HONO_l, cyc_H2O2_p, cyc_H2O2_l
      REAL     cyc_HNO4_p, cyc_HNO4_l, cyc_PAN_p,  cyc_PAN_l
      REAL     cyc_O3OH_HO2_p, cyc_O3OH_HO2_l
C
C     temporary sums for budget calculations
      REAL   OHwHC, HO2new, RO2new
      REAL   OHnew_isop, HO2new_isop, RO2new_isop
      REAL   OHterm, other_OH_prop, tmpprop, prod_hcho_from_isop
C
C***********************************************************************
C
C  { calculation of these PA outputs is currently hardwired below }
C
C     Ox Production
C     Ox Chem Destruction
C
C  { Radical Initiation }

C     New OH from O1D+H2O
c     New OH from H2O2, HNO3, HONO, PAA, OP1, OP2, O3+HC (except isoprene)
C     New HO2 from HCHO
C     New HO2 Production (Total)
C     New RO2 Production (Total)
C     New HOx from isoprene
C
C  { Radical Propagation }
C
C     OH+CO and OH+CH4
C     OH+ISOP
C     Ox+ISOP
C     OH reacted with VOC
C     other OH propagation rxns
C     Total HO2 Production
C     Total RO2 Prod
c     NO2 from HO2
C     OH from HO2
C     NO2 from RO2
C
C  { Radical Termination }
C
C     Total OH production
C     OH termination
C     HO2 termination   
C     RO2 termination  
C
C  {  NOz chemistry }
C
C     OH+NO2=HNO3
C     NO3+HC=HNO3
C     N2O5+H2O= HNO3
C     HNO3 reacted
C     net PAN Prod
C     net PAN Loss
C     production of ONIT
C
C***********************************************************************
C
c
C
C.....{ Ox Prod }
      nn =  1
      ptname(nn)  = 'OxProd'
C
C { Ox =O3, NO2, NO3*2, O3P, O1D, HNO3, PAN, HNO4, N2O5*3, ORNIT }
C
      PA(nn) = 2.*R(20)  ! NO+NO=2*NO2
     &       +    R(24)  ! OH+HONO=NO2
     &       +    R(25)  ! HONO+HONO=NO+NO2 
     &       +    R(27)  ! OH+HNO3=NO3
     &       +    R(28)  ! HO2+NO=NO2+OH
     &       +    R(46)  ! C2O3+NO
     &       +    R(64)  ! TO2+NO=.9*NO2+.1*ORNIT
     &       +    R(79)  ! XO2+NO=NO2
     &       +    R(81) !{XO2N+NO}
C
C
C.....{ Ox Loss }
      nn = nn + 1
      ptname(nn)  = 'OxLoss'
C
      PA(nn) =2.0*R(4)   ! O3P+NO2=NO
     &       +    R(11)  ! O1D+H2O=2*OH
     &       +    R(12)  ! O3+OH=HO2
     &       +    R(13)  ! O3+HO2=OH
     &       +.22*R(14)  ! NO3=NO
     &       +2.0*R(16)  ! NO3+NO2=NO+NO2
     &       +    R(18)  ! N2O5+H2O=2*HNO3
     &       +    R(21)  ! no+no2+h2o=2.000*hono
     &       +    R(40)  ! form+o=oh+ho2+co 
     &       +    R(41)  ! HCHO+NO3=HO2+HNO3+CO
     &       +    R(42)  ! ALD2+O
     &       +    R(44)  ! ALD2+NO3=ACO3+HNO3+.4*CP
     &       +    R(56)  ! OLE+O
     &       +    R(58)  ! OLT+O3
     &       +    R(59)  ! OLE+NO3=products
     &       +    R(60)  ! ETH+O
     &       +    R(62)  ! OL2+O3=HCHO+.40*ORA1+.42*CO+.12*HO2
     &       +    R(67)  ! CSL+NO3=HNO3+XNO2+.50*CSL
     &       +    R(71)  ! OPEN+O3
     &       +    R(75)  ! O+ISOP
     &       +    R(77)  ! ISO+O3
     &       +    R(78)  ! ISO+NO3=0.8*ONIT+0.2*NO2
     &       +    R(93)  ! ISPD+O3
     &       +    R(94)  ! ISPD+NO3=.85*NTR+.15*HNO3
     &       + .2*R(96)  ! ISOP+NO2=0.8*NTR
C
C Note:  We will consider NTR a form of Ox so these
C        reactions are not sinks for Ox:
C     &       +    R(55)  ! ROR+NO2=NTR
C     &       +    R(68)  ! CRO+NO2=NTR
C
C
C  { Radical Reservoirs: HONO, HNO4, H2O2, PAN                              }
C  { we will call these termination rxns if there is net production         }
C  { during a time step, or initiation if there is net release of HOx.      }
C  { It may be useful to distinguish this initiation from NO3, O3 + HC rxns }
C
C  { 23} HONO=OH+NO
C  { 22} NO+OH=HONO
C
      sum =  R(22) - R(23)          ! HONO cycle
      if (sum.GT.0.) then
         cyc_HONO_p = sum
         cyc_HONO_l = 0.0
      else
         cyc_HONO_p =  0.0
         cyc_HONO_l = -sum
      end if
C 
C  { 29} HO2+NO2=HNO4
C  { 30} HNO4=HO2+NO2
C
      sum =  R(29) - R(30)           ! HNO4 cycle
      if (sum.GT.0.) then
         cyc_HNO4_p = sum
         cyc_HNO4_l = 0.0
      else
         cyc_HNO4_p =  0.0
         cyc_HNO4_l = -sum
      end if
C
C { 32} HO2+HO2=H2O2
C { 33} HO2+HO2+H2O=H2O2
C { 34} H2O2=2.00*OH
      sum =  R(32) + R(33) - R(34)   ! H2O2 cycle
      if (sum.GT.0.) then
         cyc_H2O2_p = 2.*sum
         cyc_H2O2_l = 0.0
      else
         cyc_H2O2_p =  0.0
         cyc_H2O2_l = -2.0*sum
      end if
C
C  { 47} ACO3+NO2=PAN
C  { 48} PAN=ACO3+NO2
C
      sum =  R(47) - R(48)            ! PAN Cycle
      if (sum.GT.0.) then
         cyc_PAN_p = sum
         cyc_PAN_l = 0.0
      else
         cyc_PAN_p =  0.0
         cyc_PAN_l = -sum
      end if
C                                           { O3-HOx Cycle }
C {12} O3+OH=HO2
C {13} O3+HO2=OH
      sum = R(12) - R(13) 
      if (sum.GT.0.) then
         cyc_O3OH_HO2_p = sum
         cyc_O3OH_HO2_l = 0.0
      else
         cyc_O3OH_HO2_p =  0.0
         cyc_O3OH_HO2_l = -sum
      end if
C
C
C
C
C      {      R a d i c a l    I n i t i a t i o n       }
C
C
C
C...  New OH from O1D+H2O
      nn = nn + 1
      ptname(nn)  = 'newOH_O1D'
C
      PA(nn) = 2.*R(11)   !O1D+H2O->2*OH
C
C
C
C     { calculate total rad init from isoprene to be used below }
      OHnew_isop =
     &       +0.266*R(77)  !O3+ISO =.066*HO2+.266*OH
     &       +0.268*R(93)  !O3+ISPD =.154*HO2+.268*OH+.114*C2O3
C
      HO2new_isop =
     &       + .25*R(75)    !ISO+O
     &       +.066*R(77)    !ISO+O3
     &       + .80*R(78)    !ISO+NO3
     &       +.154*R(93)    !ISPD+O3
     &       +.925*R(94)    !ISPD+NO3
     &       +1.033*R(95)   !ISPD
     &       + .80*R(96)    !ISO+NO2
C
      RO2new_isop =
     &        + 0.25*R(75) !ISOP+O
     &        + 0.20*R(77) !ISOP+O3
     &        + .114*R(93) !O3+ISPD =.154*HO2+.268*OH+.114*C2O3
     &        + .075*R(94) !ISPD+NO3
     &        + .967*R(95) !ISPD
C
C
C
C
C...  New OH from H2O2, HNO3, HONO, HNO4,  O3+HC, etc
      nn = nn + 1
      ptname(nn)  = 'newOHother'
C
      PA(nn) = cyc_HONO_l  
     &       + cyc_H2O2_l
     &       +      R(40)  !O+HCHO =OH+HO2+CO
     &       +      R(42)  !O+ALD  =C2O3+OH
     &       + 0.20*R(56)  !O+OLE  =.38*HO2+.2*OH
     &       + 0.10*R(58)  !O3+OLE =.44*HO2+.10*OH
     &       + 0.30*R(60)  !O+ETH  =HCHO+.70*XO2+CO+1.70*HO2+.30*OH
     &       + 0.08*R(71)  !OPEN+O3=.62*C2O3+.08*OH
     &       + OHnew_isop
C
C
C
C...  New HO2 from HCHO
      nn = nn + 1
      ptname(nn)  = 'nwHO2_HCHO'
C
      PA(nn) = 2.*R(38) !HCHO=2*HO2+CO
     &       +    R(40) !HCHO+O
     &       +    R(41) !HCHO+NO3
C
C
C
C...  New HO2 Production (total)
      nn = nn + 1
      ptname(nn)  = 'newHO2tot'
C
      HO2new = cyc_HNO4_l
     &       +2.00*R(38)    !HCHO=2.00*HO2
     &       +     R(40)    !HCHO+O=OH+HO2
     &       +     R(41)    !HCHO+NO3=HNO3+HO2
     &       +2.00*R(45)    !ALD=XO2+2.00*HO2
     &       + .38*R(56)    !O+OLE
     &       + .44*R(58)    !O3+OLE
     &       + 1.7*R(60)    !O+ETH
     &       + .12*R(62)    !O3+ETH
     &       +     R(69)    !OPEN=C2O3+HO2
     &       +2.00*R(70)    !OPEN+OH=2.00*HO2+C2O3 {count the HO2 as new and the C2O3 as prop}
     &       + .76*R(71)    !OPEN+O3
     &       +     R(74)    !MGLY=C2O3+HO2
     &       +  HO2new_isop
      PA(nn) = HO2new
C
C
C
C...  New RO2 Production
      nn = nn + 1
      ptname(nn)  = 'newRO2tot'
C
      RO2new = cyc_PAN_l
     &       +      R(42) !ALD+O
     &       +      R(44) !ALD+NO3
     &       +      R(69) !OPEN  
     &       + 0.62*R(71) !OPEN+O3
     &       +      R(74) !MGLY
     &       +      RO2new_isop
      PA(nn) = RO2new 
C
C
C
C...  New HOx from isoprene reactions
      nn = nn + 1
      ptname(nn)  = 'nHOx_isop'
C
      PA(nn) = OHnew_isop + HO2new_isop + RO2new_isop
C
C
C
C
C        {    R a d i c a l     P r o p a g a t i o n    }
C
C
C
C...  OH+CO, OH+CH4
      nn = nn + 1
      ptname(nn)  = 'OHwCO_CH4'
C
      PA(nn) =  R(36) !OH+CO
     &       +  R(51) !OH+CH4
C
C
C...  OH+ISO
      nn = nn + 1
      ptname(nn)  = 'ISOPwOH'
C
      PA(nn) = R(76)  !ISOP+OH

C
C...  other reactions of isoprene with Ox
      nn = nn + 1
      ptname(nn)  = 'ISOPwOx'
C
      PA(nn) = R(75) !ISOP+O3P
     &       + R(77) !ISOP+O3
     &       + R(78) !ISOP+NO3
     &       + R(96) !ISOP+NO2
C
C
C...  OH+HC
      nn = nn + 1
      ptname(nn)  = 'OHw_all_HC'
C
      OHwHC  = R(36)   !OH+CO
     &       + R(37)   !HCHO
     &       + R(43)   !ALD
     &       + R(51)   !OH+CH4
     &       + R(52)   !PAR
     &       + R(57)   !OLE  
     &       + R(61)   !ETH
     &       + R(63)   !TOL
     &       + R(66)   !CRES
     &       + R(70)   !OPEN
     &       + R(72)   !XYL
     &       + R(73)   !MGLY
     &       + R(76)   !ISOP+OH
     &       + R(84)   !{MEOH+OH}
     &       + R(85)   !{ETOH+OH}
     &       + R(92)   !OH+ISPD 
      PA(nn) = OHwHC
C
C  {NOTE:  OH+ISPD is included here because we want OHwISOP  }
C  {       to represent only the amount of isoprene reacted. }
C
C
C
c...  other OH propagation rxns 
      nn = nn + 1
      ptname(nn)  = 'OHpropmisc'
C
      other_OH_prop = cyc_O3OH_HO2_p   !net OH+O3=HO2   
     &               + R(35)            !OH+H2O2=HO2     
     &               + R(82)            !SO2+OH=SULF+HO2 
      PA(nn) = other_OH_prop
C
C
C...  Total HO2 Production = new + prop
      nn = nn + 1
      ptname(nn)  = 'HO2TotProd'
C
C     { sum all rad prop rxns that produce HO2}
C
      PA(nn)  =      HO2new
     &        +      other_OH_prop
     &        +      R(36)  !OH+CO 
     &        +      R(37)  !OH+HCHO
     &        +      R(46)  !C2O3+NO
     &        + 2.00*R(49)  !C2O3+C2O3
     &        +      R(51)  !OH+CH4
     &        + 0.11*R(52)  !OH+PAR
     &        + 0.94*R(53)  !ROR
     &        +      R(54)  !ROR
     &        +      R(57)  !OH+OLE
     &        +      R(61)  !OH+ETH
     &        + 0.44*R(63)  !OH+TOL
     &        + 0.90*R(64)  !TO2+NO
     &        +      R(65)  !TO2
     &        + 0.60*R(66)  !OH+CRES
     &        + 0.70*R(72)  !OH+XYL{don't count R70 because they are both new}
     &        +0.912*R(76)  !OH+ISO
     &        +      R(84)  !MEOH+OH
     &        +      R(85)  !ETOH+OH
     &        +0.503*R(92)  !OH+ISPD
C
C
C...  Total RO2 Prod = new + prop
      nn = nn + 1
      ptname(nn)  = 'RO2TotProd'
C
      PA(nn) =      RO2new !RO2new is summed above
     &       +      R(43)  !OH+ALD=C2O3
     &       + 0.56*R(63)  !OH+TOL=.56*TO2
     &       +      R(70)  !OH+OPEN=C2O3  + 2 new HO2
     &       + 0.30*R(72)  !OH+XYL=.30*TO2
     &       +      R(73)  !OH+MGLY=C2O3
     &       +0.498*R(92)  !OH+ISPD=.498*C2O3
C
C  Note: We won't count the PAR+OH=ROR as RO2 prod because
C        the ROR does not react with NO.
C
C
C
C...  NO2 from HO2
      nn = nn + 1
      ptname(nn)  = 'HO2_to_NO2'
C
      PA(nn) = R(28)  !HO2+NO=OH+NO2
C
C
C  Here we add up all rxns that convert HO2 to OH, including  
C  production of H2O2 followed by photolysis to produce 2 OH.
C  We exclude rxn of transported or initial H2O2 because it is
C  considered new OH, i.e., it is not propagation of HO2.
C  cyc_H2O2_l will be zero if there is net accum of H2O2.
C  If cyc_H2O2_l positive it represents the amount of 
C  H2O2 converted to OH in excess of the the amount of H2O2 produced.
C
C...  OH from HO2
      nn = nn + 1
      ptname(nn)  = 'HO2_to_OH'
C
      tmpprop = 2.*R(34) - cyc_H2O2_l
      if (tmpprop.LT.0.) tmpprop = 0.0

      PA(nn) = tmpprop        !H2O2=2 OH
     &       + cyc_O3OH_HO2_l !HO2+O3  
     &       +      R(28)     !HO2+NO
     &       + 0.79*R(50)     !HO2+C2O3

C
C
C...  NO2 from RO2
      nn = nn + 1
      ptname(nn)  = 'RO2_to_NO2'
C
      PA(nn) =     R(46)  !NO+C2O3
     &       + 0.9*R(64)  !NO+TO2
     &       +     R(79)  !NO+XO2
C
C
C
C          {    R a d i c a l    T e r m i n a t i o n    }
C
C NOTE: XO2N and CRO could be considered a radical but they always
C       terminate, so count them as termination.  GST
C
C
      OHterm  = cyc_HONO_p
     &        +      R(24)  !OH+HONO=NO2
     &        +      R(26)  !OH+NO2=HNO3
     &        +      R(27)  !OH+HNO3=NO3
     &        +      R(31)  !OH+HNO4=NO2
     &        + 0.13*R(52)  !PAR+OH  {ROR carrys a radical}
     &        + 0.40*R(66)  !OH+CRES
     &        +0.088*R(76)  !OH+ISO
     &        +      R(90)  !OH+HO2
C

C
C
C...  total OH prod
      nn = nn + 1
      ptname(nn)  = 'OH_reacted'
C
      PA(nn) =    OHterm + other_OH_prop + OHwHC 

C
C
C...  OHterm
      nn = nn + 1
      ptname(nn)  = 'OHterm'
C
      PA(nn) = OHterm
C
C
C
C...  HO2 termination rxns 
      nn = nn + 1
      ptname(nn)  = 'HO2term'
C
      PA(nn) = cyc_HNO4_p 
     &       + cyc_H2O2_p
     &       + 0.21*R(50)  !HO2+C2O3
     &       +      R(86)  !XO2+HO2=
     &       +      R(87)  !XO2N+HO2=
     &       +      R(90)  !OH+HO2

C
C
C     { organic rads are: C2O3, ROR, TO2 }
C...  RO2term
      nn = nn + 1
      ptname(nn)  = 'RO2term'
C   
      PA(nn) = cyc_PAN_p
     &        + 0.21*R(50)  !HO2+C2O3
     &        + 0.06*R(53)  !ROR=.94*HO2
     &        +      R(55)  !ROR+NO2=ORNIT
     &        + 0.10*R(64)  !TO2+NO=.90*HO2

C
C...  HCHO Production from isoprene
      nn = nn + 1
      ptname(nn)  = 'HCHOp_isop'
C
      prod_hcho_from_isop =
     &       + 0.50*R(75)  !ISOP+O
     &       +0.629*R(76)  !ISOP+OH
     &       + 0.60*R(77)  !ISOP+O3
     &       +0.167*R(92)  !ISPD+OH
     &       + 0.15*R(93)  !ISPD+O3
     &       +0.282*R(94)  !ISOP+NO3
     &       + 0.90*R(95)  !ISOP
      PA(nn) = prod_hcho_from_isop     
C
C
C
C...  Total HCHO Production
      nn = nn + 1
      ptname(nn)  = 'HCHOp_Tot'
C
      PA(nn) =      R(45)  !ALD=
     &       +      R(46)  !C2O3+NO
     &       + 2.00*R(49)  !C2O3+C2O3
     &       + 0.79*R(50)  !C2O3+HO2
     &       +      R(51)  !CH4+OH
     &       + 0.20*R(56)  !O+OLE
     &       +      R(57)  !OH+OLE
     &       + 0.74*R(58)  !O3+OLE
     &       +      R(59)  !NO3+OLE
     &       +      R(60)  !O+ETH
     &       + 1.56*R(61)  !OH+ETH
     &       +      R(62)  !O3+ETH
     &       +      R(70)  !OPEN+OH
     &       + 0.70*R(71)  !OPEN+O3
     &       +      R(84)  !MEOH+OH
     &       + prod_hcho_from_isop
C
C
C
C
C            {  N O z     C h e m i s t r y      }
C
C
C...  OH+NO2=HNO3 
      nn = nn + 1
      ptname(nn)  = 'HNO3_OHNO2'
C
      PA(nn) = R(26)  !OH+NO2
C
C
C
C...  NO3+HC=HNO3
      nn = nn + 1
      ptname(nn)  = 'HNO3_NO3HC'
C
      PA(nn) =      R(41)  !NO3+HCHO
     &       +      R(44)  !NO3+ALD
     &       +      R(67)  !NO3+CRES
     &       + 0.15*R(94)  !NO3+ISPD
C
C
C
C...  other HNO3
      nn = nn + 1
      ptname(nn)  = 'HNO3_N2O5'
C
      PA(nn) = 2.0*R(18)  !N2O5+H2O
C
C
C
C...  HNO3 reacted
      nn = nn + 1
      ptname(nn)  = 'HNO3reacte'
C
      PA(nn) = R(27)
C
C
C
C...  net PAN prod
      nn = nn + 1
      ptname(nn)  = 'PANprodNet'
C
      PA(nn) = cyc_PAN_p
C
C
C
C...  net PAN loss
      nn = nn + 1
      ptname(nn)  = 'PANlossNet'
C
      PA(nn) = cyc_PAN_l
C
C
C
C... production of ONIT
      nn = nn + 1
      ptname(nn)  = 'RNO3_prod'
C
      PA(nn) =      R(55)  !ROR+NO2
     &       + 0.10*R(64)  !TO2+NO
     &       +      R(68)  !CRO+NO2
     &       + 0.80*R(78)  !NO3+ISO
     &       +      R(81)  !XO2N+NO
     &       + 0.85*R(94)  !NO3+ISPD
     &       + 0.80*R(96)  !NO2+ISPD
C
C
C   { convert to ppb/hr}
      Do n = 1,nn
        PA(n) = 1000.*PA(n)
      End do
C
C
      If( nn .GT. MXCPA ) Then
         write(iout,'(//,a)') 'ERROR in CPAMECH3:'
         write(iout,*) 'Number of outputs requested exceeds limit.'
         write(iout,*) 'Increase parameter MXCPA and try again.'
         call camxerr()
      End if

      NPA_INIT = nn  !set num of outputs on first dummy call from pasetup

      Return
      End
