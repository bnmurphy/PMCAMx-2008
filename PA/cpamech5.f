      subroutine cpamech5( r, nr, pa, npa, npa_init )
c
c-----CAMx v4.02 030709
c
c     Copyright 2002, 2003
c     ENVIRON International Corporation
c
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
C**   01     12/02  GST  Created for use with the SAPRC99 photochemical mechanism.
C**                      This is hardwired for the CAMx SAPRC99 vers 3.1
C**                      
C**                         
C***********************************************************************
C
C
      include 'camx.prm'
      include 'filunit.com'
      include 'tracer.com'
      include 'procan.com'
C
C... Local variable declarations
C
      INTEGER  NR, NPA, NPA_INIT, nn
      REAL     R(NR), PA(NPA)

      Real     cyc_HONO_p, cyc_HONO_l, cyc_H2O2_p, cyc_H2O2_l
      Real     cyc_ROOH_p, cyc_ROOH_l, cyc_COOH_p, cyc_COOH_l
      Real     cyc_HNO4_p, cyc_HNO4_l

      Real     cyc_PAN_p,  cyc_PAN_l,  cyc_PAN2_p,  cyc_PAN2_l
      Real     cyc_PBZN_p, cyc_PBZN_l, cyc_PPN_p, cyc_PPN_l
CC
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
cc

c
c... { Ox Production }
      nn =  1
      ptname(nn)  = 'OxProd'
      PA(nn) = 
     &      + 2.000*R( 10) !{NO+NO}
     &      + 1.000*R( 23) !{HONO}
     &      + 1.000*R( 24) !{OH+HONO}
     &      + 1.000*R( 27) !{OH+HNO3}
     &      + 1.000*R( 31) !{HO2+NO}
     &      + 0.390*R( 34) !{HNO4}
     &      + 1.000*R( 46) !{CXO2+NO}
     &      + 1.000*R( 51) !{RO2R+NO}
     &      + 1.000*R( 56) !{R2O2+NO}
     &      + 1.000*R( 62) !{RO2N+NO}
     &      + 1.000*R( 71) !{CCO3+NO}
     &      + 0.250*R( 72) !{CCO3+HO2}
     &      + 1.000*R( 81) !{RCO3+NO}
     &      + 0.250*R( 82) !{RCO3+HO2}
     &      + 1.000*R( 92) !{BZCO+NO}
     &      + 0.250*R( 93) !{BZCO+HO2}
     &      + 1.000*R(104) !{MCO3+NO}
     &      + 0.250*R(105) !{MCO3+HO2}
     &      + 1.000*R(120) !{BZNO+NO2}
     &      + 1.000*R(128) !{HCO3+NO}
c
c
c... { Ox Loss }
      nn =  nn + 1
      ptname(nn)  = 'OxLoss'
      PA(nn) = 
     &      + 2.000*R(  3) !{O+O3}
     &      + 2.000*R(  5) !{O+NO2}
     &      + 1.000*R( 13) !{N2O5+H2O}
     &      + 2.000*R( 14) !{NO2+NO3}
     &      + 2.000*R( 15) !{NO3}
     &      + 1.000*R( 19) !{O1D+H2O}
     &      + 1.000*R( 26) !{OH+NO3}
     &      + 1.000*R( 30) !{OH+O3}
     &      + 1.000*R( 36) !{HO2+O3}
     &      + 1.000*R( 39) !{NO3+HO2}
     &      + 2.000*R( 40) !{NO3+NO3}
     &      + 1.000*R( 48) !{CXO2+NO3}
     &      + 1.000*R( 53) !{RO2R+NO3}
     &      + 1.000*R( 58) !{R2O2+NO3}
     &      + 1.000*R( 65) !{RO2N+NO3}
     &      + 1.000*R( 73) !{CCO3+NO3}
     &      + 1.000*R( 83) !{RCO3+NO3}
     &      + 1.000*R( 94) !{BZCO+NO3}
     &      + 1.000*R(106) !{MCO3+NO3}
     &      + 1.000*R(117) !{BZO+NO2}
     &      + 1.000*R(129) !{HCHO+NO3}
     &      + 1.000*R(132) !{CCHO+NO3}
     &      + 1.000*R(135) !{RCHO+NO3}
     &      + 1.000*R(148) !{GLY+NO3}
     &      + 1.000*R(151) !{MGLY+NO3}
     &      + 1.000*R(154) !{PHEN+NO3}
     &      + 1.000*R(156) !{CRES+NO3}
     &      + 1.000*R(157) !{NPHE+NO3}
     &      + 1.000*R(160) !{BALD+NO3}
     &      + 1.000*R(162) !{METH+O3}
     &      + 1.000*R(163) !{METH+NO3}
     &      + 1.000*R(164) !{METH+O}
     &      + 1.000*R(167) !{MVK+O3}
     &      + 1.000*R(168) !{MVK+O}
     &      + 1.000*R(171) !{ISPD+O3}
     &      + 1.000*R(172) !{ISPD+NO3}
     &      + 1.000*R(179) !{DCB1+O3}
     &      + 1.000*R(186) !{ETHE+O3}
     &      + 1.000*R(187) !{ETHE+NO3}
     &      + 1.000*R(188) !{ETHE+O}
     &      + 1.000*R(190) !{ISOP+O3}
     &      + 1.000*R(191) !{ISOP+NO3}
     &      + 1.000*R(192) !{ISOP+O}
     &      + 1.000*R(194) !{TERP+O3}
     &      + 1.000*R(195) !{TERP+NO3}
     &      + 1.000*R(196) !{TERP+O}
     &      + 1.000*R(205) !{OLE1+O3}
     &      + 1.000*R(206) !{OLE1+NO3}
     &      + 1.000*R(207) !{OLE1+O}
     &      + 1.000*R(209) !{OLE2+O3}
     &      + 1.000*R(210) !{OLE2+NO3}
     &      + 1.000*R(211) !{OLE2+O}
c
C
C
C  { Radical Reservoirs: HONO, HNO4, H2O2  etc                              }
C  { we will call these termination rxns if there is net production         }
C  { during a time step, or initiation if there is net release of HOx.      }
C  { It may be useful to distinguish this initiation from NO3, O3 + HC rxns }
C
C  { 21} NO+OH=HONO
C  { 22} HONO=OH+NO
C  { 23} HONO=0.100*HO2+0.100*NO2
C
        sum =  R(21) - R(20) - 0.1*R(23)  ! HONO cycle
        if (sum.GT.0.) then
          cyc_HONO_p = sum
          cyc_HONO_l = 0.0
         else
          cyc_HONO_p =  0.0
          cyc_HONO_l = -sum
         end if
C
C  { 32} HO2+NO2=HNO4
C  { 33} HNO4=HO2+NO2
C  { 34} HNO4=0.610*HO2+0.610*NO2+0.390*OH+0.390*NO3
        sum =  R(32) - R(33)  - R(34)    ! HNO4 cycle
        if (sum.GT.0.) then
          cyc_HNO4_p = sum
          cyc_HNO4_l = 0.0
         else
          cyc_HNO4_p =  0.0
          cyc_HNO4_l = -sum
         end if
C
C
C{ 37} HO2+HO2=H2O2
C{ 38} HO2+HO2+M=H2O2
C{ 41} H2O2=2.00*OH
C
        sum =  2.*(  R(37) + R(38) - R(41)  )   ! H2O2 cycle
        if (sum.GT.0.) then
          cyc_H2O2_p = sum
          cyc_H2O2_l = 0.0
         else
          cyc_H2O2_p =  0.0
          cyc_H2O2_l = -sum
         end if
C
C
C{ ..} RO2+HO2=OOH
C{ 69} OOH=HO2+OH
C
        sum =  R( 52)   ! ROOH cycle
     &      +  R( 63)   ! ROOH cycle
     &      -  R(144)   ! ROOH cycle
     &      +  R( 47)   ! COOH cycle
     &      -  R(142)   ! COOH cycle
        if (sum.GT.0.) then
          cyc_ROOH_p = sum
          cyc_ROOH_l = 0.0
         else
          cyc_ROOH_p =  0.0
          cyc_ROOH_l = -sum
         end if
C
C
C  { 69} CCO3+NO2=PAN
C  { 70} PAN=CCO3+NO2
C
        sum = R(69) - R(70)    ! PAN Cycle
        if (sum.GT.0.) then
          cyc_PAN_p = sum
          cyc_PAN_l = 0.0
        else
          cyc_PAN_p =  0.0
          cyc_PAN_l = -sum
        end if
C
C
C  { 79} RCO3+NO2=PAN2
C  { 80} PAN2=RCO3+NO2
C
        sum = R(79) - R(80)    ! PAN2 Cycle
        if (sum.GT.0.) then
          cyc_PAN2_p = sum
          cyc_PAN2_l = 0.0
        else
          cyc_PAN2_p =  0.0
          cyc_PAN2_l = -sum
        end if
C
C  { 90} BZCO_O2+NO2=PBZN
C  { 91} PBZN=BZCO3+NO2                   ! PBZN Cycle
        sum = R(90) - R(91)    
        if (sum.GT.0.) then
          cyc_PBZN_p = sum
          cyc_PBZN_l = 0.0
        else
          cyc_PBZN_p =  0.0
          cyc_PBZN_l = -sum
        end if
C
C  {102} MCO3+NO2=MPAN
C  {103} MPAN=MCO3+NO2
C
        sum =  R(102) - R(103)   ! MPAN Cycle
        if (sum.GT.0.) then
          cyc_MPAN_p = sum
          cyc_MPAN_l = 0.0
        else
          cyc_MPAN_p =  0.0
          cyc_MPAN_l = -sum
        end if
C
C
C
c
c  {   R a d i c a l     I n i t i a t i o n   }
c
C
C
C...  New OH from O1D+H2O
      nn = nn + 1
      ptname(nn)  = 'newOH_O1D'
C
      PA(nn) = 2.*R(19)   !O1D+H2O->2*OH
C
c
c...{ Rad initiation for OH       }
      nn =  nn + 1
      ptname(nn)  = 'newOHother'
      PA(nn) = cyc_HONO_l  
     &       + cyc_H2O2_l
     &       + cyc_ROOH_l             !also some new OH from HNO4 but small so count as HO2
     &      + 1.000*R( 28) !{HNO3}
     &      + 0.390*R( 34) !{HNO4}
     &      + 0.208*R(162) !{METH+O3}
     &      + 0.330*R(165) !{METH}
     &      + 0.164*R(167) !{MVK+O3}
     &      + 0.285*R(171) !{ISPD+O3}
     &      + 0.500*R(179) !{DCB1+O3}
     &      + 0.120*R(186) !{ETHE+O3}
     &      + 0.266*R(190) !{ISOP+O3}
     &      + 0.567*R(194) !{TERP+O3}
     &      + 0.155*R(205) !{OLE1+O3}
     &      + 0.378*R(209) !{OLE2+O3}
C
C
C...  New HO2 from HCHO
      nn = nn + 1
      ptname(nn)  = 'nwHO2_HCHO'
C
      PA(nn) = 2.*R(123) !HCHO=2*HO2+CO
     &       +    R(129) !HCHO+NO3
C
C
Cc
c
c...{ Rad initiation for HO2      }
      nn =  nn + 1
      ptname(nn)  = 'newHO2tot '
      HO2new = cyc_ROOH_l 
     &      +  cyc_HNO4_l !count the new OH from HNO4 here
     &      + 2.000*R(123) !{HCHO}
     &      + 1.000*R(127) !{HCO3}
     &      + 1.000*R(128) !{HCO3+NO}
     &      + 1.000*R(129) !{HCHO+NO3}
     &      + 1.000*R(131) !{CCHO}
     &      + 1.000*R(134) !{RCHO}
     &      + 2.000*R(145) !{GLY}
     &      + 0.630*R(148) !{GLY+NO3}
     &      + 1.000*R(149) !{MGLY}
     &      + 0.008*R(162) !{METH+O3}
     &      + 0.340*R(165) !{METH}
     &      + 0.064*R(167) !{MVK+O3}
     &      + 0.400*R(171) !{ISPD+O3}
     &      + 1.233*R(173) !{ISPD}
     &      + 0.341*R(177) !{RNO3}
     &      + 1.500*R(179) !{DCB1+O3}
     &      + 0.500*R(181) !{DCB2}
     &      + 0.500*R(183) !{DCB3}
     &      + 0.120*R(186) !{ETHE+O3}
     &      + 0.500*R(188) !{ETHE+O}
     &      + 0.033*R(194) !{TERP+O3}
     &      + 0.056*R(205) !{OLE1+O3}
     &      + 0.003*R(209) !{OLE2+O3}
     &      + 0.013*R(211) !{OLE2+O}
      PA(nn) = HO2new
c
c
c... { begin summing new RO2 }
c
c...{ Rad initiation for CXO2     }
c      nn =  nn + 1
c      ptname(nn)  = 'newCXO2    '
c      PA(nn) = 
      RO2new =
     &      + 1.000*R(131) !{CCHO}
     &      + 1.000*R(137) !{ACET}
     &      + 0.300*R(169) !{MVK}
     &      + 0.300*R(188) !{ETHE+O}
     &      + 0.250*R(192) !{ISOP+O}
     &      + 0.076*R(205) !{OLE1+O3}
     &      + 0.197*R(209) !{OLE2+O3}
     &      + 0.030*R(210) !{OLE2+NO3}
c
c
c...{ Rad initiation for CCO3     }
c      nn =  nn + 1
c      ptname(nn)  = 'newCCO3    '
c      PA(nn) = 
      RO2new = RO2new
     &      + 1.000*R( 70) !{PAN}
     &      + 1.000*R(132) !{CCHO+NO3}
     &      + 1.000*R(137) !{ACET}
     &      + 1.000*R(139) !{MEK}
     &      + 1.000*R(149) !{MGLY}
     &      + 1.000*R(151) !{MGLY+NO3}
     &      + 2.000*R(152) !{BACL}
     &      + 0.670*R(165) !{METH}
     &      + 0.467*R(173) !{ISPD}
     &      + 0.667*R(175) !{PROD}
     &      + 0.500*R(181) !{DCB2}
     &      + 0.500*R(183) !{DCB3}
     &      + 0.123*R(194) !{TERP+O3}
     &      + 0.137*R(209) !{OLE2+O3}
c
c
c...{ Rad initiation for RCO3     }
c      nn =  nn + 1
c      ptname(nn)  = 'newRCO3    '
c      PA(nn) = 
      RO2new = RO2new
     &      + 1.000*R( 80) !{PAN2}
     &      + 1.000*R(135) !{RCHO+NO3}
     &      + 0.370*R(148) !{GLY+NO3}
     &      + 0.100*R(162) !{METH+O3}
     &      + 0.050*R(167) !{MVK+O3}
     &      + 0.048*R(171) !{ISPD+O3}
     &      + 0.300*R(173) !{ISPD}
     &      + 0.333*R(175) !{PROD}
     &      + 0.201*R(194) !{TERP+O3}
     &      + 0.006*R(209) !{OLE2+O3}
c
c
c...{ Rad initiation for RO2R     }
c      nn =  nn + 1
c      ptname(nn)  = 'newRO2R    '
c      PA(nn) = 
      RO2new = RO2new
     &      + 1.000*R(134) !{RCHO}
     &      + 1.000*R(139) !{MEK}
     &      + 0.100*R(162) !{METH+O3}
     &      + 0.500*R(163) !{METH+NO3}
     &      + 0.330*R(165) !{METH}
     &      + 0.050*R(167) !{MVK+O3}
     &      + 0.048*R(171) !{ISPD+O3}
     &      + 0.799*R(172) !{ISPD+NO3}
     &      + 0.960*R(175) !{PROD}
     &      + 0.564*R(177) !{RNO3}
     &      + 1.000*R(181) !{DCB2}
     &      + 1.000*R(183) !{DCB3}
     &      + 1.000*R(187) !{ETHE+NO3}
     &      + 0.200*R(188) !{ETHE+O}
     &      + 0.066*R(190) !{ISOP+O3}
     &      + 0.749*R(191) !{ISOP+NO3}
     &      + 0.031*R(194) !{TERP+O3}
     &      + 0.276*R(195) !{TERP+NO3}
     &      + 0.022*R(205) !{OLE1+O3}
     &      + 0.824*R(206) !{OLE1+NO3}
     &      + 0.033*R(209) !{OLE2+O3}
     &      + 0.442*R(210) !{OLE2+NO3}
     &      + 0.012*R(211) !{OLE2+O}
c
c
c...{ Rad initiation for RO2N     }
c      nn =  nn + 1
c      ptname(nn)  = 'newRO2N    '
c      PA(nn) = 
      RO2new = RO2new
     &      + 0.051*R(172) !{ISPD+NO3}
     &      + 0.040*R(175) !{PROD}
     &      + 0.095*R(177) !{RNO3}
     &      + 0.008*R(190) !{ISOP+O3}
     &      + 0.064*R(191) !{ISOP+NO3}
     &      + 0.010*R(192) !{ISOP+O}
     &      + 0.180*R(194) !{TERP+O3}
     &      + 0.250*R(195) !{TERP+NO3}
     &      + 0.001*R(205) !{OLE1+O3}
     &      + 0.176*R(206) !{OLE1+NO3}
     &      + 0.002*R(209) !{OLE2+O3}
     &      + 0.136*R(210) !{OLE2+NO3}
     &      + 0.001*R(211) !{OLE2+O}
c
c
c...{ Rad initiation for BZCO     }
c      nn =  nn + 1
c      ptname(nn)  = 'newBZCO    '
c      PA(nn) = 
      RO2new = RO2new
     &      + 1.000*R( 91) !{PBZN}
     &      + 1.000*R(160) !{BALD+NO3}
c
c
c...{ Rad initiation for MCO3     }
c      nn =  nn + 1
c      ptname(nn)  = 'newMCO3    '
c      PA(nn) = 
      RO2new = RO2new
     &      + 1.000*R(103) !{MPAN}
     &      + 0.500*R(163) !{METH+NO3}
     &      + 0.330*R(165) !{METH}
     &      + 0.300*R(169) !{MVK}
     &      + 0.150*R(172) !{ISPD+NO3}
     &      + 0.192*R(190) !{ISOP+O3}
     &      + 0.240*R(192) !{ISOP+O}
c
c
c...{ Rad initiation for TBUO is zero  }
c
c...{ Total RO2 Rad initiation  }
c
      nn =  nn + 1
      ptname(nn)  = 'newRO2tot'
      PA(nn) =   RO2new 
c
c
c

C...  OH+CO, OH+CH4
      nn = nn + 1
      ptname(nn)  = 'OHwCO_CH4'
C
      PA(nn) =  R( 29) !OH+CO
     &       +  R(184) !OH+CH4
C
C
C...  OH+ISO
      nn = nn + 1
      ptname(nn)  = 'ISOPwOH'
C
      PA(nn) = R(189)  !ISOP+OH

C
C...  other reactions of isoprene with Ox
      nn = nn + 1
      ptname(nn)  = 'ISOPwOx'
C
      PA(nn) = 
     &       + R(190) !ISOP+O3
     &       + R(191) !ISOP+NO3
     &       + R(192) !ISOP+O
C
Cc
c
c... { OH rxns w/Organics }
      nn =  nn + 1
      ptname(nn)  = 'OHwHC'
      OHwHC = 
     &      +R( 29) !{OH+CO}
     &      +R(125) !{HCHO+OH}
     &      +R(130) !{CCHO+OH}
     &      +R(133) !{RCHO+OH}
     &      +R(136) !{ACET+OH}
     &      +R(138) !{MEK+OH}
     &      +R(140) !{MEOH+OH}
     &      +R(141) !{COOH+OH}
     &      +R(143) !{ROOH+OH}
     &      +R(147) !{GLY+OH}
     &      +R(150) !{MGLY+OH}
     &      +R(153) !{PHEN+OH}
     &      +R(155) !{CRES+OH}
     &      +R(158) !{BALD+OH}
     &      +R(161) !{METH+OH}
     &      +R(166) !{MVK+OH}
     &      +R(170) !{ISPD+OH}
     &      +R(174) !{PROD+OH}
     &      +R(176) !{RNO3+OH}
     &      +R(178) !{DCB1+OH}
     &      +R(180) !{DCB2+OH}
     &      +R(182) !{DCB3+OH}
     &      +R(184) !{CH4+OH}
     &      +R(185) !{ETHE+OH}
     &      +R(189) !{ISOP+OH}
     &      +R(193) !{TERP+OH}
     &      +R(197) !{ALK1+OH}
     &      +R(198) !{ALK2+OH}
     &      +R(199) !{ALK3+OH}
     &      +R(200) !{ALK4+OH}
     &      +R(201) !{ALK5+OH}
     &      +R(202) !{ARO1+OH}
     &      +R(203) !{ARO2+OH}
     &      +R(204) !{OLE1+OH}
     &      +R(208) !{OLE2+OH}
      PA(nn) = OHwHC
C
C  {NOTE:  OH+ISPD is included here because we want OHwISOP  }
C  {       to represent only the amount of isoprene reacted. }
C
C
C
c  {other OH propagation rxns }
c
c
c... { other OH prop rxns }
      nn =  nn + 1
      ptname(nn)  = 'OHpropmisc'

      other_OH_prop = 
     &      +R( 26) !{OH+NO3}
     &      +R( 30) !{OH+O3}  !need to reconsisder whether to treat this as a cycle...gst
     &      +R( 42) !{HO2H+OH}
     &      +R( 44) !{OH+SO2}
     &      +R( 45) !{OH+H2}
      PA(nn) = other_OH_prop
c
c
c
c
c... { HCHO reactions }
c      nn =  nn + 1
c      ptname(nn)  = 'HCHOreacted'
c      PA(nn) = 
c     &      +R(123) !{HCHO}
c     &      +R(124) !{HCHO}
c     &      +R(125) !{HCHO+OH}
c     &      +R(126) !{HCHO+HO2}
c     &      +R(129) !{HCHO+NO3}
c
C
c
c
c... { total HO2 Production }
      nn =  nn + 1
      ptname(nn)  = 'HO2prod'
      PA(nn) = HO2new      !HO2new is summed above
     &      + 1.000*R( 26) !{OH+NO3}
     &      + 1.000*R( 29) !{OH+CO}
     &      + 1.000*R( 30) !{OH+O3}
     &      + 1.000*R( 42) !{HO2H+OH}
     &      + 1.000*R( 44) !{OH+SO2}
     &      + 1.000*R( 45) !{OH+H2}
     &      + 1.000*R( 46) !{CXO2+NO}
     &      + 1.000*R( 48) !{CXO2+NO3}
     &      + 2.000*R( 50) !{CXO2+CXO2}
     &      + 1.000*R( 51) !{RO2R+NO}
     &      + 1.000*R( 53) !{RO2R+NO3}
     &      + 1.000*R( 54) !{RO2R+CXO2}
     &      + 1.000*R( 55) !{RO2R+RO2R}
     &      + 1.000*R( 57) !{R2O2+HO2}
     &      + 1.000*R( 64) !{RO2N+CXO2}
     &      + 1.000*R( 65) !{RO2N+NO3}
     &      + 1.000*R( 66) !{RO2N+RO2R}
     &      + 1.000*R( 68) !{RO2N+RO2N}
     &      + 1.000*R(125) !{HCHO+OH}
     &      + 1.000*R(140) !{MEOH+OH}
     &      + 0.630*R(147) !{GLY+OH}
     &      + 0.379*R(174) !{PROD+OH}
     &      + 0.113*R(176) !{RNO3+OH}
     &      + 0.121*R(198) !{ALK2+OH}
     &      + 0.224*R(202) !{ARO1+OH}
     &      + 0.187*R(203) !{ARO2+OH}
c
c
c... { total RO2 Production }
      nn =  nn + 1
      ptname(nn)  = 'RO2prod'
      PA(nn) =  RO2new     !RO2new is summed above
     &      + 1.000*R( 59) !{R2O2+CXO2}
     &      + 1.000*R( 60) !{R2O2+RO2R}
     &      + 1.000*R( 67) !{RO2N+R2O2}
     &      + 1.000*R( 71) !{CCO3+NO}
     &      + 1.000*R( 73) !{CCO3+NO3}
     &      + 1.000*R( 76) !{CCO3+R2O2}
     &      + 2.000*R( 78) !{CCO3+CCO3}
     &      + 1.000*R( 81) !{RCO3+NO}
     &      + 1.000*R( 83) !{RCO3+NO3}
     &      + 1.000*R( 86) !{RCO3+R2O2}
     &      + 2.000*R( 88) !{RCO3+CCO3}
     &      + 2.000*R( 89) !{RCO3+RCO3}
     &      + 1.000*R( 97) !{BZCO+R2O2}
     &      + 1.000*R( 99) !{BZCO+CCO3}
     &      + 1.000*R(100) !{BZCO+RCO3}
     &      + 1.000*R(104) !{MCO3+NO}
     &      + 1.000*R(106) !{MCO3+NO3}
     &      + 1.000*R(109) !{MCO3+R2O2}
     &      + 2.000*R(111) !{MCO3+CCO3}
     &      + 2.000*R(112) !{MCO3+RCO3}
     &      + 1.000*R(113) !{MCO3+BZCO}
     &      + 2.000*R(114) !{MCO3+MCO3}
     &      + 1.000*R(116) !{TBUO}
     &      + 1.000*R(130) !{CCHO+OH}
     &      + 1.000*R(133) !{RCHO+OH}
     &      + 1.000*R(136) !{ACET+OH}
     &      + 1.000*R(138) !{MEK+OH}
     &      + 0.650*R(141) !{COOH+OH}
     &      + 0.340*R(143) !{ROOH+OH}
     &      + 0.370*R(147) !{GLY+OH}
     &      + 1.000*R(150) !{MGLY+OH}
     &      + 0.760*R(153) !{PHEN+OH}
     &      + 0.760*R(155) !{CRES+OH}
     &      + 1.000*R(158) !{BALD+OH}
     &      + 1.000*R(161) !{METH+OH}
     &      + 1.000*R(166) !{MVK+OH}
     &      + 1.000*R(170) !{ISPD+OH}
     &      + 0.621*R(174) !{PROD+OH}
     &      + 0.549*R(176) !{RNO3+OH}
     &      + 1.000*R(178) !{DCB1+OH}
     &      + 1.000*R(180) !{DCB2+OH}
     &      + 1.000*R(182) !{DCB3+OH}
     &      + 1.000*R(184) !{CH4+OH}
     &      + 1.000*R(185) !{ETHE+OH}
     &      + 1.000*R(189) !{ISOP+OH}
     &      + 1.000*R(193) !{TERP+OH}
     &      + 1.000*R(197) !{ALK1+OH}
     &      + 0.633*R(198) !{ALK2+OH}
     &      + 1.000*R(199) !{ALK3+OH}
     &      + 1.000*R(200) !{ALK4+OH}
     &      + 1.000*R(201) !{ALK5+OH}
     &      + 0.776*R(202) !{ARO1+OH}
     &      + 0.813*R(203) !{ARO2+OH}
     &      + 1.000*R(204) !{OLE1+OH}
     &      + 1.000*R(208) !{OLE2+OH}
c
c
c
c
c... { OH termination rxns }
      OHterm = cyc_HONO_p
     &      + 1.000*R( 24) !{OH+HONO}
     &      + 1.000*R( 25) !{OH+NO2}
     &      + 1.000*R( 27) !{OH+HNO3}
     &      + 1.000*R( 35) !{HNO4+OH}
     &      + 1.000*R( 43) !{OH+HO2}
     &      + 0.240*R(153) !{PHEN+OH}
     &      + 0.240*R(155) !{CRES+OH}
     &      + 0.338*R(176) !{RNO3+OH}
c
C
C
C...  total OH prod
      nn = nn + 1
      ptname(nn)  = 'OH_reacted'
C
      PA(nn) =    OHterm + other_OH_prop + OHwHC 
C
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
C... { HO2 termination rxns }
      nn =  nn + 1
      ptname(nn)  = 'HO2term'
      PA(nn) = 
     &      + 1.000*R( 32) !{HO2+NO2}
     &      + 2.000*R( 37) !{HO2+HO2}
     &      + 2.000*R( 38) !{HO2+HO2}
     &      + 0.200*R( 39) !{NO3+HO2}
     &      + 1.000*R( 43) !{OH+HO2}
     &      + 1.000*R( 47) !{CXO2+HO2}
     &      + 1.000*R( 52) !{RO2R+HO2}
     &      + 1.000*R( 63) !{RO2N+HO2}
     &      + 1.000*R( 72) !{CCO3+HO2}
     &      + 1.000*R( 82) !{RCO3+HO2}
     &      + 1.000*R( 93) !{BZCO+HO2}
     &      + 1.000*R(105) !{MCO3+HO2}
     &      + 1.000*R(118) !{BZO+HO2}
     &      + 1.000*R(121) !{BZNO+HO2}
     &      + 1.000*R(126) !{HCHO+HO2}
c
c
c... { RO2 termination rxns }
      nn =  nn + 1
      ptname(nn)  = 'RO2term'
      PA(nn) = 
     &      + 1.000*R( 47) !{CXO2+HO2}
     &      + 2.000*R( 49) !{CXO2+CXO2}
     &      + 1.000*R( 52) !{RO2R+HO2}
     &      + 1.000*R( 54) !{RO2R+CXO2}
     &      + 1.000*R( 55) !{RO2R+RO2R}
     &      + 1.000*R( 62) !{RO2N+NO}
     &      + 1.000*R( 63) !{RO2N+HO2}
     &      + 1.000*R( 64) !{RO2N+CXO2}
     &      + 1.000*R( 66) !{RO2N+RO2R}
     &      + 1.000*R( 68) !{RO2N+RO2N}
     &      + 1.000*R( 69) !{CCO3+NO2}
     &      + 1.000*R( 72) !{CCO3+HO2}
     &      + 2.000*R( 74) !{CCO3+CXO2}
     &      + 2.000*R( 75) !{CCO3+RO2R}
     &      + 2.000*R( 77) !{CCO3+RO2N}
     &      + 1.000*R( 79) !{RCO3+NO2}
     &      + 1.000*R( 82) !{RCO3+HO2}
     &      + 2.000*R( 84) !{RCO3+CXO2}
     &      + 2.000*R( 85) !{RCO3+RO2R}
     &      + 2.000*R( 87) !{RCO3+RO2N}
     &      + 1.000*R( 90) !{BZCO+NO2}
     &      + 1.000*R( 92) !{BZCO+NO}
     &      + 1.000*R( 93) !{BZCO+HO2}
     &      + 1.000*R( 94) !{BZCO+NO3}
     &      + 2.000*R( 95) !{BZCO+CXO2}
     &      + 2.000*R( 96) !{BZCO+RO2R}
     &      + 2.000*R( 98) !{BZCO+RO2N}
     &      + 1.000*R( 99) !{BZCO+CCO3}
     &      + 1.000*R(100) !{BZCO+RCO3}
     &      + 2.000*R(101) !{BZCO+BZCO}
     &      + 1.000*R(102) !{MCO3+NO2}
     &      + 1.000*R(105) !{MCO3+HO2}
     &      + 2.000*R(107) !{MCO3+CXO2}
     &      + 2.000*R(108) !{MCO3+RO2R}
     &      + 2.000*R(110) !{MCO3+RO2N}
     &      + 1.000*R(113) !{MCO3+BZCO}
     &      + 1.000*R(115) !{TBUO+NO2}



c  { NO2 from  HO2  }
c
c
cc... { NO2 from HO2 }
c      nn =  nn + 1
c      ptname(nn)  = 'NO2fromHO2'
c      PA(nn) = 
c     &      + 1.000*R( 31) !{HO2+NO}
c
c
cc... { NO2 from RO2 }
c      nn =  nn + 1
c      ptname(nn)  = 'NO2fromRO2'
c      PA(nn) = 
c     &      + 1.000*R( 46) !{CXO2+NO}
c     &      + 1.000*R( 51) !{RO2R+NO}
c     &      + 1.000*R( 56) !{R2O2+NO}
c     &      + 1.000*R( 71) !{CCO3+NO}
c     &      + 1.000*R( 81) !{RCO3+NO}
c     &      + 1.000*R( 92) !{BZCO+NO}
c     &      + 1.000*R(104) !{MCO3+NO}
c
C
C... { HCHO Production from isoprene }
      nn = nn + 1
      ptname(nn)  = 'HCHOp_isop'
C
      prod_hcho_from_isop =
     &      + 0.055*R(170) !{ISPD+OH}
     &      + 0.125*R(171) !{ISPD+O3}
     &      + 0.227*R(172) !{ISPD+NO3}
     &      + 0.300*R(173) !{ISPD}
     &      + 0.624*R(189) !{ISOP+OH}
     &      + 0.592*R(190) !{ISOP+O3}
     &      + 0.240*R(192) !{ISOP+O}
     &      + 0.276*R(193) !{TERP+OH}
     &      + 0.235*R(194) !{TERP+O3}  
C
C
      PA(nn) = prod_hcho_from_isop     
C
c
c... { HCHO Prod from VOC }
      nn =  nn + 1
      ptname(nn)  = 'HCHOp_Tot'
C
      PA(nn) = prod_hcho_from_isop
     &      + 1.000*R( 46) !{CXO2+NO}
     &      + 1.000*R( 48) !{CXO2+NO3}
     &      + 1.000*R( 49) !{CXO2+CXO2}
     &      + 2.000*R( 50) !{CXO2+CXO2}
     &      + 0.750*R( 54) !{RO2R+CXO2}
     &      + 0.750*R( 64) !{RO2N+CXO2}
     &      + 1.000*R( 74) !{CCO3+CXO2}
     &      + 1.000*R( 84) !{RCO3+CXO2}
     &      + 1.000*R( 95) !{BZCO+CXO2}
     &      + 1.000*R(104) !{MCO3+NO}
     &      + 1.000*R(106) !{MCO3+NO3}
     &      + 1.000*R(107) !{MCO3+CXO2}
     &      + 1.000*R(111) !{MCO3+CCO3}
     &      + 1.000*R(112) !{MCO3+RCO3}
     &      + 1.000*R(113) !{MCO3+BZCO}
     &      + 2.000*R(114) !{MCO3+MCO3}
     &      + 1.000*R(127) !{HCO3}
     &      + 1.000*R(136) !{ACET+OH}
     &      + 0.115*R(138) !{MEK+OH}
     &      + 1.000*R(140) !{MEOH+OH}
     &      + 0.350*R(141) !{COOH+OH}
     &      + 1.000*R(142) !{COOH}
     &      + 1.000*R(146) !{GLY}
     &      + 0.084*R(161) !{METH+OH}
     &      + 0.200*R(162) !{METH+O3}
     &      + 0.670*R(165) !{METH}
     &      + 0.300*R(166) !{MVK+OH}
     &      + 0.100*R(167) !{MVK+O3}
     &      + 0.213*R(174) !{PROD+OH}
     &      + 0.506*R(175) !{PROD}
     &      + 0.010*R(176) !{RNO3+OH}
     &      + 0.134*R(177) !{RNO3}
     &      + 1.610*R(185) !{ETHE+OH}
     &      + 1.000*R(186) !{ETHE+O3}
     &      + 0.191*R(188) !{ETHE+O}
     &      + 0.039*R(198) !{ALK2+OH}
     &      + 0.026*R(199) !{ALK3+OH}
     &      + 0.024*R(200) !{ALK4+OH}
     &      + 0.026*R(201) !{ALK5+OH}
     &      + 0.732*R(204) !{OLE1+OH}
     &      + 0.500*R(205) !{OLE1+O3}
     &      + 0.244*R(208) !{OLE2+OH}
     &      + 0.269*R(209) !{OLE2+O3}
     &      + 0.079*R(210) !{OLE2+NO3}
c
C
C
C
C            {     N O z     C h e m i s t r y      }
C
C
c... { HNO3 production from OH+NO2}
      nn =  nn + 1
      ptname(nn)  = 'HNO3_OHNO2'
C
      PA(nn) =  1.000*R( 25) !{OH+NO2}
C
C
C
C...  { HNO3 production from NO3 }
      nn = nn + 1
      ptname(nn)  = 'HNO3_NO3HC'
C
       PA(nn) = 0.200*R( 39) !{NO3+HO2}
     &        + 1.000*R(129) !{HCHO+NO3}
     &        + 1.000*R(132) !{CCHO+NO3}
     &        + 1.000*R(135) !{RCHO+NO3}
     &        + 1.000*R(148) !{GLY+NO3}
     &        + 1.000*R(151) !{MGLY+NO3}
     &        + 1.000*R(154) !{PHEN+NO3}
     &        + 1.000*R(156) !{CRES+NO3}
     &        + 1.000*R(157) !{NPHE+NO3}
     &        + 1.000*R(160) !{BALD+NO3}
     &        + 0.500*R(163) !{METH+NO3}
     &        + 0.150*R(172) !{ISPD+NO3}
C
C...  
C
      nn = nn + 1
      ptname(nn)  = 'HNO3_N2O5'
C
       PA(nn) =  2.000*R( 13) !{N2O5+H2O}
C
C

C...  net PAN prod
      nn = nn + 1
      ptname(nn)  = 'PANSprdNet'
C
      PA(nn) = cyc_PAN_p + cyc_PAN2_p + cyc_MPAN_p + cyc_PBZN_p
C
C
C
C...  net PAN loss
      nn = nn + 1
      ptname(nn)  = 'PANSlosNet'
C
      PA(nn) = cyc_PAN_l + cyc_PAN2_l + cyc_MPAN_l + cyc_PBZN_l
C
CC...  
C
      nn = nn + 1
      ptname(nn)  = 'NOzreact'
C
       PA(nn) =  1.000*R( 27) !{ 27} OH+HNO3=NO3+H2O
     &        +  1.000*R( 28)  !{ 28} HNO3=OH+NO2
     &        +  0.338*R(176)  !{176} RNO3+OH=0.338*NO2+0.310*RNO3+0.352*XN
     &        +  1.000*R(177)  !{177} RNO3=NO2C
c
c
c... { RNO3 Production }
      nn =  nn + 1
      ptname(nn)  = 'RNO3XNprod'
      PA(nn) = 
     &      + 1.000*R( 62) !{RO2N+NO}
     &      + 1.000*R(115) !{TBUO+NO2}
     &      + 0.572*R(172) !{ISPD+NO3}
     &      + 0.813*R(191) !{ISOP+NO3}
     &      + 0.276*R(195) !{TERP+NO3}
     &      + 0.511*R(206) !{OLE1+NO3}
     &      + 0.321*R(210) !{OLE2+NO3}
C
C
C...  XN prod
C      nn = nn + 1
C      ptname(nn)  = 'XNprod'
C
C       PA(nn) =     ....combine XN with RNO3 to save space
     &      + 2.000*R(120) !{120} BZNO+NO2
     &      + 0.500*R(163) !{{163} METH+NO3
     &      + 0.278*R(172) !{172} ISPD+NO3=0.278*XN
     &      + 1.000*R(187) !{187} ETHE+NO3
     &      + 0.250*R(195) !{195} TERP+NO3
     &      + 0.489*R(206) !{206} OLE1+NO3
     &      + 0.288*R(210) !{210} OLE2+NO3C
c     &      + 0.352*R(176) !{176} RNO3+OH  don't count this because we already counted RNO3
C
C...  
C
c      nn = nn + 1
c      ptname(nn)  = ''
C
c       PA(nn) =    
C
C
C   { convert to ppb/hr}
      Do n = 1,nn
        PA(n) = 1000.*PA(n)
      End do
C
C
      If( nn .GT. MXCPA ) Then
         write(iout,'(//,a)') 'ERROR in CPAMECH5:'
         write(iout,*) 'Number of outputs requested exceeds limit.'
         write(iout,*) 'Increase parameter MXCPA and try again.'
         call camxerr()
      End if

      NPA_INIT = nn  !set num of outputs on first dummy call from pasetup
      Return
      End
