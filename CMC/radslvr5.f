      subroutine radslvr5(ldark,H2O,M,O2,CH4,H2,cncrad,conc,r,crold,
     &                    dt)
c
c-----CAMx v4.02 030709
c
c     RADSLVR5 initializes radical concentrations
c
c     Copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Routines Called:
c        CPIVOT
c
c     Called by:
c        TRAP
c
      include "camx.prm"
      include "chmstry.com"
      include "filunit.com"
c
      logical ldark, lusecp, lsafe
      real conc(MXSPEC+1),cncrad(MXRADCL),r(MXRXN),rate(MXSPEC+1),
     &     loss(MXSPEC+1),gain(MXSPEC+1),jac(MXRADCL,MXRADCL),
     &     crold(MXRADCL)
      integer i, j, l, kount, nstrt, nend, nendcp, ierr
      real H2O, M, O2, CH4, H2
      real tol, weight, errbig, thresh, det, err, rlim, bnstrt, bnend
      real dt
      real gnNO3, lsNO3, lsN2O5, pdNO3, tpNO3, tpN2O5, fwd, bck, denom
      data tol/0.01/
c
c
c  new group of radicals
c
c   O1D
c
c
c  this radical is formed by photolysis only
c  solved by direct substitution
c
      if (ldark) then
        cncrad(kO1D) = 0.
      else
        Loss(kO1D  )= +( 1.000)*r( 19)+( 1.000)*r( 20)
        Gain(kO1D  )= +( 1.000)*r( 18)
        cncrad(kO1D) = gain(kO1D)/loss(kO1D)*cncrad(kO1D)
        cncrad(kO1D) = amax1(bdlrad, cncrad(kO1D))
      r( 19) = rk( 19)*cncrad(kO1D)*H2O
      r( 20) = rk( 20)*cncrad(kO1D)*M
      endif
c
c  new group of radicals
c
c     O
c
c
c  this radical is formed by photolysis only
c  solved by direct substitution
c
      if (ldark) then
        cncrad(kO) = 0.
      else
        Loss(kO    )= +( 1.000)*r(  2)+( 1.000)*r(  3)+( 1.000)*r(  4)
     &                +( 1.000)*r(  5)+( 1.000)*r(  6)+( 1.000)*r(164)
     &                +( 1.000)*r(168)+( 1.000)*r(188)+( 1.000)*r(192)
     &                +( 1.000)*r(196)+( 1.000)*r(207)+( 1.000)*r(211)
        Gain(kO    )= +( 1.000)*r(  1)+( 1.000)*r( 16)+( 1.000)*r( 17)
     &                +( 1.000)*r( 20)
        cncrad(kO) = gain(kO)/loss(kO)*cncrad(kO)
        cncrad(kO) = amax1(bdlrad, cncrad(kO))
      r(  2) = rk(  2)*cncrad(kO)*O2*M
      r(  3) = rk(  3)*cncrad(kO)*conc(kO3)
      r(  4) = rk(  4)*cncrad(kO)*conc(kNO)*M
      r(  5) = rk(  5)*cncrad(kO)*conc(kNO2)
      r(  6) = rk(  6)*cncrad(kO)*conc(kNO2)
      r(164) = rk(164)*conc(kMETH)*cncrad(kO)
      r(168) = rk(168)*conc(kMVK)*cncrad(kO)
      r(188) = rk(188)*conc(kETHE)*cncrad(kO)
      r(192) = rk(192)*conc(kISOP)*cncrad(kO)
      r(196) = rk(196)*conc(kTERP)*cncrad(kO)
      r(207) = rk(207)*conc(kOLE1)*cncrad(kO)
      r(211) = rk(211)*conc(kOLE2)*cncrad(kO)
      endif
cVAK
c  new group of radicals
c
c  N2O5   NO3
c
c  Explicitly solve N2O5 and NO3 using Hertel's solution
c  Note that NO3 self-reaction is hacked
cVAK
        gnNO3  =      +( 1.000)*r(  6)
     &                +( 1.000)*r(  8)
     &                +( 1.000)*r( 27)
     &                +( 0.390)*r( 34)
        lsN2O5 =      +( 1.000)*rk( 12)
     &                +( 1.000)*rk( 13)*H2O
        lsNO3  =      +( 1.000)*rk(  9)*conc(kNO)
     &                +( 1.000)*rk( 11)*conc(kNO2)
     &                +( 1.000)*rk( 14)*conc(kNO2)
     &                +( 1.000)*rk( 15)
     &                +( 1.000)*rk( 16)
     &                +( 1.000)*rk( 26)*cncrad(kOH)
     &                +( 1.000)*rk( 39)*cncrad(kHO2)
     &                +( 2.000)*rk( 40)*cncrad(kNO3)
     &                +( 1.000)*rk( 48)*cncrad(kCXO2)
     &                +( 1.000)*rk( 53)*cncrad(kRO2R)
     &                +( 1.000)*rk( 58)*cncrad(kR2O2)
     &                +( 1.000)*rk( 65)*cncrad(kRO2N)
     &                +( 1.000)*rk( 73)*cncrad(kCCO3)
     &                +( 1.000)*rk( 83)*cncrad(kRCO3)
     &                +( 1.000)*rk( 94)*cncrad(kBZCO)
     &                +( 1.000)*rk(106)*cncrad(kMCO3)
     &                +( 1.000)*rk(129)*conc(kHCHO)
     &                +( 1.000)*rk(132)*conc(kCCHO)
     &                +( 1.000)*rk(135)*conc(kRCHO)
     &                +( 1.000)*rk(148)*conc(kGLY)
     &                +( 1.000)*rk(151)*conc(kMGLY)
     &                +( 1.000)*rk(154)*conc(kPHEN)
     &                +( 1.000)*rk(156)*conc(kCRES)
     &                +( 1.000)*rk(157)*conc(kNPHE)
     &                +( 1.000)*rk(160)*conc(kBALD)
     &                +( 1.000)*rk(163)*conc(kMETH)
     &                +( 1.000)*rk(172)*conc(kISPD)
     &                +( 1.000)*rk(187)*conc(kETHE)
     &                +( 1.000)*rk(191)*conc(kISOP)
     &                +( 1.000)*rk(195)*conc(kTERP)
     &                +( 1.000)*rk(206)*conc(kOLE1)
     &               +( 1.000)*rk(210)*conc(kOLE2)
      lsNO3 = 1.0 + dt*lsNO3
      lsN2O5 = 1.0 + dt*lsN2O5
      pdNO3 = crold(kNO3) + dt*gnNO3
      fwd = dt*rk( 12)
      bck = dt*rk( 11)*conc(kNO2)
      tpNO3 = lsN2O5*pdNO3 + fwd*crold(kN2O5)
      tpN2O5 = lsNO3*crold(kN2O5) + bck*pdNO3
      denom = (lsN2O5*lsNO3) - (fwd*bck)
      cncrad(kNO3) = tpNO3/denom
      cncrad(kN2O5) = tpN2O5/denom
      r(  9) = rk(  9)*conc(kNO)*cncrad(kNO3)
      r( 11) = rk( 11)*conc(kNO2)*cncrad(kNO3)
      r( 12) = rk( 12)*cncrad(kN2O5)
      r( 13) = rk( 13)*cncrad(kN2O5)*H2O
      r( 14) = rk( 14)*conc(kNO2)*cncrad(kNO3)
      r( 15) = rk( 15)*cncrad(kNO3)
      r( 16) = rk( 16)*cncrad(kNO3)
      r( 26) = rk( 26)*cncrad(kOH)*cncrad(kNO3)
      r( 39) = rk( 39)*cncrad(kNO3)*cncrad(kHO2)
      r( 40) = rk( 40)*cncrad(kNO3)*cncrad(kNO3)
      r( 48) = rk( 48)*cncrad(kCXO2)*cncrad(kNO3)
      r( 53) = rk( 53)*cncrad(kRO2R)*cncrad(kNO3)
      r( 58) = rk( 58)*cncrad(kR2O2)*cncrad(kNO3)
      r( 65) = rk( 65)*cncrad(kRO2N)*cncrad(kNO3)
      r( 73) = rk( 73)*cncrad(kCCO3)*cncrad(kNO3)
      r( 83) = rk( 83)*cncrad(kRCO3)*cncrad(kNO3)
      r( 94) = rk( 94)*cncrad(kBZCO)*cncrad(kNO3)
      r(106) = rk(106)*cncrad(kMCO3)*cncrad(kNO3)
      r(129) = rk(129)*conc(kHCHO)*cncrad(kNO3)
      r(132) = rk(132)*conc(kCCHO)*cncrad(kNO3)
      r(135) = rk(135)*conc(kRCHO)*cncrad(kNO3)
      r(148) = rk(148)*conc(kGLY)*cncrad(kNO3)
      r(151) = rk(151)*conc(kMGLY)*cncrad(kNO3)
      r(154) = rk(154)*conc(kPHEN)*cncrad(kNO3)
      r(156) = rk(156)*conc(kCRES)*cncrad(kNO3)
      r(157) = rk(157)*conc(kNPHE)*cncrad(kNO3)
      r(160) = rk(160)*conc(kBALD)*cncrad(kNO3)
      r(163) = rk(163)*conc(kMETH)*cncrad(kNO3)
      r(172) = rk(172)*conc(kISPD)*cncrad(kNO3)
      r(187) = rk(187)*conc(kETHE)*cncrad(kNO3)
      r(191) = rk(191)*conc(kISOP)*cncrad(kNO3)
      r(195) = rk(195)*conc(kTERP)*cncrad(kNO3)
      r(206) = rk(206)*conc(kOLE1)*cncrad(kNO3)
      r(210) = rk(210)*conc(kOLE2)*cncrad(kNO3)
c
c  new group of radicals
c
c    OH   HO2  RO2R  R2O2  RO2N  CCO3  RCO3  MCO3  BZCO  CXO2
c  HCO3  TBUO   BZO  BZNO
c
      nstrt = kOH
      nend  = kBZNO
      weight = 1.0
      errbig = 1.0
      thresh = 1.0e-15
      lsafe = .false.
      kount = 0
 14   kount = kount + 1
      if (kount.gt.500) then
        write(iout,'(//,A,//)') 'ERROR in RADSLVR5:'
        write(iout,*) 'No Convergence in RADSLVR5, errbig = ', errbig
        write(iout,*) 'statement number = ', 14
        write(iout,*) 'igrd,i, j, k = ', igrdchm,ichm,jchm,kchm
        write(iout,*) 'LDARK is set ', ldark
        do l=1,ngas
          write(iout,'(i3,2x,a7,1pe10.3)') l,spname(l),conc(l)
        enddo
        write(iout,*) 'The radicals are: '
        do l= 1 , nrad
          write(iout,'(i3,2x,a7,1pe10.3)')
     &          l,nmrad(l),cncrad(l)
        enddo
        write(iout,*) 'Currently solving ', nstrt, ' to ',nendcp
        do l=nstrt,nendcp
          write(iout,'(i3,2x,a7,1p2e10.3)')
     &          l,nmrad(l),cncrad(l),abs(rate(l))/cncrad(l)
        enddo
        call camxerr()
      endif
      do i=nstrt,nend
        do j=nstrt,nend
          jac(i,j) = 0.
        enddo
      enddo
c
c  to reduce the matrix size, HCO3 is solved first and substituted
c
        Loss(kHCO3 )= +( 1.000)*r(127)+( 1.000)*r(128)
        Gain(kHCO3 )= +( 1.000)*r(126)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kHCO3) = gain(kHCO3)/loss(kHCO3)*cncrad(kHCO3)
      cncrad(kHCO3) = amax1(bdlrad, cncrad(kHCO3))
      r(127) = rk(127)*cncrad(kHCO3)
      r(128) = rk(128)*cncrad(kHCO3)*conc(kNO)
c
c  to reduce the matrix size, TBUO is solved first and substituted
c
        Loss(kTBUO )= +( 1.000)*r(115)+( 1.000)*r(116)
        Gain(kTBUO )= +( 0.236)*r(199)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kTBUO) = gain(kTBUO)/loss(kTBUO)*cncrad(kTBUO)
      cncrad(kTBUO) = amax1(bdlrad, cncrad(kTBUO))
      r(115) = rk(115)*cncrad(kTBUO)*conc(kNO2)
      r(116) = rk(116)*cncrad(kTBUO)
c
c  to reduce the matrix size, BZO is solved first and substituted
c
        Loss(kBZO  )= +( 1.000)*r(117)+( 1.000)*r(118)+( 1.000)*r(119)
        Gain(kBZO  )= +( 1.000)*r( 92)+( 1.000)*r( 94)+( 1.000)*r( 99)
     &                +( 1.000)*r(100)+( 2.000)*r(101)+( 1.000)*r(113)
     &                +( 0.240)*r(153)+( 1.000)*r(154)+( 0.240)*r(155)
     &                +( 1.000)*r(156)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kBZO) = gain(kBZO)/loss(kBZO)*cncrad(kBZO)
      cncrad(kBZO) = amax1(bdlrad, cncrad(kBZO))
      r(117) = rk(117)*cncrad(kBZO)*conc(kNO2)
      r(118) = rk(118)*cncrad(kBZO)*cncrad(kHO2)
      r(119) = rk(119)*cncrad(kBZO)
c
c  to reduce the matrix size, BZNO is solved first and substituted
c
        Loss(kBZNO )= +( 1.000)*r(120)+( 1.000)*r(121)+( 1.000)*r(122)
        Gain(kBZNO )= +( 1.000)*r(157)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kBZNO) = gain(kBZNO)/loss(kBZNO)*cncrad(kBZNO)
      cncrad(kBZNO) = amax1(bdlrad, cncrad(kBZNO))
      r(120) = rk(120)*cncrad(kBZNO)*conc(kNO2)
      r(121) = rk(121)*cncrad(kBZNO)*cncrad(kHO2)
      r(122) = rk(122)*cncrad(kBZNO)
        Loss(kOH   )= +( 1.000)*r( 21)+( 1.000)*r( 24)+( 1.000)*r( 25)
     &                +( 1.000)*r( 26)+( 1.000)*r( 27)+( 1.000)*r( 29)
     &                +( 1.000)*r( 30)+( 1.000)*r( 35)+( 1.000)*r( 42)
     &                +( 1.000)*r( 43)+( 1.000)*r( 44)+( 1.000)*r( 45)
     &                +( 1.000)*r(125)+( 1.000)*r(130)+( 1.000)*r(133)
     &                +( 1.000)*r(136)+( 1.000)*r(138)+( 1.000)*r(140)
     &                +( 1.000)*r(141)+( 1.000)*r(143)+( 1.000)*r(147)
     &                +( 1.000)*r(150)+( 1.000)*r(153)+( 1.000)*r(155)
     &                +( 1.000)*r(158)+( 1.000)*r(161)+( 1.000)*r(166)
        Loss(kOH   ) = Loss(kOH   )
     &                +( 1.000)*r(170)+( 1.000)*r(174)+( 1.000)*r(176)
     &                +( 1.000)*r(178)+( 1.000)*r(180)+( 1.000)*r(182)
     &                +( 1.000)*r(184)+( 1.000)*r(185)+( 1.000)*r(189)
     &                +( 1.000)*r(193)+( 1.000)*r(197)+( 1.000)*r(198)
     &                +( 1.000)*r(199)+( 1.000)*r(200)+( 1.000)*r(201)
     &                +( 1.000)*r(202)+( 1.000)*r(203)+( 1.000)*r(204)
     &                +( 1.000)*r(208)+( 1.000)*r(212)+( 1.000)*r(213)
     &                +( 1.000)*r(214)+( 1.000)*r(215)+( 1.000)*r(216)
     &                +( 1.000)*r(217)+( 1.000)*r(218)+( 1.000)*r(219)
     &                +( 1.000)*r(220)+( 1.000)*r(221)+( 1.000)*r(222)
     &                +( 1.000)*r(223)+( 1.000)*r(224)+( 1.000)*r(225)
     &                +( 1.000)*r(226)+( 1.000)*r(227)+( 1.000)*r(228)
     &                +( 1.000)*r(229)+( 1.000)*r(230)+( 1.000)*r(231)
     &                +( 1.000)*r(232)+( 1.000)*r(233)+( 1.000)*r(234)
     &                +( 1.000)*r(235)+( 1.000)*r(236)+( 1.000)*r(237)
        Gain(kOH   )= +( 2.000)*r( 19)+( 1.000)*r( 22)+( 1.000)*r( 28)
     &                +( 1.000)*r( 31)+( 0.390)*r( 34)+( 1.000)*r( 36)
     &                +( 0.800)*r( 39)+( 2.000)*r( 41)+( 0.350)*r(141)
     &                +( 1.000)*r(142)+( 0.660)*r(143)+( 1.000)*r(144)
     &                +( 0.208)*r(162)+( 0.330)*r(165)+( 0.164)*r(167)
     &                +( 0.285)*r(171)+( 0.500)*r(179)+( 0.120)*r(186)
     &                +( 0.266)*r(190)+( 0.567)*r(194)+( 0.246)*r(198)
     &                +( 0.155)*r(205)+( 0.378)*r(209)
        Loss(kHO2  )= +( 1.000)*r( 31)+( 1.000)*r( 32)+( 1.000)*r( 36)
     &                +( 2.000)*r( 37)+( 2.000)*r( 38)+( 1.000)*r( 39)
     &                +( 1.000)*r( 43)+( 1.000)*r( 47)+( 1.000)*r( 52)
     &                +( 1.000)*r( 57)+( 1.000)*r( 63)+( 1.000)*r( 72)
     &                +( 1.000)*r( 82)+( 1.000)*r( 93)+( 1.000)*r(105)
     &                +( 1.000)*r(118)+( 1.000)*r(121)+( 1.000)*r(126)
        Gain(kHO2  )= +( 1.000)*r( 23)+( 1.000)*r( 26)+( 1.000)*r( 29)
     &                +( 1.000)*r( 30)+( 1.000)*r( 33)+( 0.610)*r( 34)
     &                +( 1.000)*r( 42)+( 1.000)*r( 44)+( 1.000)*r( 45)
     &                +( 1.000)*r( 46)+( 1.000)*r( 48)+( 2.000)*r( 50)
     &                +( 1.000)*r( 51)+( 1.000)*r( 53)+( 1.000)*r( 54)
     &                +( 1.000)*r( 55)+( 1.000)*r( 57)+( 1.000)*r( 64)
     &                +( 1.000)*r( 65)+( 1.000)*r( 66)+( 1.000)*r( 68)
     &                +( 2.000)*r(123)+( 1.000)*r(125)+( 1.000)*r(127)
     &                +( 1.000)*r(128)+( 1.000)*r(129)+( 1.000)*r(131)
        Gain(kHO2  ) = Gain(kHO2  )
     &                +( 1.000)*r(134)+( 1.000)*r(140)+( 1.000)*r(142)
     &                +( 1.000)*r(144)+( 2.000)*r(145)+( 0.630)*r(147)
     &                +( 0.630)*r(148)+( 1.000)*r(149)+( 0.008)*r(162)
     &                +( 0.340)*r(165)+( 0.064)*r(167)+( 0.400)*r(171)
     &                +( 1.233)*r(173)+( 0.379)*r(174)+( 0.113)*r(176)
     &                +( 0.341)*r(177)+( 1.500)*r(179)+( 0.500)*r(181)
     &                +( 0.500)*r(183)+( 0.120)*r(186)+( 0.500)*r(188)
     &                +( 0.033)*r(194)+( 0.121)*r(198)+( 0.224)*r(202)
        Gain(kHO2  ) = Gain(kHO2  )
     &                +( 0.187)*r(203)+( 0.056)*r(205)+( 0.003)*r(209)
     &                +( 0.013)*r(211)
        Loss(kRO2R )= +( 1.000)*r( 51)+( 1.000)*r( 52)+( 1.000)*r( 53)
     &                +( 1.000)*r( 54)+( 2.000)*r( 55)+( 1.000)*r( 60)
     &                +( 1.000)*r( 66)+( 1.000)*r( 75)+( 1.000)*r( 85)
     &                +( 1.000)*r( 96)+( 1.000)*r(108)
        Gain(kRO2R )= +( 1.000)*r( 60)+( 1.000)*r( 81)+( 1.000)*r( 83)
     &                +( 1.000)*r( 88)+( 2.000)*r( 89)+( 1.000)*r(100)
     &                +( 1.000)*r(112)+( 0.034)*r(133)+( 1.000)*r(134)
     &                +( 0.370)*r(138)+( 1.000)*r(139)+( 0.340)*r(143)
     &                +( 0.760)*r(153)+( 0.760)*r(155)+( 0.500)*r(161)
     &                +( 0.100)*r(162)+( 0.500)*r(163)+( 0.330)*r(165)
     &                +( 0.300)*r(166)+( 0.050)*r(167)+( 0.670)*r(170)
     &                +( 0.048)*r(171)+( 0.799)*r(172)+( 0.473)*r(174)
     &                +( 0.960)*r(175)+( 0.376)*r(176)+( 0.564)*r(177)
        Gain(kRO2R ) = Gain(kRO2R )
     &                +( 1.000)*r(178)+( 1.000)*r(181)+( 1.000)*r(183)
     &                +( 1.000)*r(185)+( 1.000)*r(187)+( 0.200)*r(188)
     &                +( 0.907)*r(189)+( 0.066)*r(190)+( 0.749)*r(191)
     &                +( 0.750)*r(193)+( 0.031)*r(194)+( 0.276)*r(195)
     &                +( 1.000)*r(197)+( 0.612)*r(198)+( 0.695)*r(199)
     &                +( 0.835)*r(200)+( 0.653)*r(201)+( 0.765)*r(202)
     &                +( 0.804)*r(203)+( 0.910)*r(204)+( 0.022)*r(205)
     &                +( 0.824)*r(206)+( 0.918)*r(208)+( 0.033)*r(209)
        Gain(kRO2R ) = Gain(kRO2R )
     &                +( 0.442)*r(210)+( 0.012)*r(211)
        Loss(kR2O2 )= +( 1.000)*r( 56)+( 1.000)*r( 57)+( 1.000)*r( 58)
     &                +( 1.000)*r( 59)+( 1.000)*r( 60)+( 2.000)*r( 61)
     &                +( 1.000)*r( 67)+( 1.000)*r( 76)+( 1.000)*r( 86)
     &                +( 1.000)*r( 97)+( 1.000)*r(109)
        Gain(kR2O2 )= +( 1.000)*r( 92)+( 1.000)*r( 94)+( 1.000)*r( 99)
     &                +( 1.000)*r(100)+( 2.000)*r(101)+( 1.000)*r(113)
     &                +( 1.000)*r(136)+( 0.616)*r(138)+( 0.675)*r(166)
     &                +( 0.515)*r(175)+( 0.596)*r(176)+( 0.152)*r(177)
     &                +( 1.000)*r(180)+( 1.000)*r(181)+( 1.000)*r(182)
     &                +( 1.000)*r(183)+( 0.079)*r(189)+( 0.126)*r(190)
     &                +( 0.187)*r(191)+( 0.240)*r(192)+( 0.500)*r(193)
     &                +( 0.729)*r(194)+( 0.750)*r(195)+( 0.559)*r(199)
     &                +( 0.936)*r(200)+( 0.948)*r(201)+( 0.205)*r(204)
        Gain(kR2O2 ) = Gain(kR2O2 )
     &                +( 0.488)*r(206)+( 0.001)*r(208)+( 0.137)*r(209)
     &                +( 0.711)*r(210)
        Loss(kRO2N )= +( 1.000)*r( 62)+( 1.000)*r( 63)+( 1.000)*r( 64)
     &                +( 1.000)*r( 65)+( 1.000)*r( 66)+( 1.000)*r( 67)
     &                +( 2.000)*r( 68)+( 1.000)*r( 77)+( 1.000)*r( 87)
     &                +( 1.000)*r( 98)+( 1.000)*r(110)
        Gain(kRO2N )= +( 1.000)*r( 67)+( 0.001)*r(133)+( 0.042)*r(138)
     &                +( 0.025)*r(166)+( 0.041)*r(170)+( 0.051)*r(172)
     &                +( 0.070)*r(174)+( 0.040)*r(175)+( 0.173)*r(176)
     &                +( 0.095)*r(177)+( 0.093)*r(189)+( 0.008)*r(190)
     &                +( 0.064)*r(191)+( 0.010)*r(192)+( 0.250)*r(193)
     &                +( 0.180)*r(194)+( 0.250)*r(195)+( 0.021)*r(198)
     &                +( 0.070)*r(199)+( 0.143)*r(200)+( 0.347)*r(201)
     &                +( 0.011)*r(202)+( 0.009)*r(203)+( 0.090)*r(204)
     &                +( 0.001)*r(205)+( 0.176)*r(206)+( 0.082)*r(208)
        Gain(kRO2N ) = Gain(kRO2N )
     &                +( 0.002)*r(209)+( 0.136)*r(210)+( 0.001)*r(211)
        Loss(kCCO3 )= +( 1.000)*r( 69)+( 1.000)*r( 71)+( 1.000)*r( 72)
     &                +( 1.000)*r( 73)+( 1.000)*r( 74)+( 1.000)*r( 75)
     &                +( 1.000)*r( 76)+( 1.000)*r( 77)+( 2.000)*r( 78)
     &                +( 1.000)*r( 88)+( 1.000)*r( 99)+( 1.000)*r(111)
        Gain(kCCO3 )= +( 1.000)*r( 70)+( 1.000)*r( 76)+( 1.000)*r(104)
     &                +( 1.000)*r(106)+( 1.000)*r(111)+( 1.000)*r(112)
     &                +( 1.000)*r(113)+( 2.000)*r(114)+( 1.000)*r(130)
     &                +( 1.000)*r(132)+( 1.000)*r(136)+( 1.000)*r(137)
     &                +( 0.492)*r(138)+( 1.000)*r(139)+( 1.000)*r(149)
     &                +( 1.000)*r(150)+( 1.000)*r(151)+( 2.000)*r(152)
     &                +( 0.670)*r(165)+( 0.675)*r(166)+( 0.467)*r(173)
     &                +( 0.029)*r(174)+( 0.667)*r(175)+( 1.000)*r(180)
     &                +( 0.500)*r(181)+( 1.000)*r(182)+( 0.500)*r(183)
        Gain(kCCO3 ) = Gain(kCCO3 )
     &                +( 0.123)*r(194)+( 0.011)*r(200)+( 0.137)*r(209)
        Loss(kRCO3 )= +( 1.000)*r( 79)+( 1.000)*r( 81)+( 1.000)*r( 82)
     &                +( 1.000)*r( 83)+( 1.000)*r( 84)+( 1.000)*r( 85)
     &                +( 1.000)*r( 86)+( 1.000)*r( 87)+( 1.000)*r( 88)
     &                +( 2.000)*r( 89)+( 1.000)*r(100)+( 1.000)*r(112)
        Gain(kRCO3 )= +( 1.000)*r( 80)+( 1.000)*r( 86)+( 0.965)*r(133)
     &                +( 1.000)*r(135)+( 0.096)*r(138)+( 0.370)*r(147)
     &                +( 0.370)*r(148)+( 0.100)*r(162)+( 0.050)*r(167)
     &                +( 0.048)*r(171)+( 0.300)*r(173)+( 0.049)*r(174)
     &                +( 0.333)*r(175)+( 0.201)*r(194)+( 0.006)*r(209)
        Loss(kMCO3 )= +( 1.000)*r(102)+( 1.000)*r(104)+( 1.000)*r(105)
     &                +( 1.000)*r(106)+( 1.000)*r(107)+( 1.000)*r(108)
     &                +( 1.000)*r(109)+( 1.000)*r(110)+( 1.000)*r(111)
     &                +( 1.000)*r(112)+( 1.000)*r(113)+( 2.000)*r(114)
        Gain(kMCO3 )= +( 1.000)*r(103)+( 1.000)*r(109)+( 0.500)*r(161)
     &                +( 0.500)*r(163)+( 0.330)*r(165)+( 0.300)*r(169)
     &                +( 0.289)*r(170)+( 0.150)*r(172)+( 0.192)*r(190)
     &                +( 0.240)*r(192)
        Loss(kBZCO )= +( 1.000)*r( 90)+( 1.000)*r( 92)+( 1.000)*r( 93)
     &                +( 1.000)*r( 94)+( 1.000)*r( 95)+( 1.000)*r( 96)
     &                +( 1.000)*r( 97)+( 1.000)*r( 98)+( 1.000)*r( 99)
     &                +( 1.000)*r(100)+( 2.000)*r(101)+( 1.000)*r(113)
        Gain(kBZCO )= +( 1.000)*r( 91)+( 1.000)*r( 97)+( 1.000)*r(158)
     &                +( 1.000)*r(160)
        Loss(kCXO2 )= +( 1.000)*r( 46)+( 1.000)*r( 47)+( 1.000)*r( 48)
     &                +( 2.000)*r( 49)+( 2.000)*r( 50)+( 1.000)*r( 54)
     &                +( 1.000)*r( 59)+( 1.000)*r( 64)+( 1.000)*r( 74)
     &                +( 1.000)*r( 84)+( 1.000)*r( 95)+( 1.000)*r(107)
        Gain(kCXO2 )= +( 1.000)*r( 59)+( 1.000)*r( 71)+( 1.000)*r( 73)
     &                +( 2.000)*r( 78)+( 1.000)*r( 88)+( 1.000)*r( 99)
     &                +( 1.000)*r(111)+( 1.000)*r(116)+( 1.000)*r(131)
     &                +( 1.000)*r(137)+( 0.650)*r(141)+( 0.300)*r(169)
     &                +( 1.000)*r(184)+( 0.300)*r(188)+( 0.250)*r(192)
     &                +( 0.011)*r(200)+( 0.076)*r(205)+( 0.197)*r(209)
     &                +( 0.030)*r(210)
        Loss(kHCO3 )= +( 1.000)*r(127)+( 1.000)*r(128)
        Gain(kHCO3 )= +( 1.000)*r(126)
        Loss(kTBUO )= +( 1.000)*r(115)+( 1.000)*r(116)
        Gain(kTBUO )= +( 0.236)*r(199)
        Loss(kBZO  )= +( 1.000)*r(117)+( 1.000)*r(118)+( 1.000)*r(119)
        Gain(kBZO  )= +( 1.000)*r( 92)+( 1.000)*r( 94)+( 1.000)*r( 99)
     &                +( 1.000)*r(100)+( 2.000)*r(101)+( 1.000)*r(113)
     &                +( 0.240)*r(153)+( 1.000)*r(154)+( 0.240)*r(155)
     &                +( 1.000)*r(156)
        Loss(kBZNO )= +( 1.000)*r(120)+( 1.000)*r(121)+( 1.000)*r(122)
        Gain(kBZNO )= +( 1.000)*r(157)
 
          JAC(kOH  ,kOH  )= +( 1.000)*r( 21)+( 1.000)*r( 24)
     &                      +( 1.000)*r( 25)+( 1.000)*r( 26)
     &                      +( 1.000)*r( 27)+( 1.000)*r( 29)
     &                      +( 1.000)*r( 30)+( 1.000)*r( 35)
     &                      +( 1.000)*r( 42)+( 1.000)*r( 43)
     &                      +( 1.000)*r( 44)+( 1.000)*r( 45)
     &                      +( 1.000)*r(125)+( 1.000)*r(130)
     &                      +( 1.000)*r(133)+( 1.000)*r(136)
     &                      +( 1.000)*r(138)+( 1.000)*r(140)
          JAC(kOH  ,kOH  )=JAC(kOH  ,kOH  )
     &                      +( 1.000)*r(141)+(-0.350)*r(141)
     &                      +( 1.000)*r(143)+(-0.660)*r(143)
     &                      +( 1.000)*r(147)+( 1.000)*r(150)
     &                      +( 1.000)*r(153)+( 1.000)*r(155)
     &                      +( 1.000)*r(158)+( 1.000)*r(161)
     &                      +( 1.000)*r(166)+( 1.000)*r(170)
     &                      +( 1.000)*r(174)+( 1.000)*r(176)
     &                      +( 1.000)*r(178)+( 1.000)*r(180)
          JAC(kOH  ,kOH  )=JAC(kOH  ,kOH  )
     &                      +( 1.000)*r(182)+( 1.000)*r(184)
     &                      +( 1.000)*r(185)+( 1.000)*r(189)
     &                      +( 1.000)*r(193)+( 1.000)*r(197)
     &                      +( 1.000)*r(198)+(-0.246)*r(198)
     &                      +( 1.000)*r(199)+( 1.000)*r(200)
     &                      +( 1.000)*r(201)+( 1.000)*r(202)
     &                      +( 1.000)*r(203)+( 1.000)*r(204)
     &                      +( 1.000)*r(208)+( 1.000)*r(212)
     &                      +( 1.000)*r(213)+( 1.000)*r(214)
     &                      +( 1.000)*r(215)+( 1.000)*r(216)
     &                      +( 1.000)*r(217)+( 1.000)*r(218)
     &                      +( 1.000)*r(219)+( 1.000)*r(220)
     &                      +( 1.000)*r(221)+( 1.000)*r(222)
     &                      +( 1.000)*r(223)+( 1.000)*r(224)
     &                      +( 1.000)*r(225)+( 1.000)*r(226)
     &                      +( 1.000)*r(227)+( 1.000)*r(228)
     &                      +( 1.000)*r(229)+( 1.000)*r(230)
     &                      +( 1.000)*r(231)+( 1.000)*r(232)
     &                      +( 1.000)*r(233)+( 1.000)*r(234)
     &                      +( 1.000)*r(235)+( 1.000)*r(236)
     &                      +( 1.000)*r(237)
          JAC(kOH  ,kOH  )=JAC(kOH  ,kOH  )
          JAC(kOH  ,kHO2 )= +(-1.000)*r( 31)+(-1.000)*r( 36)
     &                      +(-0.800)*r( 39)+( 1.000)*r( 43)
          JAC(kHO2 ,kOH  )= +(-1.000)*r( 26)+(-1.000)*r( 29)
     &                      +(-1.000)*r( 30)+(-1.000)*r( 42)
     &                      +( 1.000)*r( 43)+(-1.000)*r( 44)
     &                      +(-1.000)*r( 45)+(-1.000)*r(125)
     &                      +(-1.000)*r(140)+(-0.630)*r(147)
     &                      +(-0.379)*r(174)+(-0.113)*r(176)
     &                      +(-0.121)*r(198)+(-0.224)*r(202)
     &                      +(-0.187)*r(203)
          JAC(kHO2 ,kHO2 )= +( 1.000)*r( 31)+( 1.000)*r( 32)
     &                      +( 1.000)*r( 36)+( 4.000)*r( 37)
     &                      +( 4.000)*r( 38)+( 1.000)*r( 39)
     &                      +( 1.000)*r( 43)+( 1.000)*r( 47)
     &                      +( 1.000)*r( 52)+( 1.000)*r( 57)
     &                      +(-1.000)*r( 57)+( 1.000)*r( 63)
     &                      +( 1.000)*r( 72)+( 1.000)*r( 82)
     &                      +( 1.000)*r( 93)+( 1.000)*r(105)
     &                      +( 1.000)*r(118)+( 1.000)*r(121)
          JAC(kHO2 ,kHO2 )=JAC(kHO2 ,kHO2 )
     &                      +( 1.000)*r(126)
          JAC(kHO2 ,kRO2R)= +(-1.000)*r( 51)+( 1.000)*r( 52)
     &                      +(-1.000)*r( 53)+(-1.000)*r( 54)
     &                      +(-2.000)*r( 55)+(-1.000)*r( 66)
          JAC(kHO2 ,kR2O2)= +( 1.000)*r( 57)+(-1.000)*r( 57)
          JAC(kHO2 ,kRO2N)= +( 1.000)*r( 63)+(-1.000)*r( 64)
     &                      +(-1.000)*r( 65)+(-1.000)*r( 66)
     &                      +(-2.000)*r( 68)
          JAC(kHO2 ,kCCO3)= +( 1.000)*r( 72)
          JAC(kHO2 ,kRCO3)= +( 1.000)*r( 82)
          JAC(kHO2 ,kMCO3)= +( 1.000)*r(105)
          JAC(kHO2 ,kBZCO)= +( 1.000)*r( 93)
          JAC(kHO2 ,kCXO2)= +(-1.000)*r( 46)+( 1.000)*r( 47)
     &                      +(-1.000)*r( 48)+(-4.000)*r( 50)
     &                      +(-1.000)*r( 54)+(-1.000)*r( 64)
          JAC(kHO2 ,kHCO3)= +(-1.000)*r(127)+(-1.000)*r(128)
          JAC(kHO2 ,kBZO )= +( 1.000)*r(118)
          JAC(kHO2 ,kBZNO)= +( 1.000)*r(121)
          JAC(kRO2R,kOH  )= +(-0.034)*r(133)+(-0.370)*r(138)
     &                      +(-0.340)*r(143)+(-0.760)*r(153)
     &                      +(-0.760)*r(155)+(-0.500)*r(161)
     &                      +(-0.300)*r(166)+(-0.670)*r(170)
     &                      +(-0.473)*r(174)+(-0.376)*r(176)
     &                      +(-1.000)*r(178)+(-1.000)*r(185)
     &                      +(-0.907)*r(189)+(-0.750)*r(193)
     &                      +(-1.000)*r(197)+(-0.612)*r(198)
     &                      +(-0.695)*r(199)+(-0.835)*r(200)
          JAC(kRO2R,kOH  )=JAC(kRO2R,kOH  )
     &                      +(-0.653)*r(201)+(-0.765)*r(202)
     &                      +(-0.804)*r(203)+(-0.910)*r(204)
     &                      +(-0.918)*r(208)
          JAC(kRO2R,kHO2 )= +( 1.000)*r( 52)
          JAC(kRO2R,kRO2R)= +( 1.000)*r( 51)+( 1.000)*r( 52)
     &                      +( 1.000)*r( 53)+( 1.000)*r( 54)
     &                      +( 4.000)*r( 55)+( 1.000)*r( 60)
     &                      +(-1.000)*r( 60)+( 1.000)*r( 66)
     &                      +( 1.000)*r( 75)+( 1.000)*r( 85)
     &                      +( 1.000)*r( 96)+( 1.000)*r(108)
          JAC(kRO2R,kR2O2)= +( 1.000)*r( 60)+(-1.000)*r( 60)
          JAC(kRO2R,kRO2N)= +( 1.000)*r( 66)
          JAC(kRO2R,kCCO3)= +( 1.000)*r( 75)+(-1.000)*r( 88)
          JAC(kRO2R,kRCO3)= +(-1.000)*r( 81)+(-1.000)*r( 83)
     &                      +( 1.000)*r( 85)+(-1.000)*r( 88)
     &                      +(-4.000)*r( 89)+(-1.000)*r(100)
     &                      +(-1.000)*r(112)
          JAC(kRO2R,kMCO3)= +( 1.000)*r(108)+(-1.000)*r(112)
          JAC(kRO2R,kBZCO)= +( 1.000)*r( 96)+(-1.000)*r(100)
          JAC(kRO2R,kCXO2)= +( 1.000)*r( 54)
          JAC(kR2O2,kOH  )= +(-1.000)*r(136)+(-0.616)*r(138)
     &                      +(-0.675)*r(166)+(-0.596)*r(176)
     &                      +(-1.000)*r(180)+(-1.000)*r(182)
     &                      +(-0.079)*r(189)+(-0.500)*r(193)
     &                      +(-0.559)*r(199)+(-0.936)*r(200)
     &                      +(-0.948)*r(201)+(-0.205)*r(204)
     &                      +(-0.001)*r(208)
          JAC(kR2O2,kHO2 )= +( 1.000)*r( 57)
          JAC(kR2O2,kRO2R)= +( 1.000)*r( 60)
          JAC(kR2O2,kR2O2)= +( 1.000)*r( 56)+( 1.000)*r( 57)
     &                      +( 1.000)*r( 58)+( 1.000)*r( 59)
     &                      +( 1.000)*r( 60)+( 4.000)*r( 61)
     &                      +( 1.000)*r( 67)+( 1.000)*r( 76)
     &                      +( 1.000)*r( 86)+( 1.000)*r( 97)
     &                      +( 1.000)*r(109)
          JAC(kR2O2,kRO2N)= +( 1.000)*r( 67)
          JAC(kR2O2,kCCO3)= +( 1.000)*r( 76)+(-1.000)*r( 99)
          JAC(kR2O2,kRCO3)= +( 1.000)*r( 86)+(-1.000)*r(100)
          JAC(kR2O2,kMCO3)= +( 1.000)*r(109)+(-1.000)*r(113)
          JAC(kR2O2,kBZCO)= +(-1.000)*r( 92)+(-1.000)*r( 94)
     &                      +( 1.000)*r( 97)+(-1.000)*r( 99)
     &                      +(-1.000)*r(100)+(-4.000)*r(101)
     &                      +(-1.000)*r(113)
          JAC(kR2O2,kCXO2)= +( 1.000)*r( 59)
          JAC(kRO2N,kOH  )= +(-0.001)*r(133)+(-0.042)*r(138)
     &                      +(-0.025)*r(166)+(-0.041)*r(170)
     &                      +(-0.070)*r(174)+(-0.173)*r(176)
     &                      +(-0.093)*r(189)+(-0.250)*r(193)
     &                      +(-0.021)*r(198)+(-0.070)*r(199)
     &                      +(-0.143)*r(200)+(-0.347)*r(201)
     &                      +(-0.011)*r(202)+(-0.009)*r(203)
     &                      +(-0.090)*r(204)+(-0.082)*r(208)
          JAC(kRO2N,kHO2 )= +( 1.000)*r( 63)
          JAC(kRO2N,kRO2R)= +( 1.000)*r( 66)
          JAC(kRO2N,kR2O2)= +( 1.000)*r( 67)+(-1.000)*r( 67)
          JAC(kRO2N,kRO2N)= +( 1.000)*r( 62)+( 1.000)*r( 63)
     &                      +( 1.000)*r( 64)+( 1.000)*r( 65)
     &                      +( 1.000)*r( 66)+( 1.000)*r( 67)
     &                      +(-1.000)*r( 67)+( 4.000)*r( 68)
     &                      +( 1.000)*r( 77)+( 1.000)*r( 87)
     &                      +( 1.000)*r( 98)+( 1.000)*r(110)
          JAC(kRO2N,kCCO3)= +( 1.000)*r( 77)
          JAC(kRO2N,kRCO3)= +( 1.000)*r( 87)
          JAC(kRO2N,kMCO3)= +( 1.000)*r(110)
          JAC(kRO2N,kBZCO)= +( 1.000)*r( 98)
          JAC(kRO2N,kCXO2)= +( 1.000)*r( 64)
          JAC(kCCO3,kOH  )= +(-1.000)*r(130)+(-1.000)*r(136)
     &                      +(-0.492)*r(138)+(-1.000)*r(150)
     &                      +(-0.675)*r(166)+(-0.029)*r(174)
     &                      +(-1.000)*r(180)+(-1.000)*r(182)
     &                      +(-0.011)*r(200)
          JAC(kCCO3,kHO2 )= +( 1.000)*r( 72)
          JAC(kCCO3,kRO2R)= +( 1.000)*r( 75)
          JAC(kCCO3,kR2O2)= +( 1.000)*r( 76)+(-1.000)*r( 76)
          JAC(kCCO3,kRO2N)= +( 1.000)*r( 77)
          JAC(kCCO3,kCCO3)= +( 1.000)*r( 69)+( 1.000)*r( 71)
     &                      +( 1.000)*r( 72)+( 1.000)*r( 73)
     &                      +( 1.000)*r( 74)+( 1.000)*r( 75)
     &                      +( 1.000)*r( 76)+(-1.000)*r( 76)
     &                      +( 1.000)*r( 77)+( 4.000)*r( 78)
     &                      +( 1.000)*r( 88)+( 1.000)*r( 99)
     &                      +( 1.000)*r(111)+(-1.000)*r(111)
          JAC(kCCO3,kRCO3)= +( 1.000)*r( 88)+(-1.000)*r(112)
          JAC(kCCO3,kMCO3)= +(-1.000)*r(104)+(-1.000)*r(106)
     &                      +( 1.000)*r(111)+(-1.000)*r(111)
     &                      +(-1.000)*r(112)+(-1.000)*r(113)
     &                      +(-4.000)*r(114)
          JAC(kCCO3,kBZCO)= +( 1.000)*r( 99)+(-1.000)*r(113)
          JAC(kCCO3,kCXO2)= +( 1.000)*r( 74)
          JAC(kRCO3,kOH  )= +(-0.965)*r(133)+(-0.096)*r(138)
     &                      +(-0.370)*r(147)+(-0.049)*r(174)
          JAC(kRCO3,kHO2 )= +( 1.000)*r( 82)
          JAC(kRCO3,kRO2R)= +( 1.000)*r( 85)
          JAC(kRCO3,kR2O2)= +( 1.000)*r( 86)+(-1.000)*r( 86)
          JAC(kRCO3,kRO2N)= +( 1.000)*r( 87)
          JAC(kRCO3,kCCO3)= +( 1.000)*r( 88)
          JAC(kRCO3,kRCO3)= +( 1.000)*r( 79)+( 1.000)*r( 81)
     &                      +( 1.000)*r( 82)+( 1.000)*r( 83)
     &                      +( 1.000)*r( 84)+( 1.000)*r( 85)
     &                      +( 1.000)*r( 86)+(-1.000)*r( 86)
     &                      +( 1.000)*r( 87)+( 1.000)*r( 88)
     &                      +( 4.000)*r( 89)+( 1.000)*r(100)
     &                      +( 1.000)*r(112)
          JAC(kRCO3,kMCO3)= +( 1.000)*r(112)
          JAC(kRCO3,kBZCO)= +( 1.000)*r(100)
          JAC(kRCO3,kCXO2)= +( 1.000)*r( 84)
          JAC(kMCO3,kOH  )= +(-0.500)*r(161)+(-0.289)*r(170)
          JAC(kMCO3,kHO2 )= +( 1.000)*r(105)
          JAC(kMCO3,kRO2R)= +( 1.000)*r(108)
          JAC(kMCO3,kR2O2)= +( 1.000)*r(109)+(-1.000)*r(109)
          JAC(kMCO3,kRO2N)= +( 1.000)*r(110)
          JAC(kMCO3,kCCO3)= +( 1.000)*r(111)
          JAC(kMCO3,kRCO3)= +( 1.000)*r(112)
          JAC(kMCO3,kMCO3)= +( 1.000)*r(102)+( 1.000)*r(104)
     &                      +( 1.000)*r(105)+( 1.000)*r(106)
     &                      +( 1.000)*r(107)+( 1.000)*r(108)
     &                      +( 1.000)*r(109)+(-1.000)*r(109)
     &                      +( 1.000)*r(110)+( 1.000)*r(111)
     &                      +( 1.000)*r(112)+( 1.000)*r(113)
     &                      +( 4.000)*r(114)
          JAC(kMCO3,kBZCO)= +( 1.000)*r(113)
          JAC(kMCO3,kCXO2)= +( 1.000)*r(107)
          JAC(kBZCO,kOH  )= +(-1.000)*r(158)
          JAC(kBZCO,kHO2 )= +( 1.000)*r( 93)
          JAC(kBZCO,kRO2R)= +( 1.000)*r( 96)
          JAC(kBZCO,kR2O2)= +( 1.000)*r( 97)+(-1.000)*r( 97)
          JAC(kBZCO,kRO2N)= +( 1.000)*r( 98)
          JAC(kBZCO,kCCO3)= +( 1.000)*r( 99)
          JAC(kBZCO,kRCO3)= +( 1.000)*r(100)
          JAC(kBZCO,kMCO3)= +( 1.000)*r(113)
          JAC(kBZCO,kBZCO)= +( 1.000)*r( 90)+( 1.000)*r( 92)
     &                      +( 1.000)*r( 93)+( 1.000)*r( 94)
     &                      +( 1.000)*r( 95)+( 1.000)*r( 96)
     &                      +( 1.000)*r( 97)+(-1.000)*r( 97)
     &                      +( 1.000)*r( 98)+( 1.000)*r( 99)
     &                      +( 1.000)*r(100)+( 4.000)*r(101)
     &                      +( 1.000)*r(113)
          JAC(kBZCO,kCXO2)= +( 1.000)*r( 95)
          JAC(kCXO2,kOH  )= +(-0.650)*r(141)+(-1.000)*r(184)
     &                      +(-0.011)*r(200)
          JAC(kCXO2,kHO2 )= +( 1.000)*r( 47)
          JAC(kCXO2,kRO2R)= +( 1.000)*r( 54)
          JAC(kCXO2,kR2O2)= +( 1.000)*r( 59)+(-1.000)*r( 59)
          JAC(kCXO2,kRO2N)= +( 1.000)*r( 64)
          JAC(kCXO2,kCCO3)= +(-1.000)*r( 71)+(-1.000)*r( 73)
     &                      +( 1.000)*r( 74)+(-4.000)*r( 78)
     &                      +(-1.000)*r( 88)+(-1.000)*r( 99)
     &                      +(-1.000)*r(111)
          JAC(kCXO2,kRCO3)= +( 1.000)*r( 84)+(-1.000)*r( 88)
          JAC(kCXO2,kMCO3)= +( 1.000)*r(107)+(-1.000)*r(111)
          JAC(kCXO2,kBZCO)= +( 1.000)*r( 95)+(-1.000)*r( 99)
          JAC(kCXO2,kCXO2)= +( 1.000)*r( 46)+( 1.000)*r( 47)
     &                      +( 1.000)*r( 48)+( 4.000)*r( 49)
     &                      +( 4.000)*r( 50)+( 1.000)*r( 54)
     &                      +( 1.000)*r( 59)+(-1.000)*r( 59)
     &                      +( 1.000)*r( 64)+( 1.000)*r( 74)
     &                      +( 1.000)*r( 84)+( 1.000)*r( 95)
     &                      +( 1.000)*r(107)
          JAC(kCXO2,kTBUO)= +(-1.000)*r(116)
          JAC(kHCO3,kHO2 )= +(-1.000)*r(126)
          JAC(kHCO3,kHCO3)= +( 1.000)*r(127)+( 1.000)*r(128)
          JAC(kTBUO,kOH  )= +(-0.236)*r(199)
          JAC(kTBUO,kTBUO)= +( 1.000)*r(115)+( 1.000)*r(116)
          JAC(kBZO ,kOH  )= +(-0.240)*r(153)+(-0.240)*r(155)
          JAC(kBZO ,kHO2 )= +( 1.000)*r(118)
          JAC(kBZO ,kCCO3)= +(-1.000)*r( 99)
          JAC(kBZO ,kRCO3)= +(-1.000)*r(100)
          JAC(kBZO ,kMCO3)= +(-1.000)*r(113)
          JAC(kBZO ,kBZCO)= +(-1.000)*r( 92)+(-1.000)*r( 94)
     &                      +(-1.000)*r( 99)+(-1.000)*r(100)
     &                      +(-4.000)*r(101)+(-1.000)*r(113)
          JAC(kBZO ,kBZO )= +( 1.000)*r(117)+( 1.000)*r(118)
     &                      +( 1.000)*r(119)
          JAC(kBZNO,kHO2 )= +( 1.000)*r(121)
          JAC(kBZNO,kBZNO)= +( 1.000)*r(120)+( 1.000)*r(121)
     &                      +( 1.000)*r(122)
c
c  complete the rates and Jacobian
c
      do i=nstrt,nend
        crold(i) = cncrad(i)
        rate(i) = gain(i) - loss(i)
        do j=nstrt,nend
          jac(i,j) = jac(i,j)/cncrad(j)
        enddo
      enddo
c
c
c  Jacobian is modified due to substitution of HCO3
c
c
c    OH   HO2  RO2R  R2O2  RO2N  CCO3  RCO3  MCO3  BZCO  CXO2  HCO3  TBUO   BZO  BZNO
c
          JAC(kHO2 ,kHO2 ) = JAC(kHO2 ,kHO2 )
     &         - JAC(kHO2 ,kHCO3)*JAC(kHCO3,kHO2 )/JAC(kHCO3,kHCO3)
c
c  end of substitution of HCO3
c
c
c  Jacobian is modified due to substitution of TBUO
c
c
c    OH   HO2  RO2R  R2O2  RO2N  CCO3  RCO3  MCO3  BZCO  CXO2  HCO3  TBUO   BZO  BZNO
c
          JAC(kCXO2,kOH  ) = JAC(kCXO2,kOH  )
     &         - JAC(kCXO2,kTBUO)*JAC(kTBUO,kOH  )/JAC(kTBUO,kTBUO)
c
c  end of substitution of TBUO
c
c
c  Jacobian is modified due to substitution of BZO
c
c
c    OH   HO2  RO2R  R2O2  RO2N  CCO3  RCO3  MCO3  BZCO  CXO2  HCO3  TBUO   BZO  BZNO
c
          JAC(kHO2 ,kOH  ) = JAC(kHO2 ,kOH  )
     &         - JAC(kHO2 ,kBZO )*JAC(kBZO ,kOH  )/JAC(kBZO ,kBZO )
          JAC(kHO2 ,kHO2 ) = JAC(kHO2 ,kHO2 )
     &         - JAC(kHO2 ,kBZO )*JAC(kBZO ,kHO2 )/JAC(kBZO ,kBZO )
          JAC(kHO2 ,kCCO3) = JAC(kHO2 ,kCCO3)
     &         - JAC(kHO2 ,kBZO )*JAC(kBZO ,kCCO3)/JAC(kBZO ,kBZO )
          JAC(kHO2 ,kRCO3) = JAC(kHO2 ,kRCO3)
     &         - JAC(kHO2 ,kBZO )*JAC(kBZO ,kRCO3)/JAC(kBZO ,kBZO )
          JAC(kHO2 ,kMCO3) = JAC(kHO2 ,kMCO3)
     &         - JAC(kHO2 ,kBZO )*JAC(kBZO ,kMCO3)/JAC(kBZO ,kBZO )
          JAC(kHO2 ,kBZCO) = JAC(kHO2 ,kBZCO)
     &         - JAC(kHO2 ,kBZO )*JAC(kBZO ,kBZCO)/JAC(kBZO ,kBZO )
c
c  end of substitution of BZO
c
c
c  Jacobian is modified due to substitution of BZNO
c
c
c    OH   HO2  RO2R  R2O2  RO2N  CCO3  RCO3  MCO3  BZCO  CXO2  HCO3  TBUO   BZO  BZNO
c
          JAC(kHO2 ,kHO2 ) = JAC(kHO2 ,kHO2 )
     &         - JAC(kHO2 ,kBZNO)*JAC(kBZNO,kHO2 )/JAC(kBZNO,kBZNO)
c
c  end of substitution of BZNO
c
c
c  solve the matrix
c
      if (errbig.gt.500.0 .or. lsafe) then
        thresh = 1.0e+15
      else
        thresh = 1.0e-15
      endif
      lusecp = .false.
c
      nendcp=nend-4
      do j=nstrt,nendcp
        if (cncrad(j).le.thresh) then
          do i=nstrt,nendcp
            jac(i,j) = 0.
            jac(j,i) = 0.
          enddo
          jac(j,j) = 1.
      else
        lusecp = .true.
        endif
      enddo
c
      if (lusecp) then
        call cpivot(nstrt,nendcp,MXRADCL,jac,rate,ierr)
        if (ierr.ne.0) goto 900
      endif
c
c  update radical concentrations
c
      errbig = 0.
      weight = 1.
      rlim = 100.0
      if (kount.gt.10) then
        rlim = 10.0
        weight = 0.7
      elseif (kount.gt.50) then
        rlim = 3.0
        weight = 0.5
      endif
      do l=nstrt,nendcp
        if (cncrad(l) .le. thresh) then ! solve independantly
          cncrad(l) = crold(l)*(1.-weight) +
     &                         weight*gain(l)/loss(l)*crold(l)
          cncrad(l) = amax1(cncrad(l),crold(l)/rlim)
          cncrad(l) = amin1(cncrad(l),crold(l)*rlim)
          err = abs((cncrad(l)/crold(l))-1.0)
          err = amin1(err,thresh)
        else                            ! part of coupled solution
          rate(l) = amin1(rate(l), rlim*cncrad(l))
          rate(l) = amax1(rate(l), -cncrad(l)*(1.0-(1.0/rlim)))
          cncrad(l) = cncrad(l) + rate(l)*weight
          err = abs((cncrad(l)/crold(l))-1.0)
        endif
        cncrad(l) = amax1(cncrad(l), bdlrad)
c        cncrad(l) = amin1(cncrad(l), 10.0)
        errbig = amax1(errbig,err)
      enddo
c
      if (kount.eq.100) then
        lsafe = .true.
        write(iout,'(a)') 'WARNING:'
        write(iout,'(a)') 'Monitor slow convergence in RADSLVR5'
        write(iout,'(a3,30a9)')
     &     ' ', (nmrad(l),l=nstrt,nendcp),'rel err'
        do l=nstrt,nendcp
          cncrad(l) = 1.0e-12
        enddo
      endif
      if (kount.ge.100)
     &   write(iout,'(i3,1p30e9.2)')
     &     kount, (cncrad(l),l=nstrt,nendcp), errbig
c
      r( 21) = rk( 21)*cncrad(kOH)*conc(kNO)
      r( 24) = rk( 24)*cncrad(kOH)*conc(kHONO)
      r( 25) = rk( 25)*cncrad(kOH)*conc(kNO2)
      r( 26) = rk( 26)*cncrad(kOH)*cncrad(kNO3)
      r( 27) = rk( 27)*cncrad(kOH)*conc(kHNO3)
      r( 29) = rk( 29)*cncrad(kOH)*conc(kCO)
      r( 30) = rk( 30)*cncrad(kOH)*conc(kO3)
      r( 31) = rk( 31)*cncrad(kHO2)*conc(kNO)
      r( 32) = rk( 32)*cncrad(kHO2)*conc(kNO2)
      r( 35) = rk( 35)*conc(kHNO4)*cncrad(kOH)
      r( 36) = rk( 36)*cncrad(kHO2)*conc(kO3)
      r( 37) = rk( 37)*cncrad(kHO2)*cncrad(kHO2)
      r( 38) = rk( 38)*cncrad(kHO2)*cncrad(kHO2)*H2O
      r( 39) = rk( 39)*cncrad(kNO3)*cncrad(kHO2)
      r( 42) = rk( 42)*conc(kHO2H)*cncrad(kOH)
      r( 43) = rk( 43)*cncrad(kOH)*cncrad(kHO2)
      r( 44) = rk( 44)*cncrad(kOH)*conc(kSO2)
      r( 45) = rk( 45)*cncrad(kOH)*H2
      r( 46) = rk( 46)*cncrad(kCXO2)*conc(kNO)
      r( 47) = rk( 47)*cncrad(kCXO2)*cncrad(kHO2)
      r( 48) = rk( 48)*cncrad(kCXO2)*cncrad(kNO3)
      r( 49) = rk( 49)*cncrad(kCXO2)*cncrad(kCXO2)
      r( 50) = rk( 50)*cncrad(kCXO2)*cncrad(kCXO2)
      r( 51) = rk( 51)*cncrad(kRO2R)*conc(kNO)
      r( 52) = rk( 52)*cncrad(kRO2R)*cncrad(kHO2)
      r( 53) = rk( 53)*cncrad(kRO2R)*cncrad(kNO3)
      r( 54) = rk( 54)*cncrad(kRO2R)*cncrad(kCXO2)
      r( 55) = rk( 55)*cncrad(kRO2R)*cncrad(kRO2R)
      r( 56) = rk( 56)*cncrad(kR2O2)*conc(kNO)
      r( 57) = rk( 57)*cncrad(kR2O2)*cncrad(kHO2)
      r( 58) = rk( 58)*cncrad(kR2O2)*cncrad(kNO3)
      r( 59) = rk( 59)*cncrad(kR2O2)*cncrad(kCXO2)
      r( 60) = rk( 60)*cncrad(kR2O2)*cncrad(kRO2R)
      r( 61) = rk( 61)*cncrad(kR2O2)*cncrad(kR2O2)
      r( 62) = rk( 62)*cncrad(kRO2N)*conc(kNO)
      r( 63) = rk( 63)*cncrad(kRO2N)*cncrad(kHO2)
      r( 64) = rk( 64)*cncrad(kRO2N)*cncrad(kCXO2)
      r( 65) = rk( 65)*cncrad(kRO2N)*cncrad(kNO3)
      r( 66) = rk( 66)*cncrad(kRO2N)*cncrad(kRO2R)
      r( 67) = rk( 67)*cncrad(kRO2N)*cncrad(kR2O2)
      r( 68) = rk( 68)*cncrad(kRO2N)*cncrad(kRO2N)
      r( 69) = rk( 69)*cncrad(kCCO3)*conc(kNO2)
      r( 71) = rk( 71)*cncrad(kCCO3)*conc(kNO)
      r( 72) = rk( 72)*cncrad(kCCO3)*cncrad(kHO2)
      r( 73) = rk( 73)*cncrad(kCCO3)*cncrad(kNO3)
      r( 74) = rk( 74)*cncrad(kCCO3)*cncrad(kCXO2)
      r( 75) = rk( 75)*cncrad(kCCO3)*cncrad(kRO2R)
      r( 76) = rk( 76)*cncrad(kCCO3)*cncrad(kR2O2)
      r( 77) = rk( 77)*cncrad(kCCO3)*cncrad(kRO2N)
      r( 78) = rk( 78)*cncrad(kCCO3)*cncrad(kCCO3)
      r( 79) = rk( 79)*cncrad(kRCO3)*conc(kNO2)
      r( 81) = rk( 81)*cncrad(kRCO3)*conc(kNO)
      r( 82) = rk( 82)*cncrad(kRCO3)*cncrad(kHO2)
      r( 83) = rk( 83)*cncrad(kRCO3)*cncrad(kNO3)
      r( 84) = rk( 84)*cncrad(kRCO3)*cncrad(kCXO2)
      r( 85) = rk( 85)*cncrad(kRCO3)*cncrad(kRO2R)
      r( 86) = rk( 86)*cncrad(kRCO3)*cncrad(kR2O2)
      r( 87) = rk( 87)*cncrad(kRCO3)*cncrad(kRO2N)
      r( 88) = rk( 88)*cncrad(kRCO3)*cncrad(kCCO3)
      r( 89) = rk( 89)*cncrad(kRCO3)*cncrad(kRCO3)
      r( 90) = rk( 90)*cncrad(kBZCO)*conc(kNO2)
      r( 92) = rk( 92)*cncrad(kBZCO)*conc(kNO)
      r( 93) = rk( 93)*cncrad(kBZCO)*cncrad(kHO2)
      r( 94) = rk( 94)*cncrad(kBZCO)*cncrad(kNO3)
      r( 95) = rk( 95)*cncrad(kBZCO)*cncrad(kCXO2)
      r( 96) = rk( 96)*cncrad(kBZCO)*cncrad(kRO2R)
      r( 97) = rk( 97)*cncrad(kBZCO)*cncrad(kR2O2)
      r( 98) = rk( 98)*cncrad(kBZCO)*cncrad(kRO2N)
      r( 99) = rk( 99)*cncrad(kBZCO)*cncrad(kCCO3)
      r(100) = rk(100)*cncrad(kBZCO)*cncrad(kRCO3)
      r(101) = rk(101)*cncrad(kBZCO)*cncrad(kBZCO)
      r(102) = rk(102)*cncrad(kMCO3)*conc(kNO2)
      r(104) = rk(104)*cncrad(kMCO3)*conc(kNO)
      r(105) = rk(105)*cncrad(kMCO3)*cncrad(kHO2)
      r(106) = rk(106)*cncrad(kMCO3)*cncrad(kNO3)
      r(107) = rk(107)*cncrad(kMCO3)*cncrad(kCXO2)
      r(108) = rk(108)*cncrad(kMCO3)*cncrad(kRO2R)
      r(109) = rk(109)*cncrad(kMCO3)*cncrad(kR2O2)
      r(110) = rk(110)*cncrad(kMCO3)*cncrad(kRO2N)
      r(111) = rk(111)*cncrad(kMCO3)*cncrad(kCCO3)
      r(112) = rk(112)*cncrad(kMCO3)*cncrad(kRCO3)
      r(113) = rk(113)*cncrad(kMCO3)*cncrad(kBZCO)
      r(114) = rk(114)*cncrad(kMCO3)*cncrad(kMCO3)
      r(115) = rk(115)*cncrad(kTBUO)*conc(kNO2)
      r(116) = rk(116)*cncrad(kTBUO)
      r(117) = rk(117)*cncrad(kBZO)*conc(kNO2)
      r(118) = rk(118)*cncrad(kBZO)*cncrad(kHO2)
      r(119) = rk(119)*cncrad(kBZO)
      r(120) = rk(120)*cncrad(kBZNO)*conc(kNO2)
      r(121) = rk(121)*cncrad(kBZNO)*cncrad(kHO2)
      r(122) = rk(122)*cncrad(kBZNO)
      r(125) = rk(125)*conc(kHCHO)*cncrad(kOH)
      r(126) = rk(126)*conc(kHCHO)*cncrad(kHO2)
      r(127) = rk(127)*cncrad(kHCO3)
      r(128) = rk(128)*cncrad(kHCO3)*conc(kNO)
      r(130) = rk(130)*conc(kCCHO)*cncrad(kOH)
      r(133) = rk(133)*conc(kRCHO)*cncrad(kOH)
      r(136) = rk(136)*conc(kACET)*cncrad(kOH)
      r(138) = rk(138)*conc(kMEK)*cncrad(kOH)
      r(140) = rk(140)*conc(kMEOH)*cncrad(kOH)
      r(141) = rk(141)*conc(kCOOH)*cncrad(kOH)
      r(143) = rk(143)*conc(kROOH)*cncrad(kOH)
      r(147) = rk(147)*conc(kGLY)*cncrad(kOH)
      r(150) = rk(150)*conc(kMGLY)*cncrad(kOH)
      r(153) = rk(153)*conc(kPHEN)*cncrad(kOH)
      r(155) = rk(155)*conc(kCRES)*cncrad(kOH)
      r(158) = rk(158)*conc(kBALD)*cncrad(kOH)
      r(161) = rk(161)*conc(kMETH)*cncrad(kOH)
      r(166) = rk(166)*conc(kMVK)*cncrad(kOH)
      r(170) = rk(170)*conc(kISPD)*cncrad(kOH)
      r(174) = rk(174)*conc(kPROD)*cncrad(kOH)
      r(176) = rk(176)*conc(kRNO3)*cncrad(kOH)
      r(178) = rk(178)*conc(kDCB1)*cncrad(kOH)
      r(180) = rk(180)*conc(kDCB2)*cncrad(kOH)
      r(182) = rk(182)*conc(kDCB3)*cncrad(kOH)
      r(184) = rk(184)*CH4*cncrad(kOH)
      r(185) = rk(185)*conc(kETHE)*cncrad(kOH)
      r(189) = rk(189)*conc(kISOP)*cncrad(kOH)
      r(193) = rk(193)*conc(kTERP)*cncrad(kOH)
      r(197) = rk(197)*conc(kALK1)*cncrad(kOH)
      r(198) = rk(198)*conc(kALK2)*cncrad(kOH)
      r(199) = rk(199)*conc(kALK3)*cncrad(kOH)
      r(200) = rk(200)*conc(kALK4)*cncrad(kOH)
      r(201) = rk(201)*conc(kALK5)*cncrad(kOH)
      r(202) = rk(202)*conc(kARO1)*cncrad(kOH)
      r(203) = rk(203)*conc(kARO2)*cncrad(kOH)
      r(204) = rk(204)*conc(kOLE1)*cncrad(kOH)
      r(208) = rk(208)*conc(kOLE2)*cncrad(kOH)
      r(212) = rk(212)*cncrad(kOH)*conc(kCPO1)
      r(213) = rk(213)*cncrad(kOH)*conc(kCPO2)
      r(214) = rk(214)*cncrad(kOH)*conc(kCPO3)
      r(215) = rk(215)*cncrad(kOH)*conc(kCPO4)
      r(216) = rk(216)*cncrad(kOH)*conc(kCPO5)
      r(217) = rk(217)*cncrad(kOH)*conc(kCPO6)
      r(218) = rk(218)*cncrad(kOH)*conc(kCPO7)
      r(219) = rk(219)*cncrad(kOH)*conc(kCPO8)
      r(220) = rk(220)*cncrad(kOH)*conc(kCPO9)
      r(221) = rk(221)*cncrad(kOH)*conc(kCOO1)
      r(222) = rk(222)*cncrad(kOH)*conc(kCOO2)
      r(223) = rk(223)*cncrad(kOH)*conc(kCOO3)
      r(224) = rk(224)*cncrad(kOH)*conc(kCOO4)
      r(225) = rk(225)*cncrad(kOH)*conc(kCOO5)
      r(226) = rk(226)*cncrad(kOH)*conc(kCOO6)
      r(227) = rk(227)*cncrad(kOH)*conc(kCOO7)
      r(228) = rk(228)*cncrad(kOH)*conc(kCOO8)
      r(229) = rk(229)*cncrad(kOH)*conc(kCOO9)
      r(230) = rk(230)*cncrad(kOH)*conc(kCBS1)
      r(231) = rk(231)*cncrad(kOH)*conc(kCBS2)
      r(232) = rk(232)*cncrad(kOH)*conc(kCBS3)
      r(233) = rk(233)*cncrad(kOH)*conc(kCBS4)
      r(234) = rk(234)*cncrad(kOH)*conc(kCAS1)
      r(235) = rk(235)*cncrad(kOH)*conc(kCAS2)
      r(236) = rk(236)*cncrad(kOH)*conc(kCAS3)
      r(237) = rk(237)*cncrad(kOH)*conc(kCAS4)
 
      if (errbig.gt.tol) goto 14
c
c  new group of radicals
c
c  HCO3
c
       Loss(kHCO3 )= +( 1.000)*r(127)+( 1.000)*r(128)
       Gain(kHCO3 )= +( 1.000)*r(126)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kHCO3) = gain(kHCO3)/loss(kHCO3)*cncrad(kHCO3)
      cncrad(kHCO3) = amax1(bdlrad, cncrad(kHCO3))
      r(127) = rk(127)*cncrad(kHCO3)
      r(128) = rk(128)*cncrad(kHCO3)*conc(kNO)
c
c  new group of radicals
c
c  TBUO
c
       Loss(kTBUO )= +( 1.000)*r(115)+( 1.000)*r(116)
       Gain(kTBUO )= +( 0.236)*r(199)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kTBUO) = gain(kTBUO)/loss(kTBUO)*cncrad(kTBUO)
      cncrad(kTBUO) = amax1(bdlrad, cncrad(kTBUO))
      r(115) = rk(115)*cncrad(kTBUO)*conc(kNO2)
      r(116) = rk(116)*cncrad(kTBUO)
c
c  new group of radicals
c
c   BZO
c
        Loss(kBZO  )= +( 1.000)*r(117)+( 1.000)*r(118)+( 1.000)*r(119)
        Gain(kBZO  )= +( 1.000)*r( 92)+( 1.000)*r( 94)+( 1.000)*r( 99)
     &                +( 1.000)*r(100)+( 2.000)*r(101)+( 1.000)*r(113)
     &                +( 0.240)*r(153)+( 1.000)*r(154)+( 0.240)*r(155)
     &                +( 1.000)*r(156)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kBZO) = gain(kBZO)/loss(kBZO)*cncrad(kBZO)
      cncrad(kBZO) = amax1(bdlrad, cncrad(kBZO))
      r(117) = rk(117)*cncrad(kBZO)*conc(kNO2)
      r(118) = rk(118)*cncrad(kBZO)*cncrad(kHO2)
      r(119) = rk(119)*cncrad(kBZO)
c
c  new group of radicals
c
c  BZNO
c
        Loss(kBZNO )= +( 1.000)*r(120)+( 1.000)*r(121)+( 1.000)*r(122)
        Gain(kBZNO )= +( 1.000)*r(157)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kBZNO) = gain(kBZNO)/loss(kBZNO)*cncrad(kBZNO)
      cncrad(kBZNO) = amax1(bdlrad, cncrad(kBZNO))
      r(120) = rk(120)*cncrad(kBZNO)*conc(kNO2)
      r(121) = rk(121)*cncrad(kBZNO)*cncrad(kHO2)
      r(122) = rk(122)*cncrad(kBZNO)
c
      return
c
 900  continue
      write(iout,'(//,A,//)') 'ERROR in RADSLVR5:'
      write(iout,*) 'Zero determinant in CPIVOT at ', ierr
      write(iout,*) 'igrd,i, j, k = ', igrdchm,ichm,jchm,kchm
      write(iout,*) 'LDARK is set ', ldark
      do l=1,ngas
        write(iout,'(i3,2x,a7,1pe10.3)') l,spname(l),conc(l)
      enddo
      write(iout,*) 'The radicals are: '
      do l= 1 , nrad
        write(iout,'(i3,2x,a7,1pe10.3)')
     &        l,nmrad(l),cncrad(l)
      enddo
      write(iout,*) 'Currently solving ', nstrt, ' to ',nendcp
      do l=nstrt,nendcp
        write(iout,'(i3,2x,a7,1p2e10.3)')
     &        l,nmrad(l),cncrad(l),abs(rate(l))/cncrad(l)
      enddo
      call camxerr()
      end
