      subroutine rxnrate3(H2O,M,O2,CH4,H2,cncrad,conc,r)
c
c-----CAMx v4.02 030709
c
c     RXNRATE computes fluxes for each reaction
c
c     Copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Routines Called:
c        none
c
c     Called by:
c        RADINIT3
c        TRAP
c
      include "camx.prm"
      include "chmstry.com"
c
      real M
      real cncrad(MXRADCL),conc(MXSPEC+1),r(MXRXN)
c
      r(  1) = rk(  1)*conc(kNO2)
      r(  2) = rk(  2)*cncrad(kO)
      r(  3) = rk(  3)*conc(kO3)*conc(kNO)
      r(  4) = rk(  4)*cncrad(kO)*conc(kNO2)
      r(  5) = rk(  5)*cncrad(kO)*conc(kNO2)
      r(  6) = rk(  6)*cncrad(kO)*conc(kNO)
      r(  7) = rk(  7)*conc(kNO2)*conc(kO3)
      r(  8) = rk(  8)*conc(kO3)
      r(  9) = rk(  9)*conc(kO3)
      r( 10) = rk( 10)*cncrad(kO1D)
      r( 11) = rk( 11)*cncrad(kO1D)*H2O
      r( 12) = rk( 12)*conc(kO3)*cncrad(kOH)
      r( 13) = rk( 13)*conc(kO3)*cncrad(kHO2)
      r( 14) = rk( 14)*cncrad(kNO3)
      r( 15) = rk( 15)*cncrad(kNO3)*conc(kNO)
      r( 16) = rk( 16)*cncrad(kNO3)*conc(kNO2)
      r( 17) = rk( 17)*cncrad(kNO3)*conc(kNO2)
      r( 18) = rk( 18)*cncrad(kN2O5)*H2O
      r( 19) = rk( 19)*cncrad(kN2O5)
      r( 20) = rk( 20)*conc(kNO)*conc(kNO)
      r( 21) = rk( 21)*conc(kNO)*conc(kNO2)*H2O
      r( 22) = rk( 22)*conc(kNO)*cncrad(kOH)
      r( 23) = rk( 23)*conc(kHONO)
      r( 24) = rk( 24)*cncrad(kOH)*conc(kHONO)
      r( 25) = rk( 25)*conc(kHONO)*conc(kHONO)
      r( 26) = rk( 26)*conc(kNO2)*cncrad(kOH)
      r( 27) = rk( 27)*cncrad(kOH)*conc(kHNO3)
      r( 28) = rk( 28)*cncrad(kHO2)*conc(kNO)
      r( 29) = rk( 29)*cncrad(kHO2)*conc(kNO2)
      r( 30) = rk( 30)*conc(kPNA)
      r( 31) = rk( 31)*cncrad(kOH)*conc(kPNA)
      r( 32) = rk( 32)*cncrad(kHO2)*cncrad(kHO2)
      r( 33) = rk( 33)*cncrad(kHO2)*cncrad(kHO2)*H2O
      r( 34) = rk( 34)*conc(kH2O2)
      r( 35) = rk( 35)*cncrad(kOH)*conc(kH2O2)
      r( 36) = rk( 36)*cncrad(kOH)*conc(kCO)
      r( 37) = rk( 37)*conc(kFORM)*cncrad(kOH)
      r( 38) = rk( 38)*conc(kFORM)
      r( 39) = rk( 39)*conc(kFORM)
      r( 40) = rk( 40)*conc(kFORM)*cncrad(kO)
      r( 41) = rk( 41)*conc(kFORM)*cncrad(kNO3)
      r( 42) = rk( 42)*conc(kALD2)*cncrad(kO)
      r( 43) = rk( 43)*conc(kALD2)*cncrad(kOH)
      r( 44) = rk( 44)*conc(kALD2)*cncrad(kNO3)
      r( 45) = rk( 45)*conc(kALD2)
      r( 46) = rk( 46)*cncrad(kC2O3)*conc(kNO)
      r( 47) = rk( 47)*cncrad(kC2O3)*conc(kNO2)
      r( 48) = rk( 48)*conc(kPAN)
      r( 49) = rk( 49)*cncrad(kC2O3)*cncrad(kC2O3)
      r( 50) = rk( 50)*cncrad(kC2O3)*cncrad(kHO2)
      r( 51) = rk( 51)*cncrad(kOH)
      r( 52) = rk( 52)*conc(kPAR)*cncrad(kOH)
      r( 53) = rk( 53)*cncrad(kROR)
      r( 54) = rk( 54)*cncrad(kROR)
      r( 55) = rk( 55)*cncrad(kROR)*conc(kNO2)
      r( 56) = rk( 56)*cncrad(kO)*conc(kOLE)
      r( 57) = rk( 57)*cncrad(kOH)*conc(kOLE)
      r( 58) = rk( 58)*conc(kO3)*conc(kOLE)
      r( 59) = rk( 59)*cncrad(kNO3)*conc(kOLE)
      r( 60) = rk( 60)*cncrad(kO)*conc(kETH)
      r( 61) = rk( 61)*cncrad(kOH)*conc(kETH)
      r( 62) = rk( 62)*conc(kO3)*conc(kETH)
      r( 63) = rk( 63)*conc(kTOL)*cncrad(kOH)
      r( 64) = rk( 64)*cncrad(kTO2)*conc(kNO)
      r( 65) = rk( 65)*cncrad(kTO2)
      r( 66) = rk( 66)*cncrad(kOH)*conc(kCRES)
      r( 67) = rk( 67)*conc(kCRES)*cncrad(kNO3)
      r( 68) = rk( 68)*cncrad(kCRO)*conc(kNO2)
      r( 69) = rk( 69)*conc(kOPEN)
      r( 70) = rk( 70)*conc(kOPEN)*cncrad(kOH)
      r( 71) = rk( 71)*conc(kOPEN)*conc(kO3)
      r( 72) = rk( 72)*cncrad(kOH)*conc(kXYL)
      r( 73) = rk( 73)*cncrad(kOH)*conc(kMGLY)
      r( 74) = rk( 74)*conc(kMGLY)
      r( 75) = rk( 75)*cncrad(kO)*conc(kISOP)
      r( 76) = rk( 76)*cncrad(kOH)*conc(kISOP)
      r( 77) = rk( 77)*conc(kO3)*conc(kISOP)
      r( 78) = rk( 78)*cncrad(kNO3)*conc(kISOP)
      r( 79) = rk( 79)*cncrad(kXO2)*conc(kNO)
      r( 80) = rk( 80)*cncrad(kXO2)*cncrad(kXO2)
      r( 81) = rk( 81)*cncrad(kXO2N)*conc(kNO)
      r( 82) = rk( 82)*conc(kSO2)*cncrad(kOH)
      r( 83) = rk( 83)*conc(kSO2)
      r( 84) = rk( 84)*conc(kMEOH)*cncrad(kOH)
      r( 85) = rk( 85)*conc(kETOH)*cncrad(kOH)
      r( 86) = rk( 86)*cncrad(kXO2)*cncrad(kHO2)
      r( 87) = rk( 87)*cncrad(kXO2N)*cncrad(kHO2)
      r( 88) = rk( 88)*cncrad(kXO2N)*cncrad(kXO2N)
      r( 89) = rk( 89)*cncrad(kXO2)*cncrad(kXO2N)
      r( 90) = rk( 90)*cncrad(kOH)*cncrad(kHO2)
      r( 91) = rk( 91)*cncrad(kCRO)
      r( 92) = rk( 92)*cncrad(kOH)*conc(kISPD)
      r( 93) = rk( 93)*conc(kO3)*conc(kISPD)
      r( 94) = rk( 94)*cncrad(kNO3)*conc(kISPD)
      r( 95) = rk( 95)*conc(kISPD)
      r( 96) = rk( 96)*conc(kNO2)*conc(kISOP)
c
      return
      end
