      subroutine ierxn4(y,ny,r,rk,nk)
c
c-----CAMx v4.02 030709
c
c     IERXN4 computes IEH solver fluxes for each reaction
c
c     Copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Routines Called:
c        none
c
c     Called by:
c        IEHSOLV
c        IERATE4
c        IEJAC4
c
      include "camx.prm"
      include "iehchem.com"
c
      real H2O, M, O2, CH4, H2
      real y(ny+6),r(nk),rk(nk)
c
c --- Entry point
c
      H2O = y(ny+2)
      M   = y(ny+3)
      O2  = y(ny+4)
      CH4 = y(ny+5)
      H2  = y(ny+6)
c
c --- Calculate reaction rates
c
      r(  1) = rk(  1)*y(iNO2)
      r(  2) = rk(  2)*y(iO)
      r(  3) = rk(  3)*y(iO3)*y(iNO)
      r(  4) = rk(  4)*y(iO)*y(iNO2)
      r(  5) = rk(  5)*y(iO)*y(iNO2)
      r(  6) = rk(  6)*y(iO)*y(iNO)
      r(  7) = rk(  7)*y(iNO2)*y(iO3)
      r(  8) = rk(  8)*y(iO3)
      r(  9) = rk(  9)*y(iO3)
      r( 10) = rk( 10)*y(iO1D)
      r( 11) = rk( 11)*y(iO1D)*H2O
      r( 12) = rk( 12)*y(iO3)*y(iOH)
      r( 13) = rk( 13)*y(iO3)*y(iHO2)
      r( 14) = rk( 14)*y(iNO3)
      r( 15) = rk( 15)*y(iNO3)*y(iNO)
      r( 16) = rk( 16)*y(iNO3)*y(iNO2)
      r( 17) = rk( 17)*y(iNO3)*y(iNO2)
      r( 18) = rk( 18)*y(iN2O5)*H2O
      r( 19) = rk( 19)*y(iN2O5)
      r( 20) = rk( 20)*y(iNO)*y(iNO)
      r( 21) = rk( 21)*y(iNO)*y(iNO2)*H2O
      r( 22) = rk( 22)*y(iNO)*y(iOH)
      r( 23) = rk( 23)*y(iHONO)
      r( 24) = rk( 24)*y(iOH)*y(iHONO)
      r( 25) = rk( 25)*y(iHONO)*y(iHONO)
      r( 26) = rk( 26)*y(iNO2)*y(iOH)
      r( 27) = rk( 27)*y(iOH)*y(iHNO3)
      r( 28) = rk( 28)*y(iHO2)*y(iNO)
      r( 29) = rk( 29)*y(iHO2)*y(iNO2)
      r( 30) = rk( 30)*y(iPNA)
      r( 31) = rk( 31)*y(iOH)*y(iPNA)
      r( 32) = rk( 32)*y(iHO2)*y(iHO2)
      r( 33) = rk( 33)*y(iHO2)*y(iHO2)*H2O
      r( 34) = rk( 34)*y(iH2O2)
      r( 35) = rk( 35)*y(iOH)*y(iH2O2)
      r( 36) = rk( 36)*y(iOH)*y(iCO)
      r( 37) = rk( 37)*y(iFORM)*y(iOH)
      r( 38) = rk( 38)*y(iFORM)
      r( 39) = rk( 39)*y(iFORM)
      r( 40) = rk( 40)*y(iFORM)*y(iO)
      r( 41) = rk( 41)*y(iFORM)*y(iNO3)
      r( 42) = rk( 42)*y(iALD2)*y(iO)
      r( 43) = rk( 43)*y(iALD2)*y(iOH)
      r( 44) = rk( 44)*y(iALD2)*y(iNO3)
      r( 45) = rk( 45)*y(iALD2)
      r( 46) = rk( 46)*y(iC2O3)*y(iNO)
      r( 47) = rk( 47)*y(iC2O3)*y(iNO2)
      r( 48) = rk( 48)*y(iPAN)
      r( 49) = rk( 49)*y(iC2O3)*y(iC2O3)
      r( 50) = rk( 50)*y(iC2O3)*y(iHO2)
      r( 51) = rk( 51)*y(iOH)
      r( 52) = rk( 52)*y(iPAR)*y(iOH)
      r( 53) = rk( 53)*y(iROR)
      r( 54) = rk( 54)*y(iROR)
      r( 55) = rk( 55)*y(iROR)*y(iNO2)
      r( 56) = rk( 56)*y(iO)*y(iOLE)
      r( 57) = rk( 57)*y(iOH)*y(iOLE)
      r( 58) = rk( 58)*y(iO3)*y(iOLE)
      r( 59) = rk( 59)*y(iNO3)*y(iOLE)
      r( 60) = rk( 60)*y(iO)*y(iETH)
      r( 61) = rk( 61)*y(iOH)*y(iETH)
      r( 62) = rk( 62)*y(iO3)*y(iETH)
      r( 63) = rk( 63)*y(iTOL)*y(iOH)
      r( 64) = rk( 64)*y(iTO2)*y(iNO)
      r( 65) = rk( 65)*y(iTO2)
      r( 66) = rk( 66)*y(iOH)*y(iCRES)
      r( 67) = rk( 67)*y(iCRES)*y(iNO3)
      r( 68) = rk( 68)*y(iCRO)*y(iNO2)
      r( 69) = rk( 69)*y(iOPEN)
      r( 70) = rk( 70)*y(iOPEN)*y(iOH)
      r( 71) = rk( 71)*y(iOPEN)*y(iO3)
      r( 72) = rk( 72)*y(iOH)*y(iXYL)
      r( 73) = rk( 73)*y(iOH)*y(iMGLY)
      r( 74) = rk( 74)*y(iMGLY)
      r( 75) = rk( 75)*y(iO)*y(iISOP)
      r( 76) = rk( 76)*y(iOH)*y(iISOP)
      r( 77) = rk( 77)*y(iO3)*y(iISOP)
      r( 78) = rk( 78)*y(iNO3)*y(iISOP)
      r( 79) = rk( 79)*y(iXO2)*y(iNO)
      r( 80) = rk( 80)*y(iXO2)*y(iXO2)
      r( 81) = rk( 81)*y(iXO2N)*y(iNO)
      r( 82) = rk( 82)*y(iSO2)*y(iOH)
      r( 83) = rk( 83)*y(iSO2)
      r( 84) = rk( 84)*y(iMEOH)*y(iOH)
      r( 85) = rk( 85)*y(iETOH)*y(iOH)
      r( 86) = rk( 86)*y(iXO2)*y(iHO2)
      r( 87) = rk( 87)*y(iXO2N)*y(iHO2)
      r( 88) = rk( 88)*y(iXO2N)*y(iXO2N)
      r( 89) = rk( 89)*y(iXO2)*y(iXO2N)
      r( 90) = rk( 90)*y(iOH)*y(iHO2)
      r( 91) = rk( 91)*y(iCRO)
      r( 92) = rk( 92)*y(iOH)*y(iISPD)
      r( 93) = rk( 93)*y(iO3)*y(iISPD)
      r( 94) = rk( 94)*y(iNO3)*y(iISPD)
      r( 95) = rk( 95)*y(iISPD)
      r( 96) = rk( 96)*y(iNO2)*y(iISOP)
      r( 97) = rk( 97)*y(iO)*y(iOLE2)
      r( 98) = rk( 98)*y(iOH)*y(iOLE2)
      r( 99) = rk( 99)*y(iO3)*y(iOLE2)
      r(100) = rk(100)*y(iNO3)*y(iOLE2)
c
      return
      end
