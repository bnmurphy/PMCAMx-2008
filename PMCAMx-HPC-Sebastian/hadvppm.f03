SUBROUTINE hadvppm(nn,dt,dx,con,vel,area,areav,flxarr,mynn)
    !this version is without STEEPEN
    IMPLICIT NONE
    INCLUDE "camx.prm"
    
    REAL :: dt, dx, x
    INTEGER :: mynn, nn, ii !is this actually necessary?
    REAL, DIMENSION(mynn) :: con,vel,area,areav,flxarr
    REAL, PARAMETER :: TWO3RDS=2./3.
    REAL, DIMENSION(MXCELLS) :: fm, fp, cm, cl, cr, dc, c6
    
    ! Initialize fluxes to zero. Either positive or negative flux will remain zero depending on the sign of the velocity zeta = dx*dx
    DO ii = 1,nn
        fm(ii) = 0.
        fp(ii) = 0.
    END DO
    
    !-----Zero order polynomial at the boundary cells
    cm(2)  = con(2)
    cm(nn) = con(nn-1)
    
    !First order polynomial at the next cells, no monotonicity constraint needed
    
    cm(3)    = (con(3) + con(2))/2.
    cm(nn-1) = (con(nn-1) + con(nn-2))/2.
    
    !Second order polynomial inside the domain
    CALL secondOrderPolynomial(nn,mynn,con,dc,cm,cr,cl)
    
    !note: for now I removed discontinuity capturing as included in CAMx, I can just put it into a function ;-)
    
    !Generate piecewise parabolic distributions    
    CALL generatePieceParDistri(nn,mynn,cr,con,cl,dc,c6)
    
    !Compute fluxes
    CALL computeFluxes(nn,mynn,dt,dx,fm,cl,dc,c6,vel,fp,cr,con)
    
    !Update concentrations
    CALL updateConcentrations(nn,mynn,flxarr,con,area,areav,fp,fm)
    
    RETURN
    END