SUBROUTINE computeFluxes(nn,mynn,dt,dx,fm,cl,dc,c6,vel,fp,cr,con)
    IMPLICIT NONE
    INCLUDE "camx.prm"
    
    INTEGER :: ii, nn, mynn
    REAL :: dt, dx, x, TWO3RDS
    REAL, DIMENSION(mynn) :: con,vel
    REAL, DIMENSION(MXCELLS) :: fm,cl,dc,c6,fp,cr
    
    !Compute fluxes from the parabolic distribution
    DO ii = 2,nn-1
        x = max(0., -vel(ii-1)*dt/dx)
        fm(ii) = x*(cl(ii) + 0.5*x*(dc(ii) + c6(ii)*(1. - TWO3RDS*x)))
        x = max(0., vel(ii)*dt/dx)
        fp(ii) = x*(cr(ii) - 0.5*x*(dc(ii) - c6(ii)*(1. - TWO3RDS*x)))
    END DO
    
    !Compute fluxes from boundary cells assuming uniform distribution
    IF (vel(1).gt.0.) THEN
        x = vel(1)*dt/dx
        fp(1) = x*con(1)
    END IF
    
    IF (vel(nn-1).lt.0.) THEN
        x = -vel(nn-1)*dt/dx
        fm(nn) = x*con(nn)
    END IF
    
    RETURN
    END