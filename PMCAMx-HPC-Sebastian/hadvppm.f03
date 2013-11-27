SUBROUTINE hadvppm(nn,dt,dx,con,vel,area,areav,flxarr)
    !this version is without STEEPEN
    !IMPLICIT NONE
    INCLUDE "camx.prm"
    
    !REAL :: dt, dx, x
    INTEGER :: ii, status !is this actually necessary?
    !REAL, DIMENSION(nn) :: con,vel,area,areav,flxarr
    !REAL, PARAMETER :: TWO3RDS=2./3.
    !REAL, DIMENSION(nn) :: fm, fp, cm, cl, cr, dc, c6
    
    real, allocatable, dimension(:) :: fm, fp, cm, cl, cr, dc, c6
    
    allocate(fm(nn),stat=status)
    allocate(fp(nn),stat=status)
    allocate(cm(nn),stat=status)
    allocate(cl(nn),stat=status)
    allocate(cr(nn),stat=status)
    allocate(dc(nn),stat=status)
    allocate(c6(nn),stat = status)
    
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
    CALL secondOrderPolynomial(nn,con,dc,cm,cr,cl)
    
    !note: for now I removed discontinuity capturing as included in CAMx, I can just put it into a function ;-)
    
    !Generate piecewise parabolic distributions    
    !CALL generatePieceParDistri(nn,cr,con,cl,dc,c6)
    
    !Compute fluxes
    !CALL computeFluxes(nn,dt,dx,fm,cl,dc,c6,vel,fp,cr,con)
    
    !Update concentrations
    !CALL updateConcentrations(nn,flxarr,con,area,areav,fp,fm)
    
    deallocate(fm,stat=status)
    deallocate(fp,stat=status)
    deallocate(cm,stat=status)
    deallocate(cl,stat=status)
    deallocate(cr,stat=status)
    deallocate(dc,stat=status)
    deallocate(c6,stat = status)
    
    RETURN
    END