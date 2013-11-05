PROGRAM unit_test_hadvppm
    !IMPLICIT NONE
    INTEGER :: ii = 0       
    INTEGER,PARAMETER :: nn = 100
    INTEGER,PARAMETER :: mynn = 100
    REAL :: dt,dx,sigma
    REAL, DIMENSION(nn) :: con,vel,area,areav,flxarr
    
    dt = 1.0
    dx = 360000
    sigma = 0.1
    
    DO ii=1,nn
        con(ii) = EXP((-(ii-9) * (ii-9)/2*sigma*sigma))
        vel(ii) = 150
        area(ii) = 1.
        areav(ii) = 1.
        flxarr(ii) = 0.
    END DO
    
    DO ii = 1,100000
        CALL hadvppm(nn,dt,dx,con,vel,area,areav,flxarr,mynn)
    END DO
    
    
    
END PROGRAM unit_test_hadvppm
