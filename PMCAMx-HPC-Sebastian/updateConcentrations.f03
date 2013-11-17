SUBROUTINE updateConcentrations(nn,flxarr,con,area,areav,fp,fm)
    IMPLICIT NONE
    INCLUDE "camx.prm"
    
    INTEGER ::  ii,nn
    REAL, DIMENSION(nn) :: flxarr,con,area,areav
    REAL, DIMENSION(nn) :: fp,fm
    
    flxarr(1) = (fp(1) - fm(2))
    DO ii = 2,nn-1
        flxarr(ii) = (fp(ii) - fm(ii+1))
        con(ii) = con(ii) - area(ii)*(areav(ii)*flxarr(ii) - areav(ii-1)*flxarr(ii-1))
    END DO
    
    RETURN
    END
    
