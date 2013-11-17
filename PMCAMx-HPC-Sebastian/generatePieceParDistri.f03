SUBROUTINE generatePieceParDistri(nn,cr,con,cl,dc,c6)
    IMPLICIT NONE
    INCLUDE "camx.prm"
    
    INTEGER ::  ii,nn
    REAL, DIMENSION(nn) :: con
    REAL, DIMENSION(nn) :: cr,cl,dc,c6
    
    DO ii = 2,nn-1
        !Monotonicity
        IF ((cr(ii) - con(ii))*(con(ii) - cl(ii)).gt.0.) THEN
            dc(ii) = cr(ii) - cl(ii)
            c6(ii) = 6.*(con(ii) - 0.5*(cl(ii) + cr(ii)))
            !Overshoot cases
            IF (dc(ii)*c6(ii) .gt. dc(ii)*dc(ii)) THEN
                cl(ii) = 3.*con(ii) - 2.*cr(ii)
            ELSE IF (-dc(ii)*dc(ii) .gt. dc(ii)*c6(ii)) THEN
                cr(ii) = 3.*con(ii) - 2.*cl(ii)
            END IF
        ELSE
            cl(ii) = con(ii)
            cr(ii) = con(ii)
        END IF
        dc(ii) = cr(ii) - cl(ii)
        c6(ii) = 6.*(con(ii) - 0.5*(cl(ii) + cr(ii)))
    END DO
    
    RETURN
    END
