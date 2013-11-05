SUBROUTINE secondOrderPolynomial(nn,mynn,con,dc,cm,cr,cl)
    !Second order polynomial inside the domain
    
    IMPLICIT NONE
    INCLUDE "camx.prm"
    
    INTEGER ii,nn,mynn
    REAL, DIMENSION(mynn) :: con
    REAL, DIMENSION(MXCELLS) :: dc,cm,cr,cl
    
    DO ii = 3,nn-2
        !Compute average slope in the ii'th cell
        dc(ii) = 0.5*(con(ii+1) - con(ii-1))
        !Guarantee that CM lies between CON(I) and CON(I+1) monotonicity constraint
        IF ((con(ii+1) - con(ii))*(con(ii) - con(ii-1)).gt.0.) THEN
            dc(ii) = sign(1.,dc(ii))*min(abs(dc(ii)), 2.*abs(con(ii+1) - con(ii)), 2.*abs(con(ii) - con(ii-1)))
        ELSE
            dc(ii) = 0.
        END IF
    END DO

    DO ii = 3,nn-3
        cm(ii+1) = con(ii) + 0.5*(con(ii+1) - con(ii)) + (dc(ii) - dc(ii+1))/6.
    END DO

    DO ii = 2,nn-1
        cr(ii) = cm(ii+1)
        cl(ii) = cm(ii)
    END DO
    
    RETURN
    END