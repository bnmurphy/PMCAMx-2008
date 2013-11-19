SUBROUTINE secondOrderPolynomial(nn,con,dc,cm,cr,cl)
    !Second order polynomial inside the domain
    
    IMPLICIT NONE
    INCLUDE "camx.prm"
    
    INTEGER jj,nn
    REAL, DIMENSION(nn) :: con
    REAL, DIMENSION(nn) :: dc,cm,cr,cl
    
    DO jj = 3,nn-2
        !Compute average slope in the jj'th cell
        dc(jj) = 0.5*(con(jj+1) - con(jj-1))
        !Guarantee that CM lies between CON(I) and CON(I+1) monotonicity constraint
        IF ((con(jj+1) - con(jj))*(con(jj) - con(jj-1)) .gt. 0.) THEN
            dc(jj) = sign(1.,dc(jj))*min(abs(dc(jj)), 2.*abs(con(jj+1) - con(jj)), 2.*abs(con(jj) - con(jj-1)))
        ELSE
            dc(jj) = 0.
        END IF
    END DO

    DO jj = 3,nn-3
        cm(jj+1) = con(jj) + 0.5*(con(jj+1) - con(jj)) + (dc(jj) - dc(jj+1))/6.
    END DO

    DO jj = 2,nn-1
        cr(jj) = cm(jj+1)
        cl(jj) = cm(jj)
    END DO
    
    RETURN
    END