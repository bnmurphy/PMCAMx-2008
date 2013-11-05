SUBROUTINE hadvppm(nn,dt,dx,con,vel,area,areav,flxarr,mynn)
    !this version is without STEEPEN
    INCLUDE "camx.prm"
    !INTEGER :: mynn !is this actually necessary?
    REAL, DIMENSION(mynn) :: con,vel,area,areav,flxarr
    REAL, PARAMETER :: TWO3RDS=2./3.
    REAL, DIMENSION(MXCELLS) :: fm, fp, cm, cl, cr, dc, c6
    ! Set all fluxes to zero. Either positive or negative flux will remain zero depending on the sign of the velocity zeta = dx*dx
      
    DO i = 1,nn
        fm(i) = 0.
        fp(i) = 0.
    END DO
    
    !-----Zero order polynomial at the boundary cells
!
      cm(2)  = con(2)
      cm(nn) = con(nn-1)
!
!-----First order polynomial at the next cells, no monotonicity constraint
!     needed
!
      cm(3)    = (con(3) + con(2))/2.
      cm(nn-1) = (con(nn-1) + con(nn-2))/2.
!
!-----Second order polynomial inside the domain
!
      do i = 3,nn-2
!
!-----Compute average slope in the i'th cell
!
        dc(i) = 0.5*(con(i+1) - con(i-1))
!      
!-----Guarantee that CM lies between CON(I) and CON(I+1)
!     monotonicity constraint
      
        if ((con(i+1) - con(i))*(con(i) - con(i-1)).gt.0.) then
          dc(i) = sign(1.,dc(i))*min(&
                                    abs(dc(i)),&
                                    2.*abs(con(i+1) - con(i)),&
                                    2.*abs(con(i) - con(i-1)))
        else
          dc(i) = 0.
        endif
      enddo
!
      do i = 3,nn-3
        cm(i+1) = con(i) + &
                 0.5*(con(i+1) - con(i)) + (dc(i) - dc(i+1))/6.
      enddo
!
      do i = 2,nn-1
        cr(i) = cm(i+1)
        cl(i) = cm(i)
      enddo
!
!REMOVED optional discontinuity capturing
!
!-----Generate piecewise parabolic distributions
!
      do i = 2,nn-1
!
!-----Monotonicity
! 
        if ((cr(i) - con(i))*(con(i) - cl(i)).gt.0.) then
          dc(i) = cr(i) - cl(i)
          c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
!
!-----Overshoot cases
!
          if (dc(i)*c6(i) .gt. dc(i)*dc(i)) then
            cl(i) = 3.*con(i) - 2.*cr(i)
          elseif (-dc(i)*dc(i) .gt. dc(i)*c6(i)) then
            cr(i) = 3.*con(i) - 2.*cl(i)
          endif
        else
          cl(i) = con(i)
          cr(i) = con(i)
        endif
        dc(i) = cr(i) - cl(i)
        c6(i) = 6.*(con(i) - 0.5*(cl(i) + cr(i)))
      enddo
!
!-----Compute fluxes from the parabolic distribution
!
      do i = 2,nn-1
        x = max(0., -vel(i-1)*dt/dx)
        fm(i) = x*(cl(i) + 0.5*x*(dc(i) + c6(i)*(1. - TWO3RDS*x)))
        x = max(0., vel(i)*dt/dx)
        fp(i) = x*(cr(i) - 0.5*x*(dc(i) - c6(i)*(1. - TWO3RDS*x)))
      enddo
!
!-----Compute fluxes from boundary cells assuming uniform distribution
!
      if (vel(1).gt.0.) then
        x = vel(1)*dt/dx
        fp(1) = x*con(1)
      endif
!
      if (vel(nn-1).lt.0.) then
        x = -vel(nn-1)*dt/dx
        fm(nn) = x*con(nn)
      endif
!
!-----Update concentrations
!
      flxarr(1) = (fp(1) - fm(2))
      do i = 2,nn-1
        flxarr(i) = (fp(i) - fm(i+1))
        con(i) = con(i) - &
               area(i)*(areav(i)*flxarr(i) - areav(i-1)*flxarr(i-1))
      enddo
!
      return
      end
