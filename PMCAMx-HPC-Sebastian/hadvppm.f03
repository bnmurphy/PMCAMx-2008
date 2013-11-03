SUBROUTINE hadvppm(NN,dt,dx,con,vel,mscl,flxarr,flux1,flux2,saflux,fc1,fc2)

!include "camx.prm"
!todo where is camx.prm=
!@ to doto fix this define MX1D here for the moment
    
REAL, DIMENSION(NN) :: con,vel,flxarr,mscl,saflux
REAL :: flux1, flux2
INTEGER, PARAMETER :: MX1D = 100
REAL, PARAMETER :: TWO3RDS = 2./3.
REAL, DIMENSION(MX1D):: fm,fp,cm,cl,cr,dc,c6,d2c,eta,etabar,cld,crd,fc1,fc2

! Set all fluxes to zero. Either positive or negative flux will remain zero depending on the sign of the velocity zeta = dx*dx
!do i = 1,NN
!fm(i) = 0.
!fp(i) = 0.
!!
!!======================== Process Analysis Begin ====================================
!!
!fc1(i) = 0.
!fc2(i) = 0.
!
fm(:) = 0.
fp(:) = 0.
fc1(:) = 0.
fc2(:)= 0.

PRINT*, "fm(100) =", fm(5)

!========================= Process Analysis End =====================================
!
!enddo
!
!-----Zero order polynomial at the boundary cells
!
cm(2)  = con(2)
cm(NN) = con(NN-1)
!
!-----First order polynomial at the next cells, no monotonicity constraint
!     needed
!
cm(3)    = (con(3) + con(2))/2.
cm(NN-1) = (con(NN-1) + con(NN-2))/2.
!
!-----Second order polynomial inside the domain
!
DO i = 3,NN-2
    !
    !-----Compute average slope in the i'th cell
    !
    dc(i) = 0.5*(con(i+1) - con(i-1))
    !      
    !-----Guarantee that CM lies between CON(I) and CON(I+1)
    !     monotonicity constraint

    IF ((con(i+1) - con(i))*(con(i) - con(i-1)).gt.0.) THEN
        dc(i) = sign(1.,dc(i))*min(abs(dc(i)),2.*abs(con(i+1) - con(i)), 2.*abs(con(i) - con(i-1)))
    ELSE
        dc(i) = 0.
    END IF
END DO
!
DO i = 3,NN-3
    cm(i+1) = con(i) + 0.5*(con(i+1) - con(i)) + (dc(i) - dc(i+1))/6.
END DO
!
do i = 2,NN-1
cr(i) = cm(i+1)
cl(i) = cm(i)
enddo
!
!-----Optional discontinuty capturing
!     This is disbaled completely in this version 
!
!***      if (STEEPEN) then
!***        do i = 2,NN-1
!***          eta(i) = 0.
!***          cld(i) = con(i)
!***          crd(i) = con(i)
!***        enddo
!***! 
!***!-----Finite diff approximation to 2nd derivative
!***!
!***        do i = 3,NN-2
!***          d2c(i) = (con(i+1) - 2.*con(i) + con(i-1))/6.
!***        enddo
!***!
!***!-----No discontinuity detection near the boundary: cells 2, 3, NN-2, NN-1
!***! 
!***        do i = 4,NN-3  
!***! 
!***!-----Compute etabars
!***! 
!***          if ((-d2c(i+1)*d2c(i-1).gt.0.) .and.
!***     &        (abs(con(i+1) - con(i-1)) -
!***     &         EPS*min(abs(con(i+1)),abs(con(i-1))).gt.0.)) then
!***            etabar(i) = -zeta*(d2c(i+1) - d2c(i-1))/
!***     &                  (con(i+1) - con(i-1))
!***          else
!***            etabar(i) = 0.
!***          endif
!***          eta(i) = max(0.,min(ETA1*(etabar(i) - ETA2),1.)) 
!***          crd(i) = con(i+1) - 0.5*dc(i+1)
!***          cld(i) = con(i-1) + 0.5*dc(i-1)
!***        enddo
!***!
!***        do i = 2,NN-1
!***          cr(i) = cm(i+1) + eta(i)*(crd(i) - cm(i+1))
!***          cl(i) = cm(i) + eta(i)*(cld(i) - cm(i))
!***        enddo
!***      endif
!
!-----Generate piecewise parabolic distributions
!
do i = 2,NN-1
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
do i = 2,NN-1
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
if (vel(NN-1).lt.0.) then
x = -vel(NN-1)*dt/dx
fm(NN) = x*con(NN)
endif
!
!-----Update concentrations
!
flxarr(1) = (fp(1) - fm(2))*dx/dt
saflux(1) = flxarr(1)*dt/dx
do i = 2,NN-1
flxarr(i) = (fp(i) - fm(i+1))*dx/dt
con(i) = con(i) - mscl(i)*(flxarr(i) - flxarr(i-1))*dt/dx
saflux(i) = flxarr(i)*dt/dx
!
!======================== Process Analysis Begin ====================================
!
fc1(i) =   mscl(i)*flxarr(i-1)*dt/dx
fc2(i) = - mscl(i)*flxarr(i)*dt/dx
!
!========================= Process Analysis End =====================================
!
enddo
flux1 = mscl(2)*flxarr(1)
flux2 = mscl(NN-1)*flxarr(NN-1)
!
return
end
