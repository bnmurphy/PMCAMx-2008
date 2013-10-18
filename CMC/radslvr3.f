      subroutine radslvr3(ldark,H2O,M,O2,CH4,H2,cncrad,conc,r)
c
c-----CAMx v4.03 031205
c
c     RADSLVR3 initializes radical concentrations
c
c     Copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Routines Called:
c        CPIVOT
c
c     Called by:
c        TRAP
c
      include "camx.prm"
      include "chmstry.com"
      include "filunit.com"
c
      logical ldark, lusecp
      real M
      real conc(MXSPEC+1),cncrad(MXRADCL),r(MXRXN),rate(MXSPEC+1),
     &     loss(MXSPEC+1),gain(MXSPEC+1),jac(MXRADCL,MXRADCL),
     &     crold(MXRADCL)
      data tol/0.01/
c
c
c  new group of radicals
c
c   O1D
c
c
c  this radical is formed by photolysis only
c  solved by direct substitution
c
      if (ldark) then
        cncrad(kO1D) = 0.
      else
        Loss(kO1D  )= +( 1.000)*r( 10)+( 1.000)*r( 11)
        Gain(kO1D  )= +( 1.000)*r(  9)
        cncrad(kO1D) = gain(kO1D)/loss(kO1D)*cncrad(kO1D)
        cncrad(kO1D) = amax1(bdlrad, cncrad(kO1D))
      r( 10) = rk( 10)*cncrad(kO1D)
      r( 11) = rk( 11)*cncrad(kO1D)*H2O
      endif
c
c  new group of radicals
c
c     O
c
c
c  this radical is formed by photolysis only
c  solved by direct substitution
c
      if (ldark) then
        cncrad(kO) = 0.
      else
        Loss(kO    )= +( 1.000)*r(  2)+( 1.000)*r(  4)+( 1.000)*r(  5)
     &                +( 1.000)*r(  6)+( 1.000)*r( 40)+( 1.000)*r( 42)
     &                +( 1.000)*r( 56)+( 1.000)*r( 60)+( 1.000)*r( 75)
        Gain(kO    )= +( 1.000)*r(  1)+( 1.000)*r(  8)+( 1.000)*r( 10)
     &                +( 0.890)*r( 14)
        cncrad(kO) = gain(kO)/loss(kO)*cncrad(kO)
        cncrad(kO) = amax1(bdlrad, cncrad(kO))
      r(  2) = rk(  2)*cncrad(kO)
      r(  4) = rk(  4)*cncrad(kO)*conc(kNO2)
      r(  5) = rk(  5)*cncrad(kO)*conc(kNO2)
      r(  6) = rk(  6)*cncrad(kO)*conc(kNO)
      r( 40) = rk( 40)*conc(kFORM)*cncrad(kO)
      r( 42) = rk( 42)*conc(kALD2)*cncrad(kO)
      r( 56) = rk( 56)*cncrad(kO)*conc(kOLE)
      r( 60) = rk( 60)*cncrad(kO)*conc(kETH)
      r( 75) = rk( 75)*cncrad(kO)*conc(kISOP)
      endif
c
c  new group of radicals
c
c  N2O5   NO3
c
      nstrt = kN2O5
      nend  = kNO3
      weight = 1.0
      errbig = 1.0
      thresh = 1.0e-15
      kount = 0
  13  kount = kount + 1
      if (kount.gt.100) then
        write(iout,'(//,A)') 'ERROR in RADSLVR3:'
        write(iout,*) 'No Convergence in RADSLVR3, errbig = ', errbig
        write(iout,*) 'statement number = ', 13
        write(iout,*) 'igrd,i, j, k = ', igrdchm,ichm,jchm,kchm
        write(iout,*) 'LDARK is set ', ldark
        do l=1,ngas
          write(iout,'(i3,2x,a7,1pe10.3)') l,spname(l),conc(l)
        enddo
        write(iout,*) 'The radicals are: '
        do l= 1 , nrad
          write(iout,'(i3,2x,a7,1pe10.3)')
     &          l,nmrad(l),cncrad(l)
        enddo
        write(iout,*) 'Currently solving ', nstrt, ' to ',nendcp
        do l=nstrt,nendcp
          write(iout,'(i3,2x,a7,1p2e10.3)')
     &          l,nmrad(l),cncrad(l),abs(rate(l))/cncrad(l)
        enddo
        call camxerr()
      endif
      do i=nstrt,nend
        do j=nstrt,nend
          jac(i,j) = 0.
        enddo
      enddo
        Loss(kN2O5 )= +( 1.000)*r( 18)+( 1.000)*r( 19)
        Gain(kN2O5 )= +( 1.000)*r( 17)
        Loss(kNO3  )= +( 1.000)*r( 14)+( 1.000)*r( 15)+( 1.000)*r( 16)
     &                +( 1.000)*r( 17)+( 1.000)*r( 41)+( 1.000)*r( 44)
     &                +( 1.000)*r( 59)+( 1.000)*r( 67)+( 1.000)*r( 78)
     &                +( 1.000)*r( 94)
        Gain(kNO3  )= +( 1.000)*r(  5)+( 1.000)*r(  7)+( 1.000)*r( 19)
     &                +( 1.000)*r( 27)

          JAC(kN2O5,kN2O5)= +( 1.000)*r( 18)+( 1.000)*r( 19)
          JAC(kN2O5,kNO3 )= +(-1.000)*r( 17)
          JAC(kNO3 ,kN2O5)= +(-1.000)*r( 19)
          JAC(kNO3 ,kNO3 )= +( 1.000)*r( 14)+( 1.000)*r( 15)
     &                      +( 1.000)*r( 16)+( 1.000)*r( 17)
     &                      +( 1.000)*r( 41)+( 1.000)*r( 44)
     &                      +( 1.000)*r( 59)+( 1.000)*r( 67)
     &                      +( 1.000)*r( 78)+( 1.000)*r( 94)
c
c  complete the rates and Jacobian
c
      do i=nstrt,nend
        crold(i) = cncrad(i)
        rate(i) = gain(i) - loss(i)
        do j=nstrt,nend
          jac(i,j) = jac(i,j)/cncrad(j)
        enddo
      enddo
c
c
c  2 X 2 matrix is solved directly
c
      nendcp = nend
      det = jac(nstrt,nstrt)*jac(nend,nend)
     &     -jac(nend,nstrt)*jac(nstrt,nend)
      bnstrt = rate(nstrt)*jac(nend,nend) - rate(nend)*jac(nstrt,nend)
      bnend = jac(nstrt,nstrt)*rate(nend) - jac(nend,nstrt)*rate(nstrt)
      rate(nstrt) = bnstrt/det
      rate(nend) = bnend/det
c
c  update radical concentrations
c
      errbig = 0.
      if (kount.gt.5) weight = 0.5
      do l=nstrt,nend-0
        if (cncrad(l) .le. thresh) then
          err = abs((cncrad(l)/crold(l))-1.0)
          err = amin1(err,thresh)
        else
          err = abs(rate(l))/cncrad(l)
        endif
        errbig = amax1(errbig,err)
        if (rate(l)*weight.gt.10.0*cncrad(l)) then
          cncrad(l) = cncrad(l)*10.
        elseif (-rate(l)*weight.gt.0.9*cncrad(l)) then
          cncrad(l) = cncrad(l)/10.
        else
          cncrad(l) = cncrad(l) + rate(l)*weight
        endif
        cncrad(l) = amax1(cncrad(l), bdlrad)
      enddo
c
      r( 14) = rk( 14)*cncrad(kNO3)
      r( 15) = rk( 15)*cncrad(kNO3)*conc(kNO)
      r( 16) = rk( 16)*cncrad(kNO3)*conc(kNO2)
      r( 17) = rk( 17)*cncrad(kNO3)*conc(kNO2)
      r( 18) = rk( 18)*cncrad(kN2O5)*H2O
      r( 19) = rk( 19)*cncrad(kN2O5)
      r( 41) = rk( 41)*conc(kFORM)*cncrad(kNO3)
      r( 44) = rk( 44)*conc(kALD2)*cncrad(kNO3)
      r( 59) = rk( 59)*cncrad(kNO3)*conc(kOLE)
      r( 67) = rk( 67)*conc(kCRES)*cncrad(kNO3)
      r( 78) = rk( 78)*cncrad(kNO3)*conc(kISOP)
      r( 94) = rk( 94)*cncrad(kNO3)*conc(kISPD)
      if (errbig.gt.tol) goto 13
c
c  new group of radicals
c
c    OH   HO2  C2O3   XO2  XO2N   TO2   ROR
c
      nstrt = kOH
      nend  = kROR
      weight = 1.0
      errbig = 1.0
      thresh = 1.0e-15
      kount = 0
  14  kount = kount + 1
      if (kount.gt.1000) then
        write(iout,'(//,A)') 'ERROR in RADSLVR3:'
        write(iout,*) 'No Convergence in RADSLVR3, errbig = ', errbig
        write(iout,*) 'statement number = ', 14
        write(iout,*) 'igrd,i, j, k = ', igrdchm,ichm,jchm,kchm
        write(iout,*) 'LDARK is set ', ldark
        do l=1,ngas
          write(iout,'(i3,2x,a7,1pe10.3)') l,spname(l),conc(l)
        enddo
        write(iout,*) 'The radicals are: '
        do l= 1 , nrad
          write(iout,'(i3,2x,a7,1pe10.3)')
     &          l,nmrad(l),cncrad(l)
        enddo
        write(iout,*) 'Currently solving ', nstrt, ' to ',nendcp
        do l=nstrt,nendcp
          write(iout,'(i3,2x,a7,1p2e10.3)')
     &          l,nmrad(l),cncrad(l),abs(rate(l))/cncrad(l)
        enddo
        call camxerr()
      endif
      do i=nstrt,nend
        do j=nstrt,nend
          jac(i,j) = 0.
        enddo
      enddo
c
c  to reduce the matrix size, TO2 is solved first and substituted
c
        Loss(kTO2  )= +( 1.000)*r( 64)+( 1.000)*r( 65)
        Gain(kTO2  )= +( 0.560)*r( 63)+( 0.300)*r( 72)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kTO2) = gain(kTO2)/loss(kTO2)*cncrad(kTO2)
      cncrad(kTO2) = amax1(bdlrad, cncrad(kTO2))
      r( 64) = rk( 64)*cncrad(kTO2)*conc(kNO)
      r( 65) = rk( 65)*cncrad(kTO2)
c
c  to reduce the matrix size, ROR is solved first and substituted
c
        Loss(kROR  )= +( 1.000)*r( 53)+( 1.000)*r( 54)+( 1.000)*r( 55)
        Gain(kROR  )= +( 0.760)*r( 52)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kROR) = gain(kROR)/loss(kROR)*cncrad(kROR)
      cncrad(kROR) = amax1(bdlrad, cncrad(kROR))
      r( 53) = rk( 53)*cncrad(kROR)
      r( 54) = rk( 54)*cncrad(kROR)
      r( 55) = rk( 55)*cncrad(kROR)*conc(kNO2)
        Loss(kOH   )= +( 1.000)*r( 12)+( 1.000)*r( 22)+( 1.000)*r( 24)
     &                +( 1.000)*r( 26)+( 1.000)*r( 27)+( 1.000)*r( 31)
     &                +( 1.000)*r( 35)+( 1.000)*r( 36)+( 1.000)*r( 37)
     &                +( 1.000)*r( 43)+( 1.000)*r( 51)+( 1.000)*r( 52)
     &                +( 1.000)*r( 57)+( 1.000)*r( 61)+( 1.000)*r( 63)
     &                +( 1.000)*r( 66)+( 1.000)*r( 70)+( 1.000)*r( 72)
     &                +( 1.000)*r( 73)+( 1.000)*r( 76)+( 1.000)*r( 82)
     &                +( 1.000)*r( 84)+( 1.000)*r( 85)+( 1.000)*r( 90)
     &                +( 1.000)*r( 92)
        Gain(kOH   )= +( 2.000)*r( 11)+( 1.000)*r( 13)+( 1.000)*r( 23)
     &                +( 1.000)*r( 28)+( 2.000)*r( 34)+( 1.000)*r( 40)
     &                +( 1.000)*r( 42)+( 0.790)*r( 50)+( 0.200)*r( 56)
     &                +( 0.100)*r( 58)+( 0.300)*r( 60)+( 0.080)*r( 71)
     &                +( 0.266)*r( 77)+( 0.268)*r( 93)
        Loss(kHO2  )= +( 1.000)*r( 13)+( 1.000)*r( 28)+( 1.000)*r( 29)
     &                +( 2.000)*r( 32)+( 2.000)*r( 33)+( 1.000)*r( 50)
     &                +( 1.000)*r( 86)+( 1.000)*r( 87)+( 1.000)*r( 90)
        Gain(kHO2  )= +( 1.000)*r( 12)+( 1.000)*r( 30)+( 1.000)*r( 35)
     &                +( 1.000)*r( 36)+( 1.000)*r( 37)+( 2.000)*r( 38)
     &                +( 1.000)*r( 40)+( 1.000)*r( 41)+( 2.000)*r( 45)
     &                +( 1.000)*r( 46)+( 2.000)*r( 49)+( 0.790)*r( 50)
     &                +( 1.000)*r( 51)+( 0.110)*r( 52)+( 0.940)*r( 53)
     &                +( 1.000)*r( 54)+( 0.380)*r( 56)+( 1.000)*r( 57)
     &                +( 0.440)*r( 58)+( 1.700)*r( 60)+( 1.000)*r( 61)
     &                +( 0.120)*r( 62)+( 0.440)*r( 63)+( 0.900)*r( 64)
     &                +( 1.000)*r( 65)+( 0.600)*r( 66)+( 1.000)*r( 69)
        Gain(kHO2  ) = Gain(kHO2  )
     &                +( 2.000)*r( 70)+( 0.760)*r( 71)+( 0.700)*r( 72)
     &                +( 1.000)*r( 74)+( 0.250)*r( 75)+( 0.912)*r( 76)
     &                +( 0.066)*r( 77)+( 0.800)*r( 78)+( 1.000)*r( 82)
     &                +( 1.000)*r( 84)+( 1.000)*r( 85)+( 0.503)*r( 92)
     &                +( 0.154)*r( 93)+( 0.925)*r( 94)+( 1.033)*r( 95)
     &                +( 0.800)*r( 96)
        Loss(kC2O3 )= +( 1.000)*r( 46)+( 1.000)*r( 47)+( 2.000)*r( 49)
     &                +( 1.000)*r( 50)
        Gain(kC2O3 )= +( 1.000)*r( 42)+( 1.000)*r( 43)+( 1.000)*r( 44)
     &                +( 1.000)*r( 48)+( 1.000)*r( 69)+( 1.000)*r( 70)
     &                +( 0.620)*r( 71)+( 1.000)*r( 73)+( 1.000)*r( 74)
     &                +( 0.250)*r( 75)+( 0.200)*r( 77)+( 0.498)*r( 92)
     &                +( 0.114)*r( 93)+( 0.075)*r( 94)+( 0.967)*r( 95)
        Loss(kXO2  )= +( 1.000)*r( 79)+( 2.000)*r( 80)+( 1.000)*r( 86)
     &                +( 1.000)*r( 89)
        Gain(kXO2  )= +( 1.000)*r( 45)+( 1.000)*r( 46)+( 2.000)*r( 49)
     &                +( 0.790)*r( 50)+( 1.000)*r( 51)+( 0.870)*r( 52)
     &                +( 0.960)*r( 53)+( 0.280)*r( 56)+( 1.000)*r( 57)
     &                +( 0.220)*r( 58)+( 0.910)*r( 59)+( 0.700)*r( 60)
     &                +( 1.000)*r( 61)+( 0.080)*r( 63)+( 0.600)*r( 66)
     &                +( 1.000)*r( 70)+( 0.030)*r( 71)+( 0.500)*r( 72)
     &                +( 1.000)*r( 73)+( 0.250)*r( 75)+( 0.991)*r( 76)
     &                +( 0.200)*r( 77)+( 1.000)*r( 78)+( 0.713)*r( 92)
     &                +( 0.064)*r( 93)+( 0.075)*r( 94)+( 0.700)*r( 95)
        Gain(kXO2  ) = Gain(kXO2  )
     &                +( 1.000)*r( 96)
        Loss(kXO2N )= +( 1.000)*r( 81)+( 1.000)*r( 87)+( 2.000)*r( 88)
     &                +( 1.000)*r( 89)
        Gain(kXO2N )= +( 0.130)*r( 52)+( 0.040)*r( 53)+( 0.020)*r( 56)
     &                +( 0.090)*r( 59)+( 0.088)*r( 76)
        Loss(kTO2  )= +( 1.000)*r( 64)+( 1.000)*r( 65)
        Gain(kTO2  )= +( 0.560)*r( 63)+( 0.300)*r( 72)
        Loss(kROR  )= +( 1.000)*r( 53)+( 1.000)*r( 54)+( 1.000)*r( 55)
        Gain(kROR  )= +( 0.760)*r( 52)

          JAC(kOH  ,kOH  )= +( 1.000)*r( 12)+( 1.000)*r( 22)
     &                      +( 1.000)*r( 24)+( 1.000)*r( 26)
     &                      +( 1.000)*r( 27)+( 1.000)*r( 31)
     &                      +( 1.000)*r( 35)+( 1.000)*r( 36)
     &                      +( 1.000)*r( 37)+( 1.000)*r( 43)
     &                      +( 1.000)*r( 51)+( 1.000)*r( 52)
     &                      +( 1.000)*r( 57)+( 1.000)*r( 61)
     &                      +( 1.000)*r( 63)+( 1.000)*r( 66)
     &                      +( 1.000)*r( 70)+( 1.000)*r( 72)
          JAC(kOH  ,kOH  )=JAC(kOH  ,kOH  )
     &                      +( 1.000)*r( 73)+( 1.000)*r( 76)
     &                      +( 1.000)*r( 82)+( 1.000)*r( 84)
     &                      +( 1.000)*r( 85)+( 1.000)*r( 90)
     &                      +( 1.000)*r( 92)
          JAC(kOH  ,kHO2 )= +(-1.000)*r( 13)+(-1.000)*r( 28)
     &                      +(-0.790)*r( 50)+( 1.000)*r( 90)
          JAC(kOH  ,kC2O3)= +(-0.790)*r( 50)
          JAC(kHO2 ,kOH  )= +(-1.000)*r( 12)+(-1.000)*r( 35)
     &                      +(-1.000)*r( 36)+(-1.000)*r( 37)
     &                      +(-1.000)*r( 51)+(-0.110)*r( 52)
     &                      +(-1.000)*r( 57)+(-1.000)*r( 61)
     &                      +(-0.440)*r( 63)+(-0.600)*r( 66)
     &                      +(-2.000)*r( 70)+(-0.700)*r( 72)
     &                      +(-0.912)*r( 76)+(-1.000)*r( 82)
     &                      +(-1.000)*r( 84)+(-1.000)*r( 85)
     &                      +( 1.000)*r( 90)+(-0.503)*r( 92)
          JAC(kHO2 ,kHO2 )= +( 1.000)*r( 13)+( 1.000)*r( 28)
     &                      +( 1.000)*r( 29)+( 4.000)*r( 32)
     &                      +( 4.000)*r( 33)+( 1.000)*r( 50)
     &                      +(-0.790)*r( 50)+( 1.000)*r( 86)
     &                      +( 1.000)*r( 87)+( 1.000)*r( 90)
          JAC(kHO2 ,kC2O3)= +(-1.000)*r( 46)+(-4.000)*r( 49)
     &                      +( 1.000)*r( 50)+(-0.790)*r( 50)
          JAC(kHO2 ,kXO2 )= +( 1.000)*r( 86)
          JAC(kHO2 ,kXO2N)= +( 1.000)*r( 87)
          JAC(kHO2 ,kTO2 )= +(-0.900)*r( 64)+(-1.000)*r( 65)
          JAC(kHO2 ,kROR )= +(-0.940)*r( 53)+(-1.000)*r( 54)
          JAC(kC2O3,kOH  )= +(-1.000)*r( 43)+(-1.000)*r( 70)
     &                      +(-1.000)*r( 73)+(-0.498)*r( 92)
          JAC(kC2O3,kHO2 )= +( 1.000)*r( 50)
          JAC(kC2O3,kC2O3)= +( 1.000)*r( 46)+( 1.000)*r( 47)
     &                      +( 4.000)*r( 49)+( 1.000)*r( 50)
          JAC(kXO2 ,kOH  )= +(-1.000)*r( 51)+(-0.870)*r( 52)
     &                      +(-1.000)*r( 57)+(-1.000)*r( 61)
     &                      +(-0.080)*r( 63)+(-0.600)*r( 66)
     &                      +(-1.000)*r( 70)+(-0.500)*r( 72)
     &                      +(-1.000)*r( 73)+(-0.991)*r( 76)
     &                      +(-0.713)*r( 92)
          JAC(kXO2 ,kHO2 )= +(-0.790)*r( 50)+( 1.000)*r( 86)
          JAC(kXO2 ,kC2O3)= +(-1.000)*r( 46)+(-4.000)*r( 49)
     &                      +(-0.790)*r( 50)
          JAC(kXO2 ,kXO2 )= +( 1.000)*r( 79)+( 4.000)*r( 80)
     &                      +( 1.000)*r( 86)+( 1.000)*r( 89)
          JAC(kXO2 ,kXO2N)= +( 1.000)*r( 89)
          JAC(kXO2 ,kROR )= +(-0.960)*r( 53)
          JAC(kXO2N,kOH  )= +(-0.130)*r( 52)+(-0.088)*r( 76)
          JAC(kXO2N,kHO2 )= +( 1.000)*r( 87)
          JAC(kXO2N,kXO2 )= +( 1.000)*r( 89)
          JAC(kXO2N,kXO2N)= +( 1.000)*r( 81)+( 1.000)*r( 87)
     &                      +( 4.000)*r( 88)+( 1.000)*r( 89)
          JAC(kXO2N,kROR )= +(-0.040)*r( 53)
          JAC(kTO2 ,kOH  )= +(-0.560)*r( 63)+(-0.300)*r( 72)
          JAC(kTO2 ,kTO2 )= +( 1.000)*r( 64)+( 1.000)*r( 65)
          JAC(kROR ,kOH  )= +(-0.760)*r( 52)
          JAC(kROR ,kROR )= +( 1.000)*r( 53)+( 1.000)*r( 54)
     &                      +( 1.000)*r( 55)
c
c  complete the rates and Jacobian
c
      do i=nstrt,nend
        crold(i) = cncrad(i)
        rate(i) = gain(i) - loss(i)
        do j=nstrt,nend
          jac(i,j) = jac(i,j)/cncrad(j)
        enddo
      enddo
c
c
c  Jacobian is modified due to substitution of TO2
c
c
c    OH   HO2  C2O3   XO2  XO2N   TO2   ROR
c
          JAC(kHO2 ,kOH  ) = JAC(kHO2 ,kOH  )
     &         - JAC(kHO2 ,kTO2 )*JAC(kTO2 ,kOH  )/JAC(kTO2 ,kTO2 )
c
c  end of substitution of TO2
c
c
c  Jacobian is modified due to substitution of ROR
c
c
c    OH   HO2  C2O3   XO2  XO2N   TO2   ROR
c
          JAC(kHO2 ,kOH  ) = JAC(kHO2 ,kOH  )
     &         - JAC(kHO2 ,kROR )*JAC(kROR ,kOH  )/JAC(kROR ,kROR )
          JAC(kXO2 ,kOH  ) = JAC(kXO2 ,kOH  )
     &         - JAC(kXO2 ,kROR )*JAC(kROR ,kOH  )/JAC(kROR ,kROR )
          JAC(kXO2N,kOH  ) = JAC(kXO2N,kOH  )
     &         - JAC(kXO2N,kROR )*JAC(kROR ,kOH  )/JAC(kROR ,kROR )
c
c  end of substitution of ROR
c
c
c  solve the matrix
c
      if (errbig.lt.500.0) then
        thresh = 1.0e-15
      else
        thresh = 1.0e+15
      endif
      lusecp = .false.
c
      nendcp=nend-2
      do j=nstrt,nendcp
        if (cncrad(j).le.thresh) then
          cncrad(j) = gain(j)/loss(j)*cncrad(j)
          do i=nstrt,nendcp
            jac(i,j) = 0.
            jac(j,i) = 0.
          enddo
          rate(j) = 0.
          jac(j,j) = 1.
      else
        lusecp = .true.
        endif
      enddo
c
      if (lusecp) then
        call cpivot(nstrt,nendcp,MXRADCL,jac,rate,ierr)
        if (ierr.ne.0) goto 900
      endif
c
c  update radical concentrations
c
      errbig = 0.
      weight = 1.
      do l=nstrt,nend-2
        if (cncrad(l) .le. thresh) then ! solve independantly
          rlim = 5.0
          weight = 0.5
          cncrad(l) = crold(l)*(1.-weight) +
     &                         weight*gain(l)/loss(l)*crold(l)
          cncrad(l) = amax1(cncrad(l),crold(l)/rlim)
          cncrad(l) = amin1(cncrad(l),crold(l)*rlim)
          err = abs((cncrad(l)/crold(l))-1.0)
          err = amin1(err,thresh)
        else                            ! part of coupled solution
          rlim = 10.0
          if (kount.gt.5) weight = 0.5
          rate(l) = amin1(rate(l), rlim*cncrad(l))
          rate(l) = amax1(rate(l), -cncrad(l)*(1.0-(1.0/rlim)))
          cncrad(l) = cncrad(l) + rate(l)*weight
          err = abs(rate(l))/cncrad(l)
        endif
        cncrad(l) = amax1(cncrad(l), bdlrad)
        cncrad(l) = amin1(cncrad(l), 10.0)
        errbig = amax1(errbig,err)
      enddo
c
      r( 12) = rk( 12)*conc(kO3)*cncrad(kOH)
      r( 13) = rk( 13)*conc(kO3)*cncrad(kHO2)
      r( 22) = rk( 22)*conc(kNO)*cncrad(kOH)
      r( 24) = rk( 24)*cncrad(kOH)*conc(kHONO)
      r( 26) = rk( 26)*conc(kNO2)*cncrad(kOH)
      r( 27) = rk( 27)*cncrad(kOH)*conc(kHNO3)
      r( 28) = rk( 28)*cncrad(kHO2)*conc(kNO)
      r( 29) = rk( 29)*cncrad(kHO2)*conc(kNO2)
      r( 31) = rk( 31)*cncrad(kOH)*conc(kPNA)
      r( 32) = rk( 32)*cncrad(kHO2)*cncrad(kHO2)
      r( 33) = rk( 33)*cncrad(kHO2)*cncrad(kHO2)*H2O
      r( 35) = rk( 35)*cncrad(kOH)*conc(kH2O2)
      r( 36) = rk( 36)*cncrad(kOH)*conc(kCO)
      r( 37) = rk( 37)*conc(kFORM)*cncrad(kOH)
      r( 43) = rk( 43)*conc(kALD2)*cncrad(kOH)
      r( 46) = rk( 46)*cncrad(kC2O3)*conc(kNO)
      r( 47) = rk( 47)*cncrad(kC2O3)*conc(kNO2)
      r( 49) = rk( 49)*cncrad(kC2O3)*cncrad(kC2O3)
      r( 50) = rk( 50)*cncrad(kC2O3)*cncrad(kHO2)
      r( 51) = rk( 51)*cncrad(kOH)
      r( 52) = rk( 52)*conc(kPAR)*cncrad(kOH)
      r( 53) = rk( 53)*cncrad(kROR)
      r( 54) = rk( 54)*cncrad(kROR)
      r( 55) = rk( 55)*cncrad(kROR)*conc(kNO2)
      r( 57) = rk( 57)*cncrad(kOH)*conc(kOLE)
      r( 61) = rk( 61)*cncrad(kOH)*conc(kETH)
      r( 63) = rk( 63)*conc(kTOL)*cncrad(kOH)
      r( 64) = rk( 64)*cncrad(kTO2)*conc(kNO)
      r( 65) = rk( 65)*cncrad(kTO2)
      r( 66) = rk( 66)*cncrad(kOH)*conc(kCRES)
      r( 70) = rk( 70)*conc(kOPEN)*cncrad(kOH)
      r( 72) = rk( 72)*cncrad(kOH)*conc(kXYL)
      r( 73) = rk( 73)*cncrad(kOH)*conc(kMGLY)
      r( 76) = rk( 76)*cncrad(kOH)*conc(kISOP)
      r( 79) = rk( 79)*cncrad(kXO2)*conc(kNO)
      r( 80) = rk( 80)*cncrad(kXO2)*cncrad(kXO2)
      r( 81) = rk( 81)*cncrad(kXO2N)*conc(kNO)
      r( 82) = rk( 82)*conc(kSO2)*cncrad(kOH)
      r( 84) = rk( 84)*conc(kMEOH)*cncrad(kOH)
      r( 85) = rk( 85)*conc(kETOH)*cncrad(kOH)
      r( 86) = rk( 86)*cncrad(kXO2)*cncrad(kHO2)
      r( 87) = rk( 87)*cncrad(kXO2N)*cncrad(kHO2)
      r( 88) = rk( 88)*cncrad(kXO2N)*cncrad(kXO2N)
      r( 89) = rk( 89)*cncrad(kXO2)*cncrad(kXO2N)
      r( 90) = rk( 90)*cncrad(kOH)*cncrad(kHO2)
      r( 92) = rk( 92)*cncrad(kOH)*conc(kISPD)
      if (errbig.gt.tol) goto 14
c
c  new group of radicals
c
c   TO2
c
        Loss(kTO2  )= +( 1.000)*r( 64)+( 1.000)*r( 65)
        Gain(kTO2  )= +( 0.560)*r( 63)+( 0.300)*r( 72)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kTO2) = gain(kTO2)/loss(kTO2)*cncrad(kTO2)
      cncrad(kTO2) = amax1(bdlrad, cncrad(kTO2))
      r( 64) = rk( 64)*cncrad(kTO2)*conc(kNO)
      r( 65) = rk( 65)*cncrad(kTO2)
c
c  new group of radicals
c
c   ROR
c
        Loss(kROR  )= +( 1.000)*r( 53)+( 1.000)*r( 54)+( 1.000)*r( 55)
        Gain(kROR  )= +( 0.760)*r( 52)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kROR) = gain(kROR)/loss(kROR)*cncrad(kROR)
      cncrad(kROR) = amax1(bdlrad, cncrad(kROR))
      r( 53) = rk( 53)*cncrad(kROR)
      r( 54) = rk( 54)*cncrad(kROR)
      r( 55) = rk( 55)*cncrad(kROR)*conc(kNO2)
c
c  new group of radicals
c
c   CRO
c
        Loss(kCRO  )= +( 1.000)*r( 68)+( 1.000)*r( 91)
        Gain(kCRO  )= +( 0.400)*r( 66)+( 1.000)*r( 67)
c
c  first order method chosen for this radical
c  solved by direct substitution
c
      cncrad(kCRO) = gain(kCRO)/loss(kCRO)*cncrad(kCRO)
      cncrad(kCRO) = amax1(bdlrad, cncrad(kCRO))
      r( 68) = rk( 68)*cncrad(kCRO)*conc(kNO2)
      r( 91) = rk( 91)*cncrad(kCRO)
c
      return
c
 900  continue
      write(iout,'(//,A)') 'ERROR in RADSLVR3:'
      write(iout,*) 'Zero determinant in CPIVOT at ', ierr
      write(iout,*) 'igrd,i, j, k = ', igrdchm,ichm,jchm,kchm
      write(iout,*) 'LDARK is set ', ldark
      do l=1,ngas
        write(iout,'(i3,2x,a7,1pe10.3)') l,spname(l),conc(l)
      enddo
      write(iout,*) 'The radicals are: '
      do l= 1 , nrad
        write(iout,'(i3,2x,a7,1pe10.3)')
     &        l,nmrad(l),cncrad(l)
      enddo
      write(iout,*) 'Currently solving ', nstrt, ' to ',nendcp
      do l=nstrt,nendcp
        write(iout,'(i3,2x,a7,1p2e10.3)')
     &        l,nmrad(l),cncrad(l),abs(rate(l))/cncrad(l)
      enddo
      call camxerr()
      end
