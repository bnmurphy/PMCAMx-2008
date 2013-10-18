      subroutine ratejac6(nstrt,neq1,cncrad,conc,r,rate,loss,jac)
c
c-----CAMx v4.02 030709
c
c     RATEJAC computes reaction rate and its Jacobian matrix for fast
c     state species
c
c     Copyright 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Routines Called:
c        none
c
c     Called by:
c        TRAP
c
      include "camx.prm"
      include "chmstry.com"
c
      real loss(MXSPEC+1),gain(MXSPEC+1),rate(MXSPEC+1),
     &     jac(MXSPEC,MXSPEC),conc(MXSPEC+1),
     &     cncrad(MXRADCL),r(MXRXN)
c
      neqtmp = neq1
      neqtmp = min0(neqtmp, nspfst)
      do l=1,neqtmp
        Loss(l) = 0.
        Gain(l) = 0.
      enddo
c
      do i=1,neqtmp
        do j=1,neqtmp
          jac(i,j) = 0.
        enddo
      enddo
c
c  decide whether NO is solved here
c
c
c      NO2, and O3 are solved together
c
      if(neq1.ge.3) then
c
c   NO2    O3
c
        Loss(kNO2  )= +( 1.000)*r(  1)+( 1.000)*r(  4)+( 1.000)*r(  5)
     &                +( 1.000)*r(  7)+( 1.000)*r( 16)+( 1.000)*r( 17)
     &                +( 1.000)*r( 21)+( 1.000)*r( 26)+( 1.000)*r( 29)
     &                +( 1.000)*r( 47)+( 1.000)*r( 55)+( 1.000)*r( 68)
     &                +( 1.000)*r( 96)
        Gain(kNO2  )= +( 1.000)*r(  3)+( 1.000)*r(  6)+( 0.890)*r( 14)
     &                +( 2.000)*r( 15)+( 1.000)*r( 16)+( 1.000)*r( 19)
     &                +( 2.000)*r( 20)+( 1.000)*r( 24)+( 1.000)*r( 25)
     &                +( 1.000)*r( 28)+( 1.000)*r( 30)+( 1.000)*r( 31)
     &                +( 1.000)*r( 46)+( 1.000)*r( 48)+( 1.000)*r( 59)
     &                +( 0.900)*r( 64)+( 0.200)*r( 78)+( 1.000)*r( 79)
     &                +( 1.000)*r(100)
        Loss(kO3   )= +( 1.000)*r(  3)+( 1.000)*r(  7)+( 1.000)*r(  8)
     &                +( 1.000)*r(  9)+( 1.000)*r( 12)+( 1.000)*r( 13)
     &                +( 1.000)*r( 58)+( 1.000)*r( 62)+( 1.000)*r( 71)
     &                +( 1.000)*r( 77)+( 1.000)*r( 93)+( 1.000)*r( 99)
        Gain(kO3   )= +( 1.000)*r(  2)

          JAC(kNO2 ,kNO2 )= +( 1.000)*r(  1)+( 1.000)*r(  4)
     &                      +( 1.000)*r(  5)+( 1.000)*r(  7)
     &                      +( 1.000)*r( 16)+(-1.000)*r( 16)
     &                      +( 1.000)*r( 17)+( 1.000)*r( 21)
     &                      +( 1.000)*r( 26)+( 1.000)*r( 29)
     &                      +( 1.000)*r( 47)+( 1.000)*r( 55)
     &                      +( 1.000)*r( 68)+( 1.000)*r( 96)
          JAC(kNO2 ,kO3  )= +(-1.000)*r(  3)+( 1.000)*r(  7)
          JAC(kO3  ,kNO2 )= +( 1.000)*r(  7)
          JAC(kO3  ,kO3  )= +( 1.000)*r(  3)+( 1.000)*r(  7)
     &                      +( 1.000)*r(  8)+( 1.000)*r(  9)
     &                      +( 1.000)*r( 12)+( 1.000)*r( 13)
     &                      +( 1.000)*r( 58)+( 1.000)*r( 62)
     &                      +( 1.000)*r( 71)+( 1.000)*r( 77)
     &                      +( 1.000)*r( 93)+( 1.000)*r( 99)
      endif
c
c     PAN is added to the coupled species
c
      if(neq1.ge.4) then
c
c   PAN
c
        Loss(kPAN  )= +( 1.000)*r( 48)
        Gain(kPAN  )= +( 1.000)*r( 47)

          JAC(kNO2 ,kPAN )= +(-1.000)*r( 48)

          JAC(kPAN ,kNO2 )= +(-1.000)*r( 47)

          JAC(kPAN ,kPAN )= +( 1.000)*r( 48)
      endif
c
c  add NO during day time or night time when NO is not zero
c
      if(nstrt.eq.1) then
        if(neq1.ge.3) then
c
c    NO
c
        Loss(kNO   )= +( 1.000)*r(  3)+( 1.000)*r(  6)+( 1.000)*r( 15)
     &                +( 2.000)*r( 20)+( 1.000)*r( 21)+( 1.000)*r( 22)
     &                +( 1.000)*r( 28)+( 1.000)*r( 46)+( 1.000)*r( 64)
     &                +( 1.000)*r( 79)+( 1.000)*r( 81)
        Gain(kNO   )= +( 1.000)*r(  1)+( 1.000)*r(  4)+( 0.110)*r( 14)
     &                +( 1.000)*r( 16)+( 1.000)*r( 23)+( 1.000)*r( 25)
     &                +( 0.200)*r( 96)

          JAC(kNO  ,kNO  )= +( 1.000)*r(  3)+( 1.000)*r(  6)
     &                      +( 1.000)*r( 15)+( 4.000)*r( 20)
     &                      +( 1.000)*r( 21)+( 1.000)*r( 22)
     &                      +( 1.000)*r( 28)+( 1.000)*r( 46)
     &                      +( 1.000)*r( 64)+( 1.000)*r( 79)
     &                      +( 1.000)*r( 81)

          JAC(kNO  ,kNO2 )= +(-1.000)*r(  1)+(-1.000)*r(  4)
     &                      +(-1.000)*r( 16)+( 1.000)*r( 21)
     &                      +(-0.200)*r( 96)
          JAC(kNO  ,kO3  )= +( 1.000)*r(  3)

          JAC(kNO2 ,kNO  )= +(-1.000)*r(  3)+(-1.000)*r(  6)
     &                      +(-2.000)*r( 15)+(-4.000)*r( 20)
     &                      +( 1.000)*r( 21)+(-1.000)*r( 28)
     &                      +(-1.000)*r( 46)+(-0.900)*r( 64)
     &                      +(-1.000)*r( 79)
          JAC(kO3  ,kNO  )= +( 1.000)*r(  3)
        endif
        if(neq1.ge.4) then


        endif
      endif
c
c  complete Jacobian terms
c
      do j=1,neqtmp
        do i=1,neqtmp
          jac(i,j) = jac(i,j)/conc(j)
        enddo
      enddo
c
c
      do l=1,neq1
        rate(l) = Gain(l) - Loss(l)
      enddo
c
      return
      end
