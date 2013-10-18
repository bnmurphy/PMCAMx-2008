      subroutine ktherm(tcell,pcell)
c
c-----CAMx v4.02 030709
c
c     KTHERM sets rate constants for thermo-chemical reactions
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        tcell               temperature (K)
c        pcell               pressure (mb)
c
c     Output arguments:
c        none
c
c     Routines Called:
c        none
c
c     Called by:
c        CHEMDRIV
c        PIGDRIVE
c        READCHM
c
      include "camx.prm"
      include "chmstry.com"
      include "filunit.com"
c
c-----Entry point
c
c-----Range check temperature and pressure
c
      if (tcell.le.TEMPLO .or. tcell.ge.TEMPHI) then
        write(iout,'(//,a)') 'ERROR in KTHERM:'
        write(iout,*) 'Temperature out of lookup table range'
        write(iout,*) 'Modify camx.prm to extend range to include'
        write(iout,'(A,F8.2)') 'Temp(K) = ', tcell
        call camxerr()
      endif
      if (pcell.le.PRESLO .or. pcell.ge.PRESHI) then
        write(iout,'(//,a)') 'ERROR in KTHERM:'
        write(iout,*) 'Pressure out of lookup table range'
        write(iout,*) 'Modify camx.prm to extend range to include'
        write(iout,'(A,F8.2)') 'Press(mb) = ', pcell
        call camxerr()
      endif
c
c-----Determine weighting factors for lookup table
c
      dtemp = (TEMPHI-TEMPLO) / (NTEMPR-1)
      j = int((tcell - tempr(1))/dtemp) + 1
      tmp = (tcell - tempr(j))/dtemp
      w1 = 1. - tmp
      dpres = (PRESHI-PRESLO) / (NPRESR-1)
      k = int((pcell - presr(1))/dpres) + 1
      tmp = (pcell - presr(k))/dpres
      w2 = 1. - tmp
c
c-----Look up rate constants
c
      do i = 1,nreact
        rk11 = rktbl(i,j,k)
        rk12 = rktbl(i,j,k+1)
        rk21 = rktbl(i,j+1,k)
        rk22 = rktbl(i,j+1,k+1)
        tmp1 = w1*rk11 + (1.-w1)*rk21
        tmp2 = w1*rk12 + (1.-w1)*rk22
        rk(i) = w2*tmp1 + (1.-w2)*tmp2
      enddo
c
      end
