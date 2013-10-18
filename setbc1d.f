      subroutine setbc1d(nnew,nold,cnew,cold)
c
c-----CAMx v4.02 030709
c
c     SETBC1D sets boundary conditions for the case in which two grids
c     share a common boundary
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        nnew                number of overlapping cells in grid receiving concs
c        nold                number of overlapping cells in grid supplying concs
c        cold                conc to supply (umol/m3)
c
c     Output arguments:
c        cnew                conc received (umol/m3)
c
c     Routines called:
c        none
c
c     Called by:
c        BCMODFY
c
      include "camx.prm"
c
      dimension cnew(nnew),cold(nold)
      dimension xnew(MX1D+1),xold(MX1D+1)
c
c-----Entry point
c
c-----Case 1: nnew = nold
c
      if (nnew.eq.nold) then
        do i=1,nnew
          cnew(i) = cold(i)
        enddo
      endif
c
      do i=1,nnew+1
        xnew(i) = i-1
      enddo
      do i=1,nold+1
        xold(i) = float(nnew)/nold*(i-1)
      enddo
c
c-----Case 2: nnew < nold (or dxnew > dxold)
c
      if (nnew.lt.nold) then
        dxbeg = 0.
        ibeg = 1
        do i1=2,nnew+1
          do i2=1,nold
            if (xnew(i1).gt.xold(i2) .and. xnew(i1).le.xold(i2+1)) then
              iend = i2
              dxend = xnew(i1) - xold(i2)
            endif
          enddo
c
          ctmp = dxbeg*cold(ibeg) + dxend*cold(iend)
          do l=ibeg+1,iend
            ctmp = ctmp + cold(l-1)*(xold(l)-xold(l-1))
          enddo
          cnew(i1-1) = ctmp/(dxbeg+dxend+xold(iend)-xold(ibeg))
          ibeg = iend
          dxbeg = xold(iend) - xnew(i1)
        enddo
      endif
c
c-----Case 3: nnew > nold (or dxnew < dxold)
c
      if (nnew.gt.nold) then
        do i1=2,nold+1
          do i2=2,nnew+1
            if (xnew(i2-1).ge.xold(i1-1) .and. 
     &          xnew(i2).le.xold(i1)) then
              cnew(i2-1) = cold(i1-1)
            endif
            if (xnew(i2-1).lt.xold(i1) .and. xnew(i2).gt.xold(i1)) then
              tmp1 = xnew(i2) - xold(i1)
              tmp2 = xold(i1) - xnew(i2-1)
              w1 = tmp1/(tmp1 + tmp2)
              cnew(i2-1) = cold(i1)*w1 + cold(i1-1)*(1.-w1)
            endif
          enddo
        enddo
      endif
c
      return
      end
