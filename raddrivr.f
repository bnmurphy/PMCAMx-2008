      subroutine raddrivr
c
c-----CAMx v4.02 030709
c
c     RADDRIVR initializes radical concentrations for the grids
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:  
c        12/11/00 - initialize radicals to 1.0e-10 - gy  
c
c     Input arguments:
c        none
c
c     Output arguments:
c        none
c
c     Routines called:
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "camx.com"
      include "bndary.com"
      include "camxfld.com"
      include "grid.com"
      include "chmstry.com"
c
c-----Entry point
c
c
c-----Initialize radicals
c
      do igrd = 1,ngrid
        do k = 1,nlay(igrd)
          do j = 2,nrow(igrd)-1
            itmp1 = 2
            itmp2 = ncol(igrd)-1
            if (igrd.eq.1) then
              if (ibeg(j).eq.-999) goto 10
              itmp1=ibeg(j)
              itmp2=iend(j)
            endif
            do i = itmp1,itmp2
              n2d = i + (j-1)*ncol(igrd)
              n3d = n2d + ncol(igrd)*nrow(igrd)*(k-1)
              do l = 1,nrad
                n4d = n3d + ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                cncrad(iptrad(igrd)-1+n4d) = 1.0e-10
              enddo
            enddo
          enddo
  10      continue
        enddo
      enddo
c
      return
      end 
