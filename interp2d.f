      subroutine interp2d(ncol,nrow,nlay,io,jo,nmesh,ncolf,nrowf,
     &                    cval,fval)
c
c-----CAMx v4.02 030709
c
c     INTERP2D horizontally interpolates a coarse grid field to a
c     fine grid
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        none
c
c     Input arguments:
c        ncol              number of columns in parent grid
c        nrow              number of rows in parent grid 
c        nlay              number of layers in parent grid
c        io                starting i index for the fine grid
c        jo                starting j index for the fine grid
c        nmesh             mesh number
c        ncolf             number of columns in fine grid
c        nrowf             number of rows in fine grid
c        cval              cell centered value on coarse grid
c
c     Output arguments:
c        fval              cell centered value on fine grid
c
c     Subroutine called:
c        none
c
c     Called by:
c        STARTUP
c        INTRPDAT
c        INTRPCNC
c
      dimension cval(ncol,nrow,nlay),fval(ncolf,nrowf,nlay)
c
c-----Entry point
c
      do 50 k = 1,nlay
        jc = jo
        yc = float(jc - jo + 1) - 0.5
        do j = 1,nrowf
          ic = io
          xc = float(ic - io + 1) - 0.5
          yf = (j - 1.5)/float(nmesh)
          if (yc.le.yf) then
              jc = jc + 1
              yc = float(jc - jo + 1) - 0.5
          endif
          do i = 1,ncolf
            xf = (i - 1.5)/float(nmesh)
            if (xc.le.xf) then
              ic = ic + 1
              xc = float(ic - io + 1) - 0.5
            endif
            dcdx1 = (cval(ic,jc-1,k) - cval(ic-1,jc-1,k))
            dcdx2 = (cval(ic,jc,k) - cval(ic-1,jc,k))
            c1 = cval(ic,jc-1,k) - dcdx1*(xc - xf)
            c2 = cval(ic,jc,k) - dcdx2*(xc - xf)
            dcdy = (c2 - c1)
            fval(i,j,k) = c2 - dcdy*(yc - yf)
          enddo
        enddo
c
  50  continue
c
      return
      end
