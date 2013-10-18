      subroutine wrfgdep(iendat,endtim)
c
c-----CAMx v4.02 030709
c
c     WRFGDEP writes the deposition fields for the fine grids.
c
c     Copyright 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        iendat              current ending date (YYJJJ)
c        endtim              current ending time (HHMM)
c
c     Output arguments:
c        none
c
c     Rountines Called:
c        none
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'camxfld.com'
      include 'chmstry.com'
      include 'grid.com'
c
      integer   iendat
      real      endtim
      integer   igrd, idx, i, j, l
      real      deptmp(MXCOLA,MXROWA)
c
c-----Entry point
c
c-----Write time span
c
      write(ifdep) endtim, iendat
c
c-----Write the concentration data for this hour
c
      do 30 igrd = 2,ngrid
c
        do l = 1,navspc
          ll = lavmap(l)
          do j = 1,nrow(igrd)
            do i = 1,ncol(igrd)
              idx = i + ncol(igrd)*(j-1) + ncol(igrd)*nrow(igrd)*(ll-1)
              deptmp(i,j) = vdep(iptrem(igrd)-1+idx)
            enddo
          enddo
          write(ifdep) ((deptmp(i,j),i=1,ncol(igrd)),j=1,nrow(igrd))
        enddo
        do 80 l = 1,3*navspc
          do 11 j = 1,nrow(igrd)
            do 12 i = 1,ncol(igrd)
              idx = i + ncol(igrd)*(j-1) + ncol(igrd)*nrow(igrd)*(l-1)
              deptmp(i,j) = depfld(iptrdp(igrd)-1+idx)
   12       continue
   11     continue
          write(ifdep) ((deptmp(i,j),i=1,ncol(igrd)),j=1,nrow(igrd))
   80   continue
c
   30 continue
c
c-----Close files and return to calling routine
c
      call flush(ifdep)
c
      return
      end
