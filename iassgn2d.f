      subroutine iassgn2d(ncol,nrow,io,jo,nmesh,ncolf,nrowf,icval,ifval)
c
c-----CAMx v4.02 030709
c
c     IASSGN2D assigns integer fine grid values from coarse grid
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        none
c
c     Input arguments:
c        ncol              number of columns in the parent grid
c        nrow              number of rows in the parent grid
c        io                starting i index for the fine grid
c        jo                starting j index for the fine grid
c        nmesh             mesh number
c        ncolf             number of columns in fine grid
c        nrowf             number of rows in fine grid
c        icval             cell centered value on coarse grid
c
c     Output arguments:
c        ifval             cell centered value on coarse grid
c
c     Subroutine called:
c        none
c
c     Called by:
c        CAMx
c        STARTUP
c
      dimension icval(ncol,nrow),ifval(ncolf,nrowf)
c
c-----Entry point
c
        do 40 jfin = 2,nrowf-1
          j = (jfin - 2)/nmesh + jo
          do 30 ifin = 2,ncolf-1
            i = (ifin - 2)/nmesh + io
            ifval(ifin,jfin) = icval(i,j)
  30      continue
  40    continue
c
      return
      end
