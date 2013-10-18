      subroutine rassgn3d(ncol,nrow,nlay,io,jo,nmesh,ncolf,nrowf,
     &                                                 rcval,rfval)
c
c-----CAMx v4.02 030709
c
c     RASSGN3D assigns real fine grid values from coarse grid. This
c     version is for a 3-D array (rows,columns,layers)
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
c        rcval             cell centered value on coarse grid
c
c     Output arguments:
c        rfval             cell centered value on coarse grid
c
c     Subroutine called:
c        none
c
c     Called by:
c        STARTUP
c        READINP
c
      dimension rcval(ncol,nrow,nlay),rfval(ncolf,nrowf,nlay)
c
c-----Entry point
c
      do ilay=1,nlay
        do 40 jfin = 1,nrowf
          j = (jfin - 2)/nmesh + jo
          do 30 ifin = 1,ncolf
            i = (ifin - 2)/nmesh + io
                rfval(ifin,jfin,ilay) = rcval(i,j,ilay)
  30      continue
  40    continue
      enddo
c
      return
      end
