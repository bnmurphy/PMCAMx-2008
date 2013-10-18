      subroutine rassgn4d(ncol,nrow,nlay,nspec,io,jo,nmesh,nmshv,ncolf,
     &                                         nrowf,nlayf,rcval,rfval)
c
c-----CAMx v4.02 030709
c
c     RASSGN4D assigns real fine grid values from coarse grid. This
c     version is for a 4-D array (rows,columns,layers,species)
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
c        nlay              number of layers in parent grid
c        io                starting i index for the fine grid
c        jo                starting j index for the fine grid
c        nmesh             mesh number
c        ncolf             number of columns in fine grid
c        nrowf             number of rows in fine grid
c        nlayf             number of layers in fine grid
c        nmshv             vertical mesh number
c        rcval             cell centered value on coarse grid
c
c     Output arguments:
c        rfval             cell centered value on coarse grid
c
c     Subroutine called:
c        EXPNDLAY
c
c     Called by:
c        RDGGCON
c
      integer   nmshv(nlay)
      dimension rcval(ncol,nrow,nlay,nspec)
      dimension rfval(ncolf,nrowf,nlayf,nspec)
c
c-----Entry point
c
c----Loop over species---
c
      do ispc=1,nspec
         kfin = 0
c
c----Loop over layers in parent grid---
c
         do ilay=1,nlay
c
c----Loop over fine grid layers in this parent layer---
c
           do ilayf=1,nmshv(ilay)
c
c----Update pointer to fine grid layer---
c
              kfin = kfin + 1
c
c----Load in the grid---
c
              do 40 jfin = 1,nrowf
                j = (jfin - 2)/nmesh + jo
                do 30 ifin = 1,ncolf
                   i = (ifin - 2)/nmesh + io
                   rfval(ifin,jfin,kfin,ispc) = rcval(i,j,ilay,ispc)
  30            continue
  40          continue
c
c----Next layer---
c
           enddo
         enddo
c
c----Next species----
c
      enddo
c
      return
      end
