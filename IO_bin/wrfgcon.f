      subroutine wrfgcon(iendat,endtim)
c
c-----CAMx v4.02 030709
c
c     WRFGCON writes the instantaneous and average concentration fields
c     for the fine grids.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
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
      include 'flags.com'
      include 'filunit.com'
      include 'camxfld.com'
      include 'chmstry.com'
      include 'grid.com'
c
      integer   iendat
      real      endtim
      integer      ifgptr(MXGRID), ifglvl(MXGRID)
      integer      iunit, igrd, idx, i, j, k, l
      real         cnctmp(MXCOLA,MXROWA)
c
c-----Entry point
c
c-----Set the fine grid pointers
c
      do i = 1,ngrid
        ifglvl(i) = 0
      enddo
      do i = 1,ngrid
        do j = 1,nchdrn(i)
          idch = idchdrn(j,i)
          ifgptr(idch) = i - 1 
          ifglvl(idch) = ifglvl(idch) + 1
          do k = 1,nchdrn(idch)
            ifglvl(idchdrn(k,idch)) = ifglvl(idchdrn(k,idch)) + 1
          enddo
        enddo
      enddo
c
c-----Set the unit number for the instananeous file
c
      if (MOD(INT(endtim/100.),2).EQ.1) then
          iunit = ifconc(1)
      else
          iunit = ifconc(2)
      endif
c
c-----Rewind file and write the file description header
c
      rewind(iunit)
      write(iunit) runmsg
      write(iunit) ngrid-1, nspec
      write(iunit) (spname(i),i=1,nspec)
c
      do 20 igrd = 2,ngrid
        write(iunit) inst1(igrd),jnst1(igrd),inst2(igrd),
     &               jnst2(igrd),meshold(igrd),meshold(igrd),
     &               ncol(igrd),nrow(igrd),nlay(igrd),
     &               ifgptr(igrd),ifglvl(igrd)
   20 continue
c
c-----Write time span
c
      write(iunit) endtim, iendat
      write(ifavg) endtim, iendat
c
c-----Write the concentration data for this hour
c
      do 30 igrd = 2,ngrid
c
        do 40 l = 1,nspec
          do 50 k = 1,nlay(igrd)
            do 60 j = 1,nrow(igrd)
              do 70 i = 1,ncol(igrd)
                idx = i + ncol(igrd)*(j-1) + 
     &                ncol(igrd)*nrow(igrd)*(k-1) +
     &                ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                cnctmp(i,j) = conc(iptr4d(igrd)-1+idx)
   70         continue
   60       continue
            write(iunit) ((cnctmp(i,j),i=1,ncol(igrd)),j=1,nrow(igrd))
   50     continue
   40   continue
c
c
c-----Set the number of layers in the average file
c
        if (l3davg) then
          nlayer = nlay(igrd)
        else
          nlayer = 1
        endif
        do 80 l = 1,navspc
          do 90 k = 1,nlayer
            do 11 j = 1,nrow(igrd)
              do 12 i = 1,ncol(igrd)
                idx = i + ncol(igrd)*(j-1) + 
     &                ncol(igrd)*nrow(igrd)*(k-1) +
     &                ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                cnctmp(i,j) = avcnc(iptr4d(igrd)-1+idx)
   12         continue
   11       continue
            write(ifavg) ((cnctmp(i,j),i=1,ncol(igrd)),j=1,nrow(igrd))
   90     continue
   80   continue
c
   30 continue
c
c-----Close files and return to calling routine
c
      call flush(iunit)
      call flush(ifavg)
c
      return
      end
