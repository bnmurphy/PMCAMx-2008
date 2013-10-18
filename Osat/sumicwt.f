      subroutine sumicwt(ncols,nrows,nlays,ibeg,iend,deltax,depth,
     &                   congrd,consum)
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.02 030709
c
c     SUMICWT sums up species concentrations from initial condition
c
c     Input arguments:
c        ncols           the number of coarse grid columns
c        nrows           the number of coarse grid rows
c        nlays           the number of coarse grid layers
c        ibeg            begining row of the domain
c        iend            ending row of the domain
c        depth           layer depth
c        congrd          species concentration for the grid
c
c     Output argumnets:
c        consum          sum of the species concentrations
c
      real   depth(ncols,nrows,nlays),congrd(ncols,nrows,nlays)
      dimension ibeg(nrows),iend(nrows),deltax(nrows)
c
c   --- calculate the average over all vertical cells,
c       weight by layer depth and cell width ---
c
      do 80 j=1,nrows
        if( ibeg(j) .GT. 0 ) then
           do 90 i=ibeg(j),iend(j)
              sumthk = 0
              sumcon = 0
              do 11 izcl=1,nlays
                 sumthk = sumthk + depth(i,j,izcl) * deltax(j)
                 sumcon = sumcon + 
     &              congrd(i,j,izcl) * depth(i,j,izcl) * deltax(j)
   11         continue
              if( sumthk .GT. 0. ) consum = consum + sumcon / sumthk
   90      continue
        endif
   80 continue
c
      return
      end
