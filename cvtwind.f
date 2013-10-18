      subroutine cvtwind(ncol,nrow,nlay,windu,windv)
c
c-----CAMx v4.02 030709
c
c     CVTWIND converts the wind field from cell-centered to staggered on cell
c     interfaces (Arakawa C grid):
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        windu               cell centered windspeed in x-direction (m/s)
c        windv               cell centered windspeed in y-direction (m/s)
c
c     Output arguments:
c        windu               staggered windspeed in x-direction (m/s)
c        windv               staggered windspeed in y-direction (m/s)
c
c     Routines Called:
c        none
c
c     Called by:
c        READZPWT
c        READINP
c
      dimension windu(ncol,nrow,nlay),windv(ncol,nrow,nlay)
c
      do 10 k=1,nlay
c
c-----x-direction
c
        do j=1,nrow
          do i=1,ncol-1
            windu(i,j,k) = (windu(i,j,k) + windu(i+1,j,k))/2.
          enddo
        enddo
c
c-----y-direction
c
        do j=1,nrow-1
          do i=1,ncol
            windv(i,j,k) = (windv(i,j,k) + windv(i,j+1,k))/2.
          enddo
        enddo
c
  10  continue
c
      return
      end
