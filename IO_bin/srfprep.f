      subroutine srfprep(igrd,ncol,nrow,fsurf)
c
c-----CAMx v4.03 031205
c
c     SRFPREP reads the landuse files for all grids and initializes the
c     landuse field arrays.  Landuse is mapped to all nested grids that
c     are not supplied with a surface file
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        6/6/03       check input data for consistency
c        7/7/03       now skips everthing if file not supplied
c
c     Input arguments:
c        igrd                grid index
c        ncol                number of columns
c        nrow                number of rows
c
c     Output arguments:
c        fsurf               fractional landuse field
c
c     Routines Called:
c        none
c
c     Called by:
c        STARTUP
c
      include 'camx.prm'
      include 'filunit.com'
      include 'bndary.com'
c
      dimension fsurf(ncol,nrow,NLU)
c
c-----Entry point
c
c  --- skip if file not provided ---
c
      iunit = isurf(igrd)
      if( iunit .LE. 0 ) goto 9999 
c
c  --- read the data, all in one record ---
c
      read(iunit) (((fsurf(i,j,l),i=1,ncol),j=1,nrow),l=1,NLU)
c
c-----Check that the data are reasonable and adjust out
c     minor inconsistencies
c
      do 10 j = 2,nrow-1
         i1 = 2
         i2 = ncol - 1
         if( igrd .eq. 1 ) then
           if( ibeg(j) .EQ. -999 ) goto 10
           i1 = ibeg(j)
           i2 = iend(j)
         endif
         do i = i1,i2
           areatot = 0.0
           do l = 1,NLU
             areatot = areatot + fsurf(i,j,l)
           enddo
           if (areatot.lt.0.95 .or. areatot.gt.1.05) goto 900
           do l = 1,NLU
             fsurf(i,j,l) = fsurf(i,j,l)/areatot
           enddo
         enddo
 10   continue
      return
c
c-----Error in landuse fractions
c
 900  write(iout,'(//,A)') 'ERROR in SRFPREP'
      write(iout,'(/,2A)') 'Sum of landuse fractions differs from 1.0 ',
     &                     'by more than 5%'
      write(iout,'(A,i3,a,2i4)') 
     &       'Grid = ', igrd, '  Cell(i,j) = ', i, j
      write(iout,'(A)') 'Table of input landuse data follows:'
      write(iout,'(/,A)') ' Class    Fraction'
      write(iout,'(i6,F10.3)') (l,fsurf(i,j,l),l=1,NLU)
      write(iout,'(A6,F10.3)') 'Total', areatot
      write(iout,'(/,2A)') 'You can use the SRFLND program to check',
     &                     ' landuse data'
      call camxerr()
c
 9999 continue
      return
      end
