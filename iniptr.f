      subroutine iniptr()
c
c-----CAMx v4.02 030709
c
c     This routine calculates the pointers into the concentration
c     and meterology arrays for each grid.  The arrays are vectors
c     so this routine calculates the pointer that stores the 
c     first element of the 2-D or 3-D or 4-D field.  Some checks
c     are also made to ensure there is no array overflow of the
c     vectors.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        10/8/99   Fixed a bug in assigning pointers for OSAT arrays for
c                  grids >= 3
c         9/3/02   Removed IPTRCL
c        1/13/03   Added IPTRDP for deposition output fields
c
c     Input arguments:
c        none
c
c     Output arguments:
c        ngcol   -- maximum number of columns for each grid
c        ngrow   -- maximum number of columns for each grid
c        nglay   -- maximum number of columns for each grid
c
c     Subroutine called:
c        none
c
c     Called by:
c        STARTUP
c
c-----Include files
c
      include 'camx.prm'
      include 'grid.com'
      include 'filunit.com'
c
c========================= Source Apportion End ========================
c
      include 'tracer.com'
c
c======================== Source Apportion Begin =======================
c
c
c-----Argument declarations
c
c
c-----Local variables
c
      integer   ngrow(MXGRID), ngcol(MXGRID), nglay(MXGRID)
      integer   ngmax
c
c-----Entry point
c
c   --- load the paramters into local arrays,  this just makes
c       the code for checking paramters a little cleaner ---
c
      ngrow(1) = MXROW1
      ngcol(1) = MXCOL1
      nglay(1) = MXLAY1
c
      if( MXGRID .GT. 1 ) then
         ngrow(2) = MXROW2
         ngcol(2) = MXCOL2
         nglay(2) = MXLAY2
      endif
c
      if( MXGRID .GT. 2 ) then
         ngrow(3) = MXROW3
         ngcol(3) = MXCOL3
         nglay(3) = MXLAY3
      endif
c
      if( MXGRID .GT. 3 ) then
         ngrow(4) = MXROW4
         ngcol(4) = MXCOL4
         nglay(4) = MXLAY4
      endif
c
      if( MXGRID .GT. 4 ) then
         ngrow(5) = MXROW5
         ngcol(5) = MXCOL5
         nglay(5) = MXLAY5
      endif
c
      if( MXGRID .GT. 5 ) then
         ngrow(6) = MXROW6
         ngcol(6) = MXCOL6
         nglay(6) = MXLAY6
      endif
c
      if( MXGRID .GT. 6 ) then
         ngrow(7) = MXROW7
         ngcol(7) = MXCOL7
         nglay(7) = MXLAY7
      endif
c
      if( MXGRID .GT. 7 ) then
         ngrow(8) = MXROW8
         ngcol(8) = MXCOL8
         nglay(8) = MXLAY1
      endif
c
      if( MXGRID .GT. 8 ) then
         ngrow(9) = MXROW9
         ngcol(9) = MXCOL9
         nglay(9) = MXLAY9
      endif
c
      if( MXGRID .GT. 9 ) then
         ngrow(10) = MXROW10
         ngcol(10) = MXCOL10
         nglay(10) = MXLAY10
      endif
c
c   --- check that the parameters for each grid are large
c       enough to hold the grid dimensions --- 
c
      ngmax = -9
      do i=1,ngrid
        ngmax  = MAX(ngrow(i),ngmax)
        ngmax  = MAX(ngcol(i),ngmax)
        ngmax  = MAX(nglay(i),ngmax)
        if( ngrow(i) .LT. nrow(i) .OR. ngcol(i) .LT. ncol(i) 
     &                                .OR. nglay(i) .LT. nlay(i) ) then
           write(iout,'(//,a)') 'ERROR in INIPTR:'
           write(iout,'(1X,A,I2,A)') 'Parameters for grid ',i,
     &                               ' are not large enough.'
           write(iout,'(10X,A,5X,A)') 'Parameters','Grid Definition'
           write(iout,'(A10,5X,I5,10X,I5)') 'Rows    :',ngrow(i),nrow(i)
           write(iout,'(A10,5X,I5,10X,I5)') 'Columns :',ngcol(i),ncol(i)
           write(iout,'(A10,5X,I5,10X,I5)') 'Layers  :',nglay(i),nlay(i)
           write(iout,*) 'Increase the parameters and recompile.'
           call camxerr()
        endif
      enddo
c
c   ---- check the maximum of the 1-Dimensions ----
c
      if( MX1D .LT. ngmax ) then
        write(iout,'(//,a)') 'ERROR in INIPTR:'
        write(iout,*) 'Parameter for maximum of 1-dimensional ',
     &                'arrays is not large enough.'
        write(iout,'(1X,2A,I5,A)') 'Increase parameter MX1D to ',
     &                              'at least: ',ngmax,' and recompile.'   
        call camxerr()
      endif
c
c  ---- everything checks out, set the pointers ---
c
      iptr2d(1) = 1
      iptr3d(1) = 1
      iptr4d(1) = 1
      iptrem(1) = 1
      iptrad(1) = 1
      iptrlu(1) = 1
      iptrdp(1) = 1
      ipsa2d(1) = 1
      ipsa3d(1) = 1
      do i=2,ngrid
         iptr2d(i) = iptr2d(i-1) + ncol(i-1)*nrow(i-1)
         iptr3d(i) = iptr3d(i-1) + ncol(i-1)*nrow(i-1)*nlay(i-1)
         iptr4d(i) = iptr4d(i-1) + ncol(i-1)*nrow(i-1)*nlay(i-1)*MXSPEC
         iptrem(i) = iptrem(i-1) + ncol(i-1)*nrow(i-1)*MXSPEC
         iptrad(i) = iptrad(i-1) + ncol(i-1)*nrow(i-1)*nlay(i-1)*MXRADCL
         iptrlu(i) = iptrlu(i-1) + ncol(i-1)*nrow(i-1)*NLU
         iptrdp(i) = iptrdp(i-1) + ncol(i-1)*nrow(i-1)*MXSPEC*3
         ipsa2d(i) = ipsa2d(i-1) + ncol(i-1)*nrow(i-1)*MXTRSP
         ipsa3d(i) = ipsa3d(i-1) + ncol(i-1)*nrow(i-1)*nlay(i-1)*MXTRSP
      enddo
c
c  --- check for array overflow here, should be caught in above 
c      checks but it can't hurt ----
c
c
      ngmax = iptr2d(ngrid) + ncol(ngrid)*nrow(ngrid) - 1
      if( MXVEC2D .LT. ngmax ) then
         write(iout,'(//,a)') 'ERROR in INIPTR:'
         write(iout,*) 'Parameter for dimensioning vectors for',
     &                 ' 2-D fields is not large enough.'
         write(iout,*) 'Make sure the parameter MXVEC2D is properly ',
     &                 'set and recompile.'
         call camxerr()
      endif
c
      ngmax = iptr3d(ngrid) + ncol(ngrid)*nrow(ngrid)*nlay(ngrid) - 1
      if( MXVEC3D .LT. ngmax ) then
         write(iout,'(//,a)') 'ERROR in INIPTR:'
         write(iout,*) 'Parameter for dimensioning vectors for',
     &                 ' 3-D fields is not large enough.'
         write(iout,*) 'Make sure the parameter MXVEC3D is properly ',
     &                 'set and recompile.'
         call camxerr()
      endif
c
      ngmax = iptr4d(ngrid) + 
     &              ncol(ngrid)*nrow(ngrid)*nlay(ngrid)*MXSPEC - 1
      if( MXVEC4D .LT. ngmax ) then
         write(iout,'(//,a)') 'ERROR in INIPTR:'
         write(iout,*) 'Parameter for dimensioning vectors for',
     &                 ' 4-D fields is not large enough.'
         write(iout,*) 'Make sure the parameter MXVEC2D is properly ',
     &                 'set and recompile.'
         call camxerr()
      endif
c
c----Return point
c
      return
      end
