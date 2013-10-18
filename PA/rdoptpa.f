c**** RDOPTPA
c
      subroutine rdoptpa(iounit,ctlfil,irec)
c  
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c      This routine reads all of the user options and flags for the
c      process analysis algorithm.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c       Argument description:      
c           iounit  I  unit number for input
c           ctlfil  C  name of file being read
c           irec    I  current record number (input and output)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     06/08/00   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'filunit.com'
      include 'grid.com'
      include 'procan.com'
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Argument declaration:
c-----------------------------------------------------------------------
c
      integer       iounit
      character*(*) ctlfil
      integer       irec
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*80 inprec
      integer      igrd, nox, noy
      logical      lerror
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. lproca ) goto 9999
c
c   ---- Read Process Analysis domains
c
      irec = irec + 1
      read(iounit,'(A)',ERR=7001) inprec 
      read(inprec(21:),'(I10)',ERR=7002) npadom
c
c   ---- Check for array overflow ----
c
      if( npadom .GT. MXPADOM ) goto 7003
      if( npadom .LE. 0 ) goto 7004
c
c   ---- Read the domain definitions ---
c
      do n = 1,npadom
c
c   --- the grid number ---
c
         irec = irec + 1
         read(iounit,'(A)',ERR=7001) inprec 
         read(inprec(21:),*,ERR=7002) ipagrd(n)
         if( ipagrd(n) .LE. 0 .OR. ipagrd(n) .GT. ngrid ) goto 7005
c
c   --- the cell indexes in the X-driection
c
         irec = irec + 1
         read(iounit,'(A)',ERR=7001) inprec 
         read(inprec(21:),*,ERR=7002) i_sw(n), i_ne(n)
c
c   --- the cell indexes in the Y-driection
c
         irec = irec + 1
         read(iounit,'(A)',ERR=7001) inprec 
         read(inprec(21:),*,ERR=7002) j_sw(n), j_ne(n)
c
c   --- the cell indexes in the Z-driection
c
         irec = irec + 1
         read(iounit,'(A)',ERR=7001) inprec 
         read(inprec(21:),*,ERR=7002) b_lay(n), t_lay(n)
c
c   --- check for valid cell indexes ---
c
         lerror = .FALSE.
         igrd = ipagrd(n)
         if( igrd .GT. 1 ) then
            nox = (inst2(igrd) - inst1(igrd) + 1 ) * meshold(igrd) + 2
            noy = (jnst2(igrd) - jnst1(igrd) + 1 ) * meshold(igrd) + 2
         else
            nox = ncol(1)
            noy = nrow(1)
         endif
         if( i_sw(n) .LE. 0 .OR. i_sw(n) .GT. i_ne(n) ) lerror = .TRUE. 
         if( i_ne(n) .GT. nox ) lerror = .TRUE. 
         if( j_sw(n) .LE. 0 .OR. j_sw(n) .GT. j_ne(n) ) lerror = .TRUE. 
         if( j_ne(n) .GT. noy ) lerror = .TRUE. 
         if( b_lay(n) .LE. 0 .OR. b_lay(n) .GT. t_lay(n) ) 
     &                                                  lerror = .TRUE.
         if( t_lay(n) .GT. nlay(igrd) ) lerror = .TRUE.
         if( lerror ) goto 7006
      enddo
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDOPTPA:'
      write(iout,'(1X,3A,I3)') 'Reading control file: ',
     &                     ctlfil(:istrln(ctlfil)),' at line: ',irec
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDOPTPA:'
      write(iout,'(1X,3A,I3)') 'Reading control file: ',
     &                     ctlfil(:istrln(ctlfil)),' at line: ',irec
      write(iout,*) '     Last line read:'
      write(iout,*) inprec(:istrln(inprec))
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in RDOPTPA:'
      write(iout,'(1X,3A,I3)') 'Reading control file: ',
     &                     ctlfil(:istrln(ctlfil)),' at line: ',irec
      write(iout,'(1X,A,I3,2A)') 
     &           'Number of Process Analysis domains selected: ',
     &                                        npadom,' exceeds max. '
      write(iout,'(1X,A)') 'Increase the parameter MXPADOM.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in RDOPTPA:'
      write(iout,'(1X,3A,I3)') 'Reading control file: ',
     &                     ctlfil(:istrln(ctlfil)),' at line: ',irec
      write(iout,'(1X,2A,I3)') 'Invalid value for number of ',
     &                         'Process Annalysis domains: ',npadom
      write(iout,*) '     Line read:'
      write(iout,*) inprec(:istrln(inprec))
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in RDOPTPA:'
      write(iout,'(1X,3A,I3)') 'Reading control file: ',
     &                     ctlfil(:istrln(ctlfil)),' at line: ',irec
      write(iout,'(1X,2A,I3)') 'Invalid value for grid number for ',
     &                            'Process Analysis domain: ',ipagrd(n)
      write(iout,*) '     Line read:'
      write(iout,*) inprec(:istrln(inprec))
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in RDOPTPA:'
      write(iout,'(1X,3A,I3)') 'Reading control file: ',
     &                     ctlfil(:istrln(ctlfil)),' at line: ',irec
      write(iout,'(1X,A,I3,A)') 'Cell indexes for PA sub-domain: ',n,
     &                                                 ' are invalid.'
      write(iout,'(10X,A,I3,10X,A)') 'Sub-Domain #',n,'Grid Definition'
      write(iout,'(14X,2A,18X,A,I3)') 'Beg ',' End','Grid #',igrd
      write(iout,'(A10,1X,2I5,20X,I5)') 'Rows    :',i_sw(n),i_ne(n),nox
      write(iout,'(A10,1X,2I5,20X,I5)') 'Columns :',j_sw(n),j_ne(n),noy
      write(iout,'(A10,1X,2I5,20X,I5)') 'Layers  :',
     &                                    b_lay(n),t_lay(n),nlay(igrd)
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
