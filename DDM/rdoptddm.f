c**** RDOPTDDM
c
      subroutine rdoptddm(iounit,ctlfil,irec)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c      This routine reads all of the user options and flags for the
c      source apportionment algorithm.  This version is for the
c      DDM technology type, which has different options.
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
c     03/23/99   --gwilson--    Original development
c     10/20/00   --cemery --    Split read of LSFCFL into 2 separate records
c     07/19/01   --gyarwood-    Initialize ngroup, nregin, lbndry
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'filunit.com'
      include 'tracer.com'
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
      character*80 inprec, action
      integer      nlines
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. lddm ) goto 9999
c
c   --- flags for opening coarse and fine grid output surface files
c
      irec = irec + 1
      action = 'Reading coarse grid output flag.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),*,ERR=7003) lsfcfl(1)
      irec = irec + 1
      action = 'Reading fine grid output flag.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),*,ERR=7003) lsfcfl(2)
c
c   --- get the number of Initial conditions species ---
c
      irec = irec + 1
      action = 'Reading the number of IC species.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),*,ERR=7003) nicddm
      if( nicddm .GT. 0 ) then
         if( MOD(nicddm,6) .EQ. 0 ) then
            nlines = nicddm/6
         else
            nlines = INT(nicddm/6) + 1
         endif
         do i=1,nlines
           irec = irec + 1
           ibeg = (i - 1) * 6 + 1
           action = 'Reading IC species names.'
           read(iounit,'(A)',ERR=7002) inprec 
           read(inprec(21:),'(6A10)',ERR=7002) 
     &                                 (icddmsp(l),l=ibeg,ibeg+5)
         enddo
         do i=1,nicddm
           call toupper( icddmsp(i) )
           call jstlft( icddmsp(i) )
         enddo
      endif
c
c   --- get the number of boundary conditions species ---
c
      lbndry = .FALSE.
      irec = irec + 1
      action = 'Reading the number of BC species.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),*,ERR=7003) nbcddm
      if( nbcddm .GT. 0 ) then
         if( MOD(nbcddm,6) .EQ. 0 ) then
            nlines = nbcddm/6
         else
            nlines = INT(nbcddm/6) + 1
         endif
         do i=1,nlines
           irec = irec + 1
           ibeg = (i - 1) * 6 + 1
           action = 'Reading BC species names.'
           read(iounit,'(A)',ERR=7002) inprec 
           read(inprec(21:),'(6A10)',ERR=7003) 
     &                                 (bcddmsp(l),l=ibeg,ibeg+5)
         enddo
         do i=1,nbcddm
           call toupper( bcddmsp(i) )
           call jstlft( bcddmsp(i) )
         enddo
c
c   --- flag for stratifying the boundary by edge ---
c
         irec = irec + 1
         action = 'Reading flag for stratified BCs.'
         read(iounit,'(A)',ERR=7002) inprec
         read(inprec(21:),'(L10)',ERR=7003) lbndry
      endif
c
c   --- get the number of emissions species ---
c
      ngroup = 0
      nregin = 0
      irec = irec + 1
      action = 'Reading the number of emissions species.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),*,ERR=7003) nemddm
c
c   --- read records that depend on nemddm ---
c
      if( nemddm .GT. 0 ) then
         if( MOD(nemddm,6) .EQ. 0 ) then
            nlines = nemddm/6
         else
            nlines = INT(nemddm/6) + 1
         endif
         do i=1,nlines
           irec = irec + 1
           ibeg = (i - 1) * 6 + 1
           action = 'Reading emissions species names.'
           read(iounit,'(A)',ERR=7002) inprec
           read(inprec(21:),'(6A10)',ERR=7003)
     &                                 (emddmsp(l),l=ibeg,ibeg+5)
         enddo
         do i=1,nemddm
           call toupper( emddmsp(i) )
           call jstlft( emddmsp(i) )
         enddo
c
c   --- number of source regions ---
c
        irec = irec + 1
        action = 'Reading number of emissions source regions.'
        read(iounit,'(A)',ERR=7002) inprec
        read(inprec(21:),'(I10)',ERR=7003) nregin
        if( nregin .EQ. 0 ) then
           write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
           write(iout,'(/,1X,2A,/,9X,A)') 'ERROR:  When requesting ',
     &                'emissions species for DDM you must ',
     &                     'provide at least 1 emissions region'
           call camxerr()
        endif
c
c   --- number of source emissions groupings ---
c
        irec = irec + 1
        action = 'Reading number of emissions source groups.'
        read(iounit,'(A)',ERR=7002) inprec
        read(inprec(21:),'(I10)',ERR=7003) ngroup
        if( ngroup .EQ. 0 ) then
           write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
           write(iout,'(/,1X,2A,/,9X,A)') 'ERROR:  When requesting ',
     &                'emissions species for DDM you must ',
     &                     'provide at least 1 emissions group.'
           call camxerr()
        endif
      endif
c
c   --- leftover is always false with DDM ---
c
      leftovr = .FALSE.
c
c   --- check for array overflow ---
c
      if( ngroup .GT. MXTEMF-1 ) then
         write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
         write(iout,'(/,1X,A,I4,A)') 
     &                  'ERROR: Number of source groupings ',ngroup,
     &                            ' exceeds maximum.  Increase MXTEMF.'
         call camxerr()
      endif
      if( ngroup .LT. 0 ) then
         write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
         write(iout,'(1X,A,I4,A)') 
     &                  'ERROR: Number of emissions groups ',ngroup,
     &                                                 ' is invalid.'
         call camxerr()
      endif
c
c   --- number of timing realeases is always zero for DDM ---
c
      ntrtim = 0
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
      write(iout,*) 'ERROR in RDOPTDDM reading control file (',
     &                     ctlfil(:istrln(ctlfil)),') at line: ',irec
      write(iout,*) action(:istrln(action))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in RDOPTDDM:'
      write(iout,*) 'ERROR in RDOPTDDM reading control file (',
     &                     ctlfil(:istrln(ctlfil)),') at line: ',irec
      write(iout,*) action(:istrln(action))
      write(iout,*) '     Last line read:'
      write(iout,*) inprec(:istrln(inprec))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
