c**** RDOPTSA
c
      subroutine rdoptsa(iounit,ctlfil,irec)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c      This routine reads all of the user options and flags for the
c      source apportionment algorithm.  
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
c     05/16/96   --gwilson--    Original development
c     02/10/97   --cemery--     Added read of fileroot for SA output files
c     04/14/97   --gwilson--    Changed the way number of source groups
c                               is specified
c     04/28/97   --gwilson--    Added flag for OPPAT
c     05/22/97   --gwilson--    Added flag for APCA
c     10/20/00   --cemery --    Split read of LSFCFL into 2 separate records
c     01/10/02   --cemery --    Eliminated the read for fine grid
c                               flag if no nests 
c     03/21/03   --gwilson--    Removed the OSAT technology type OPPAT
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'filunit.com'
      include 'grid.com'
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
      integer      nemiss
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. ltrace ) goto 9999
c
c   --- flags for opening coarse and fine grid output surface files
c
      irec = irec + 1
      action = 
     &   'Reading flag for coarse grid output for source apportionment.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),*,ERR=7003) lsfcfl(1)
      if( ngrid .GT. 1 ) then
          irec = irec + 1
          action = 
     &     'Reading flag for fine grid output for source apportionment.'
          read(iounit,'(A)',ERR=7002) inprec
          read(inprec(21:),*,ERR=7003) lsfcfl(2)
      endif
c
c   --- flag for stratifying the boundary by edge ---
c
      irec = irec + 1
      action = 
     &  'Reading flag for stratifying boundary in source apportionment.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),'(L10)',ERR=7003) lbndry
c
c   --- number of source regions ---
c
      irec = irec + 1
      action = 
     &     'Reading number of source regions in source apportionment.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),'(I10)',ERR=7003) nregin
c
c   --- number of source emissions groupings ---
c
      irec = irec + 1
      action = 
     &     'Reading number of emissions groups in source apportionment.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),'(I10)',ERR=7003) nemiss
      if( tectyp .EQ. GOAT .AND. nemiss .NE. 1 ) goto 7001
c
c   --- flags for determining if the leftover group should be used ---
c
      irec = irec + 1
      action = 
     &      'Reading flag for leftover group in source apportionment.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),'(L10)',ERR=7003) leftovr
c
c   --- set the number of source groups from the emissions groups ---
c
      if( nemiss .EQ. 1 .AND. leftovr ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,2A)') 'Cannot have leftover group ',
     &                           'with only one source group.'
         write(iout,'(1X,2A)') 'Set number of groups to 2 or turn ',
     &                         'off leftover group.'
         call camxerr()
      endif
      if( leftovr .OR. nemiss .EQ. 1 ) then
          ngroup = nemiss - 1
      else
          ngroup = nemiss 
      endif
c
c   --- Need the biogenics group if doing APCA ---
c
      if( tectyp .EQ. APCA .AND. ngroup .EQ. 0 ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,3A)')'Need biogenic sources as a seperate ',
     &                          'emissions group when doing APCA.'
         call camxerr()
      endif
c
c   --- check for array overflow ---
c
      if( ngroup .GT. MXTEMF-1 ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(/,1X,A,I4,A)')'Number of source groupings ',
     &                     ngroup,' exceeds maximum.  Increase MXTEMF.'
         call camxerr()
      endif
      if( ngroup .LT. 0 ) then
         write(iout,'(//,a)') 'ERROR in RDOPTSA:'
         write(iout,'(1X,A,I4,A)')'Number of emissions groups ',nemiss,
     &                            ' is invalid.'
         call camxerr()
      endif
c
c   --- number of timing releases per day ----
c
      irec = irec + 1
      action = 
     &      'Reading number of timing tracers in source apportionment.'
      read(iounit,'(A)',ERR=7002) inprec
      read(inprec(21:),'(I10)',ERR=7003) ntrtim
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
      write(iout,'(//,a)') 'ERROR in RDOPTSA:'
      write(iout,'(/,1X,2A)') 'The technology type option: ',
     &               ' GOAT requires that there be 1 emissions group.'
      write(iout,'(10X,A,I3,A)') 'You supplied: ',nemiss,
     &                                            ' emissions groups.' 
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDOPTSA:'
      write(iout,*) 'Reading control file ( ',
     &                     ctlfil(:istrln(ctlfil)),' ) at line: ',irec
      write(iout,'(A)') action(:istrln(action))
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in RDOPTSA:'
      write(iout,*) 'Reading control file ( ',
     &                     ctlfil(:istrln(ctlfil)),' ) at line: ',irec
      write(iout,'(A)') action(:istrln(action))
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
