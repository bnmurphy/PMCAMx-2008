c**** RDRCPRT
c
      subroutine rdrcprt()
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the receptor definition file and stores all of
c   the receptor data in common arrays to be used by the RTRAC
c   version of the CAMx.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Argument description:
c      Inputs:
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/27/96   --gwilson--    Original development
c     04/02/03   --gwilson--    Added grid number to recptors defined by
c                               cell index
c         
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'tracer.com'
      include 'rtracchm.com'
      include 'procan.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 line
      character*15  rcptmp
      integer       irec, ierr, idtmp, icel, jcel, igrd, n2d
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      irec = 0
      lwalls = .FALSE.
c
c  --- open the file ---
c
      open(unit=IORRCP,file=rcpfil,ERR=7000,status='UNKNOWN')
c
c  --- initialize the number of receptors ----
c
      nrcprt = 0
c
c  --- read the data as character strings since we don't know what
c      type of receptor it is ---
c
  111 continue
      irec = irec + 1
      read(IORRCP,'(A)',IOSTAT=ierr,END=222) line
      if( ierr .NE. 0 ) then
         write(iout,'(//,a)') 'ERROR in RDRCPRT:'
         write(iout,'(/,1X,A,A,I5)') 'Reading the receptor ',
     &                    'definition file for RTRAC at record: ',irec
         call camxerr()
      endif
c
c  ---- find out what kind of receptor is being provided ---
c
      rcptmp = line(1:15)
      if( rcptmp .EQ. CDCEL ) then
          idtmp = IDCEL 
      else
          write(iout,'(//,a)') 'ERROR in RDRCPRT:'
          write(iout,'(/,1X,A,A,A)') 'Invalid receptor ',
     &       'code found in RTRAC receptor definition file: ',idtmp
          write(iout,'(1X,A,A)') 'Line read: ',line
          write(iout,'(1X,A)') 'Code must be one of the following:'
          write(iout,'(10X,A)') CDCEL
          call camxerr()
      endif
c
c   --- if it is a SINGLE CELL type just read the cell indexes ---
c
      if( idtmp .EQ. IDCEL ) then
          read(line(31:60),'(3I10)',IOSTAT=ierr) igrd, icel, jcel
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RDRCPRT:'
             write(iout,'(/,1X,A,A,A)') 'Reading cell indexes ',
     &                      'from receptor definition file.  Line read:'
             write(iout,'(A)') line(:istrln(line))
             call camxerr()
          endif
c
c  --- check for validity of receptor ---
c
          if( igrd .LT. 0 .OR. igrd .GT. ngrid ) then
             write(iout,'(//,A)') 'ERROR in RDRCPRT:'
             write(iout,'(1X,2A)') 
     &          'Invalid grid index specified in receptor definition.',
     &                                                  '  Line read: '
             write(iout,'(A)') line(:istrln(line))
          endif
          if( icel .LE. 0 .OR. icel .GT. ncol(igrd) .OR. 
     &                    jcel .LE. 0 .OR. jcel .GT. nrow(igrd) ) then
             write(iout,'(//,a)') 'ERROR in RDRCPRT:'
             write(iout,'(1X,A)') 
     &          'Invalid cell index specified in receptor definition.'
             write(iout,'(1X,2A)') 'Cell is outside specified domain. ',
     &                                                    'Line read: '  
             write(iout,'(A)') line(:istrln(line))
          endif
          n2d = icel + (jcel-1)*ncol(igrd)
          if( idfin( iptr2d(igrd)-1+n2d ) .GT. 0 ) then
             write(iout,'(//,a)') 'ERROR in RDRCPRT:'
             write(iout,'(1X,A)') 
     &          'Invalid cell index specified in receptor definition.'
             write(iout,'(1X,2A)')'Cell cannot be contained in a nest.',
     &                                                    'Line read: ' 
             write(iout,'(A)') line(:istrln(line))
          endif
c
c  --- cell is valid, fill array with entire column ---
c
          do k=1,nlay(igrd)
             nrcprt = nrcprt + 1
             if( nrcprt .GT. MXRTCEL ) then
                write(iout,'(//,a)') 'ERROR in RDRCPRT:'
                write(iout,'(/,1X,A,A,I5)') 'Number of receptors in ',
     &              'receptor definition file exceeds maximum: ',MXRTCEL
                write(iout,'(1X,2A)') 'Note that each specified ',
     &                'receptor represents an entire column of cells.'
                write(iout,'(1X,A)') 'Increase the paramater MXRTCEL.'
                call camxerr()
             endif
             n3d = icel + ncol(igrd)*(jcel-1) + 
     &                               ncol(igrd)*nrow(igrd)*(k-1)
             ipacl_3d ( iptr3d(igrd)-1+n3d ) = nrcprt
             idomrt(nrcprt) = igrd
             ircprt(nrcprt) = icel
             jrcprt(nrcprt) = jcel
             krcprt(nrcprt) = k
          enddo
      endif
c      
c  --- get next record ---
c
      goto 111
c      
c  --- close file and return to calling routine ---
c
  222 continue
      close(IORRCP)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDRCPRT:'
      write(iout,'(/,1X,2A)') 'Cannot open input file: ',
     &                                         rcpfil(:istrln(rcpfil))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
