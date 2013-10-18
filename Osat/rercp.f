c**** RERCP
c
      subroutine rercp()
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the receptor definition file and stores all of
c   the receptor data in common arrays to be used by the OSAT
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
      integer       irec, ierr, idtmp, icl, jcl, icrs, jcrs
      integer       i
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
c  --- initialize the hourly peak receptor ----
c
      nrecep = 1
      idrcp(nrecep) = IDCEL
      rcpnam(nrecep) = 'CoarsePeak'
c
c  --- read the data as character strings since we don't know what
c      type of receptor it is ---
c
  111 continue
      irec = irec + 1
      read(IORRCP,'(A)',IOSTAT=ierr,END=222) line
      if( ierr .NE. 0 ) then
         write(iout,'(//,a)') 'ERROR in RERCP:'
         write(iout,'(/,1X,2A,I5)') 'Reading the receptor ',
     &                    'definition file for OSAT at record: ',irec
         call camxerr()
      endif
c
c  ---- find out what kind of receptor is being provided ---
c
      rcptmp = line(1:15)
      if( rcptmp .EQ. CDPNT ) then
          idtmp = IDPNT
      else if( rcptmp .EQ. CDCEL ) then
          idtmp = IDCEL 
      else if( rcptmp .EQ. CDAVG ) then
          idtmp = IDAVG
      else if( rcptmp .EQ. CDWAL ) then
          idtmp = IDWAL
      else
          write(iout,'(//,a)') 'ERROR in RERCP:'
          write(iout,'(/,1X,3A)') 'Invalid receptor ',
     &                  'code found in receptor definition file: ',idtmp
          write(iout,'(1X,2A)') 'Line read: ',line
          write(iout,'(1X,A)') 'Code must be one of the following:'
          write(iout,'(10X,A)') CDPNT
          write(iout,'(10X,A)') CDCEL
          write(iout,'(10X,A)') CDAVG
          write(iout,'(10X,A)') CDWAL
          call camxerr()
      endif
c
c   --- incrament counter and check for overflow ---
c
      nrecep = nrecep + 1
      if( nrecep .GT. MXRECP ) then
          write(iout,'(//,a)') 'ERROR in RERCP:'
          write(iout,'(/,1X,2A,I5)') 'Number of receptors in ',
     &              'receptor definition file exceeds maximum: ',MXRECP
          write(iout,'(1X,A)') 'Increase the paramater MXRECP.'
          call camxerr()
      endif
c
c   --- if it is a POINT type of receptor, just read the coordinates ---
c
      if( idtmp .EQ. IDPNT ) then
          idrcp(nrecep) = IDPNT
          rcpnam(nrecep) = line(21:30)
          read(line(31:50),'(2F10.0)',IOSTAT=ierr) recepx(nrecep), 
     &                                                   recepy(nrecep)
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Reading coordinates ',
     &                    'from receptor definition file.  Line read: '
             write(iout,'(A)') line
             call camxerr()
          endif
          if( recepy(nrecep) .GE. yorg ) then
              jcrs = INT( (recepy(nrecep) - yorg) / 
     &                                   (deltay(1)/1000.0) ) + 1
          else 
              jcrs = -1
          endif
          if( jcrs .GT. 0 .AND. jcrs .LE. nycell(1) .AND. 
     &                               recepx(nrecep) .GE. xorg ) then
             icrs = INT( (recepx(nrecep) - xorg) / 
     &                                   (deltax(jcrs,1)/1000.0) ) + 1
          else
              icrs = -1
          endif
          if( icrs .LE. 0 .OR. icrs .GT. nxcell(1)
     &                 .OR. jcrs .LE. 0 .OR. jcrs .GT. nycell(1) ) then
             write(iout,'(//,a)') 'WARNING in RERCP:'
             write(iout,'(1X,3A,F8.2,A,F8.2,A)')
     &                  'WARNING: Cell for receptor: ',rcpnam(nrecep),
     &                      'not in domain: (',
     &                       recepx(nrecep),',',recepy(nrecep),')'
             write(iout,'(1X,A)') '   Skipping receptor ...'
             nrecep = MAX0(0,nrecep-1)
          endif
c
c  --- find the grid in which this point is contained ---
c
          igrd = 1
          icl = icrs
          jcl = jcrs
          do ig = 2,ngrid
            if(icrs .GE. inst1(ig) .AND. icrs .LE. inst2(ig) .AND.
     &             jcrs .GE. jnst1(ig) .AND. jcrs .LE. jnst2(ig)) then
               xloctmp = (recepx(nrecep)-xorg) - (inst1(ig)-1)*delx
               yloctmp = (recepy(nrecep)-yorg) - (jnst1(ig)-1)*dely
               if( xloctmp .GT. 0. .AND. yloctmp .GT. 0. ) then
                  ii = 2 + INT( xloctmp/delx * FLOAT( meshold(ig) ) )
                  jj = 2 + INT( yloctmp/dely * FLOAT( meshold(ig) ) )       
                  igrd = ig
                  icl = ii
                  jcl = jj
               endif
             endif
          enddo
          igrdrcp(nrecep) = igrd
          irecep(nrecep,1) = icl
          jrecep(nrecep,1) = jcl
      endif
c
c   --- if it is a SINGLE CELL type just read the cell indexes ---
c
      if( idtmp .EQ. IDCEL ) then
          idrcp(nrecep) = IDCEL
          rcpnam(nrecep) = line(21:30)
          nclrcp(nrecep) = 1
          read(line(31:60),'(3I10)',IOSTAT=ierr) igrdrcp(nrecep),
     &                            irecep(nrecep,1), jrecep(nrecep,1)
          if( ierr .NE. 0 ) then
             write(iout,'(//,A)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Reading cell indexes ',
     &                      'from receptor definition file.  Line read:'
             write(iout,'(A)') line
             write(iout,'(1X,2A,/,A)') 'This version of CAMx requires',
     &           ' the grid number to be specified',' at each receptor.'
             call camxerr()
          endif
          igrd = igrdrcp(nrecep)
          if( igrd .LE. 0 .OR. igrd .GT. ngrid ) then
             write(iout,'(//,A)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Invalid grid number read',
     &                     ' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             write(iout,'(1X,2A,/,A)') 'This version of CAMx requires',
     &           ' the grid number to be specified',' at each receptor.'
             call camxerr()
          endif
          if( jrecep(nrecep,1) .EQ. 0 ) then
             write(iout,'(//,A)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Invalid cell indexes read',
     &                     ' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             write(iout,'(1X,2A,/,A)') 'This version of CAMx requires',
     &           ' the grid number to be specified',' at each receptor.'
             call camxerr()
          endif
          if(irecep(nrecep,1) .LE. 0 .OR. 
     &             irecep(nrecep,1).GT. nxcell(igrd) 
     &                     .OR. jrecep(nrecep,1) .LE. 0 .OR. 
     &                        jrecep(nrecep,1) .GT. nycell(igrd) ) then
             write(iout,'(//,a)') 'WARNING in RERCP:'
             write(iout,'(1X,3A,I5,A,I5,A)') 
     &                  'WARNING: Cell for receptor: ',rcpnam(nrecep),
     &                      'not in specified grid: (',
     &                       irecep(nrecep,1),',',jrecep(nrecep,1),')'
             write(iout,'(1X,A)') '   Skipping receptor ...'
             nrecep = MAX0(0,nrecep-1)
           endif
      endif
c
c   --- if it is a CELL AVERAGE type read the number of cells ---
c
      if( idtmp .EQ. IDAVG ) then
          idrcp(nrecep) = IDAVG
          rcpnam(nrecep) = line(21:30)
          read(line(31:50),'(2I10)',IOSTAT=ierr) igrdrcp(nrecep), 
     &                                              nclrcp(nrecep)
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Reading number of ',
     &               'cells from receptor definition file.  Line read:'
             write(iout,'(1X,2A)') 'This version of CAMx requires the',
     &                  ' grid number to be specified at each receptor.'
             write(iout,'(A)') line
             call camxerr()
          endif
          igrd = igrdrcp(nrecep)
          if( igrd .LE. 0 .OR. igrd .GT. ngrid ) then
             write(iout,'(//,A)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Invalid grid number read',
     &                     ' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             write(iout,'(1X,2A)') 'This version of CAMx requires the',
     &                  ' grid number to be specified at each receptor.'
             call camxerr()
          endif
          if( nclrcp(nrecep) .LE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,A,A,I5)') 'Number of cells in ',
     &                        ' receptor definition file is invalid.'
             write(iout,'(1X,A,A)') 'Line read: ',line
             write(iout,'(1X,2A)') 'This version of CAMx requires the',
     &                  ' grid number to be specified at each receptor.'
             call camxerr()
          endif
          if( nclrcp(nrecep) .GT. MXCELR ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,A,A,I5)') 'Number of cells in ',
     &             ' receptor definition file exceeds maximum: ',MXCELR
             write(iout,'(1X,A,A)') 'Line read: ',line
             write(iout,'(1X,A)') 'Increase the paramater MXCELR.'
             call camxerr()
          endif
          do 10 i=1,nclrcp(nrecep)
              irec = irec + 1
              read(IORRCP,'(A)',IOSTAT=ierr) line
              if( ierr .NE. 0 ) then
                 write(iout,'(//,a)') 'ERROR in RERCP:'
                 write(iout,'(/,1X,A,A,I5)') 'Reading the ',
     &                      'receptor definition file at record: ',irec
                 call camxerr()
              endif
              read(line,'(2I10)',IOSTAT=ierr) irecep(nrecep,i), 
     &                                                jrecep(nrecep,i)
              if( ierr .NE. 0 ) then
                 write(iout,'(//,a)') 'ERROR in RERCP:'
                 write(iout,'(/,1X,A,A,A)') 'Reading cell ',
     &              'indexes from receptor definition file.  Line read:'
                 write(iout,'(A)') line
                 call camxerr()
              endif
              if( irecep(nrecep,i) .LE. 0 .OR. 
     &                 irecep(nrecep,i) .GT. nxcell(igrd) 
     &                     .OR. jrecep(nrecep,i) .LE. 0 .OR. 
     &                         jrecep(nrecep,i) .GT. nycell(igrd) ) then
                 write(iout,'(//,a)') 'WARNING in RERCP:'
                 write(iout,'(1X,3A,I5,A,I5,A)') 
     &                  'WARNING: Cell for receptor: ',rcpnam(nrecep),
     &                      'not in grid: (',
     &                       irecep(nrecep,i),',',jrecep(nrecep,i),')'
                 write(iout,'(1X,A)') '   Skipping receptor ...'
                 nrecep = MAX0(0,nrecep-1)
               endif
   10     continue
      endif
c
c   --- if it is a WALL OF CELLS type read the extent in each direction ---
c
      if( idtmp .EQ. IDWAL ) then
          idrcp(nrecep) = IDWAL
          rcpnam(nrecep) = line(21:30)
          nclrcp(nrecep) = 1
          lwalls = .TRUE.
c
c   --- read the column indexes ---
c
          read(line(31:60),'(3I10)',IOSTAT=ierr) igrdrcp(nrecep), 
     &                                  iwalbg(nrecep), iwalnd(nrecep)
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,A,A,A)') 
     &            'Reading column indexes of ',
     &             CDWAL,' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             write(iout,'(1X,2A)') 'This version of CAMx requires the',
     &                  ' grid number to be specified at each receptor.'
             call camxerr()
          endif
          igrd = igrdrcp(nrecep)
          if( igrd .LE. 0 .OR. igrd .GT. ngrid ) then
             write(iout,'(//,A)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Invalid grid number read',
     &                     ' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             write(iout,'(1X,2A)') 'This version of CAMx requires the',
     &                  ' grid number to be specified at each receptor.'
             call camxerr()
          endif
          if( iwalnd(nrecep) .EQ. 0 ) then
             write(iout,'(//,A)') 'ERROR in RERCP:'
             write(iout,'(/,1X,3A)') 'Invalid cell indexes read',
     &                     ' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             write(iout,'(1X,A)') 'This version of CAMx requires the',
     &                  ' grid number to be specified at each receptor.'
             call camxerr()
          endif
          if( iwalbg(nrecep) .LE. 0 .OR. 
     &                 iwalnd(nrecep) .GT. nxcell(igrd) ) then
            write(iout,'(//,a)') 'WARNING in RERCP:'
            write(iout,'(1X,3A,I5,A,I5,A)') 
     &         'WARNING: Columns for receptor: ',rcpnam(nrecep),
     &         ' not in domain: (',iwalbg(nrecep),',',iwalnd(nrecep),')'
            write(iout,'(1X,A)') '   Skipping receptor ...'
            nrecep = MAX0(0,nrecep-1)
          endif
c
c   --- read the row indexes ---
c
          irec = irec + 1
          read(IORRCP,'(A)',IOSTAT=ierr) line
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,A,A,I5)') 'Reading the ',
     &                      'receptor definition file at record: ',irec
             call camxerr()
          endif
          read(line(41:60),'(2I10)',IOSTAT=ierr) jwalbg(nrecep), 
     &                                                 jwalnd(nrecep)
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,A,A,A)') 'Reading row indexes of ',
     &             CDWAL,' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             call camxerr()
          endif
          if( jwalbg(nrecep) .LE. 0 .OR. 
     &                 jwalnd(nrecep) .GT. nycell(igrd) ) then
            write(iout,'(//,a)') 'WARNING in RERCP:'
            write(iout,'(1X,3A,I5,A,I5,A)') 
     &         'WARNING: Rows for receptor: ',rcpnam(nrecep),
     &         ' not in domain: (',jwalbg(nrecep),',',jwalnd(nrecep),')'
            write(iout,'(1X,A)') '   Skipping receptor ...'
            nrecep = MAX0(0,nrecep-1)
          endif
c
c   --- read the layer indexes ---
c
          irec = irec + 1
          read(IORRCP,'(A)',IOSTAT=ierr) line
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,A,A,I5)') 'Reading the ',
     &                      'receptor definition file at record: ',irec
             call camxerr()
          endif
          read(line(41:60),'(2I10)',IOSTAT=ierr) kwalbg(nrecep), 
     &                                                 kwalnd(nrecep)
          if( ierr .NE. 0 ) then
             write(iout,'(//,a)') 'ERROR in RERCP:'
             write(iout,'(/,1X,A,A,A)') 'Reading layer indexes of ',
     &             CDWAL,' from receptor definition file.  Line read:'
             write(iout,'(A)') line
             call camxerr()
          endif
          if( kwalbg(nrecep) .LE. 0 .OR. 
     &                 kwalnd(nrecep) .GT. nlay(igrd) ) then
            write(iout,'(//,a)') 'WARNING in RERCP:'
            write(iout,'(1X,3A,I5,A,I5,A)') 
     &         'WARNING: Layers for receptor: ',rcpnam(nrecep),
     &         ' not in domain: (',kwalbg(nrecep),',',kwalnd(nrecep),')'
            write(iout,'(1X,A)') '   Skipping receptor ...'
            nrecep = MAX0(0,nrecep-1)
          endif
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
      write(iout,'(//,a)') 'ERROR in RERCP:'
      write(iout,'(/,1X,2A)') 'Cannot open input file: ',
     &                                          rcpfil(:istrln(rcpfil))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
