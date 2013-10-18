c***** RDPTRT.F
c
      subroutine rdptrt(ndate,ttime)
c
c-----CAMx v4.02 030709
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions for the RTRAC process
c   and fills the approproate arrays.  The emissions file for one grid
c   but each emissions groups is read.
c    Argument descriptions:
c     Outputs:
c     Inputs:
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/24/01   --gwilson--   Originial development
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'flags.com'
      include 'chmstry.com'
      include 'grid.com'
      include 'bndary.com'
      include 'tracer.com'
      include 'rtracchm.com'
c       
c-----------------------------------------------------------------------
c   Argument declarations:
c-----------------------------------------------------------------------
c
      integer ndate
      real    ttime
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  cname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, idx, iseg
      integer       i, j, ispc
      integer       icel, jcel, npoint, idum, imod
      real          emspnt(MXPTSRC), bgtim, edtim, xloctmp, yloctmp
      real          emstmp
      logical       lfound, luse, lpass
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c  --- initialize emissions to zero ---
c
      do i=1,MXPTSRC
         do j=1,MXTRSP
            sapnts(i,j) = 0.
         enddo
      enddo
c
c   --- skip if filename not supplied ---
c
      if( .NOT. ltptfl(1) .OR. .NOT. lptsrc ) goto 9999
c
c   --- set the unit number for surface emissions file ---
c
      iounit = IORTPT
      fname = tptfil(1)
c
c   --- read the date and time, again ---
c
      lfound = .FALSE.
      lpass = .FALSE.
  111 continue
      read(iounit,END=222) ibgdat, bgtim, iendat, edtim
      bgtim = bgtim*100.0
      edtim = edtim*100.0
c
c   --- read the emissions for this hour ---
c
      read(iounit,ERR=7000,END=7001) iseg, npoint
      if( npoint .GT. MXPTSRC ) goto 7002
      read(iounit,ERR=7000,END=7001) idum
      do 20 ispc=1,nspcpt(1)
          read(iounit,ERR=7000) iseg, (iname(i),i=1,10), 
     &                                        (emspnt(i),i=1,npoint)
          write(cname,'(10A1)') iname
c
c   --- if date and time does not match this hour, skip this record ---
c
           if( le1day ) then
               if( bgtim .NE. ttime ) goto 20
           else
               if( ndate .NE. ibgdat .OR. bgtim .NE. ttime ) goto 20
           endif
           lfound  = .TRUE.
c
c   --- if the species is a not modelled or not used for RTRAC ----
c
           idx = idxpts(1,ispc)
           if( idx .LE. 0 ) goto 20
c
c  --- find this species in the modeled species list ---
c
           luse = .FALSE.
           do imod=1,nspec
              if( cname .EQ. ptname(imod) ) then
                 luse = .TRUE.
                 do 30 i=1,npoint
                    if( llatlon ) then
                       xloctmp = xlocpt(i) - xorg
                       yloctmp = ylocpt(i) - yorg
                    else
                       xloctmp = xlocpt(i)/1000. - xorg
                       yloctmp = ylocpt(i)/1000. - yorg
                    endif
                    icel = 1 + INT( xloctmp/delx )
                    jcel = 1 + INT( yloctmp/dely )
                    if( icel .LE. 0 .OR. icel .GT. ncol(1) ) goto 30
                    if( jcel .LE. 0 .OR. jcel .GT. nrow(1) ) goto 30
                    if( ibeg(jcel) .EQ. -999 ) goto 30
                    if( icel .LT. ibeg(jcel) .OR. 
     &                                   icel .GT. iend(jcel) ) goto 30
c
c  --- convert to emissions time and put into array ---
c
                    emstmp = emspnt(i)/(60.*dtems)
                    if( emstmp .GT. rtlbnd(imod) )
     &                                      sapnts(i,imod) = emstmp
c
c  --- next point ---
c
   30            continue
              endif
c
c  --- next RTRAC emssions species ---
c
           enddo
c
c  --- if not found in model species list, write a message ----     
c   
           if( .NOT. luse ) then
                 write(idiag,'(1X,4A)') 'Species in RTRAC ',
     &              'point source file: ',cname(:istrln(cname)),
     &                      ' not found in species list ... Skipping.'
           endif
c
c   --- next species ---
c
  20  continue
c
c   --- if the correct hour has not been found, 
c       go back and read some more else read next file ---
c
      if( .NOT. lfound ) then
         goto 111
      else
         goto 9999
      endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222 continue
      if( le1day ) then
          rewind(iounit)
          read(iounit)
          read(iounit)
          read(iounit)
          read(iounit)
          if( lpass ) goto 7001
          if( .NOT. lpass ) lpass = .TRUE.
          goto 111
      else
          goto 7001
      endif
c
c  --- get the next file ---
c
      goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in RDPTRT:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'ERROR: Reading emissions ',
     &    'after hour ',ibgdat, bgtim,' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in RDPTRT:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)')
     &              'ERROR: Premature end-of-file',
     &              ' in point source file after hour ',ibgdat, bgtim,
     &                                ' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in RDPTRT:'
      write(iout,'(/,1X,A,I10,2A)') 'ERROR:  Number of points: ',npoint,
     &                   ' exceeds max in file: ',fname(:istrln(fname))
      write(iout,'(1X,A)') 'Increase the parameter MXPTSRC.'
      call camxerr()
c
c----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
