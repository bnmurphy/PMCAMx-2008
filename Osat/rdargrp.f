      subroutine rdargrp(igrid,ndate,ttime,nox,noy,dum1,dum2,
     &                   emsvoc,emsnox)
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c      10/24/01  Removed BSWAP and converted integer strings to character*4
c      05/01/03  Time span of emissions must now match emiss update interval
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'flags.com'
      include 'grid.com'
      include 'bndary.com'
      include 'tracer.com'
c
      character*200 fname
      integer       ndate, ispc
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, igroup, idx, iseg
      integer       i, j
      real          ttime, bgtim, edtim
      real          emsvoc(0:MXTEMF,MXCOLA,MXROWA)
      real          emsnox(0:MXTEMF,MXCOLA,MXROWA)
      real          emsgrd(MXCOLA,MXROWA)
      logical       lfound
c
c   --- loop over all of the groups ----
c
      do 10 igroup=0,ngroup
c
c   --- skip if filename not supplied ---
c
         if( .NOT. ltemfl(igrid,igroup) .OR. .NOT. larsrc ) goto 10
c
c   --- set the unit number for surface emissions file ---
c
         if( igroup .EQ. 0 ) then
             iounit = iarem(1,igrid)
             write(fname,'(A,I3)') 'EMISSIONS -- UNIT ',iarem(1,igrid)
c
c    --- if emissions is regular model emissions file, backup one hour ---
c
             do 90 i=1,nspcem(igrid,igroup)
                 backspace(iounit)
   90        continue
             backspace(iounit)
         else
             iounit = IORTEM + (igrid-1)*ngroup + igroup
             fname = temfil(igrid,igroup)
         endif
c
c   --- read the date and time, again ---
c
        lfound = .FALSE.
  111   continue
        read(iounit,END=222) ibgdat, bgtim, iendat, edtim
        ichktm1 = NINT( 1000*(bgtim) )
        ichktm2 = NINT( 1000*(edtim) )
        if( NINT(edtim) .EQ. 0 ) ichktm2 = 24000
        ichkems = NINT( 1000*(dtems/60.) )
        if( (ichktm2 - ichktm1) .NE. ichkems ) then
            write(iout,'(//,a)')'ERROR in RDARGRP:'
            write(iout,*) 'Time interval in surface emissions file does'
            write(iout,*)  ' not match emissions update time interval.'
            write(iout,*) '   Beginning Date/Time (HHMM): ',ibgdat,bgtim
            write(iout,*) '   Ending Date/Time    (HHMM): ',iendat,edtim
            write(iout,*) '   Emiss Input interval (min): ',dtems
            call camxerr()
        endif
        if(NINT(edtim) .EQ. 0) then
          edtim = 24.
          iendat = iendat - 1
        endif
c
c   --- read the emissions for this hour ---
c
         do 20 ispc=1,nspcem(igrid,igroup)
             read(iounit,ERR=7000) iseg, (iname(i),i=1,10), 
     &                         ((emsgrd(i,j),i=1,nox),j=1,noy)
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
             idx = idxems(igrid,igroup,ispc)
             if( idx .LE. 0 ) goto 20
             if( .NOT. lvocsp(idx) .AND. .NOT. lnoxsp(idx) ) goto 20
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
c   --- convert to PPM (PPMC) and add to totals for VOC or NOx ---
c
             do j=2,noy-1
               do i=2,nox-1
                 emsgrd(i,j) = emsgrd(i,j)/(60.*dtems)
               enddo
             enddo
c
             call fillar(igrid,igroup,MXTEMF,MXCOLA,MXROWA,nox,noy,
     &                   ibeg,iend,lvocsp(idx),lnoxsp(idx),
     &                   crbnum(idx),emsgrd,emsvoc,emsnox)
c
c   --- next species ---
c
  20     continue
c
c   --- if the correct hour has not been found, 
c       go back and read some more else read next file ---
c
         if( .NOT. lfound ) then
            goto 111
         else
            goto 10
         endif
c
c   --- if using 1 day emissions, we need to rewind the file to
c       get the current hour ---
c
  222    continue
         if( le1day ) then
             rewind(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             goto 111
         else
             goto 7001
         endif
c
c  --- get the next file ---
c
  10   continue
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
      return
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDARGRP:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'Reading emissions ',
     &   'after hour ',ibgdat, bgtim,' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDARGRP:'
      write(iout,'(/,1X,2A,I8.5,F4.1,2A)') 'Premature end-of-file',
     &              ' in emissions file after hour ',ibgdat, bgtim,
     &              ' in file: ',fname(:istrln(fname))
      call camxerr()
      end
