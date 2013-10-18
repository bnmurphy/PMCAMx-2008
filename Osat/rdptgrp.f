      subroutine rdptgrp(ndate,ttime,emsvoc,emsnox,izcel)
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.02 030709
c
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     08/18/99  --gwilson--   Added code to implement the override
c                             flag for the source area of point sources
c     10/24/01  Removed BSWAP and converted integer strings to character*4
c     07/05/02  Changed to account for new type of the PiG flag
c     05/01/03  Time span of emissions must now match emiss update interval
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'flags.com'
      include 'grid.com'
      include 'ptemiss.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, igroup, idx, iseg
      integer       npoint, ispc, i, idum, izcel(MXPTSRC)
      integer       ndate
      real          bgtim, edtim, rdum
      real          ttime
      real          emsvoc(0:MXTEMF,MXPTSRC), emsnox(0:MXTEMF,MXPTSRC)
      real          emspts(MXPTSRC)
      logical       lfound, lpig
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c
      lpig = .FALSE.
      if( ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG ) lpig = .TRUE.
c
c   --- loop over all of the groups ----
c
      do 10 igroup=0,ngroup
c
c   --- skip if filename not supplied ---
c
         if( .NOT. ltptfl(igroup) ) goto 10
c
c   --- set the unit number for file ---
c
         if( igroup .EQ. 0 ) then
             iounit = iptem
             write(fname,'(A,I3)') 'PTSOURCE -- UNIT ',iptem
             do 20 i=1,nspcpt(igroup)
                 backspace(iounit)
   20        continue
             backspace(iounit)
             backspace(iounit)
             backspace(iounit)
         else
             iounit = IORTPT + igroup
             fname = tptfil(igroup)
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
            write(iout,'(//,a)')'ERROR in RDPTGRP:'
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
c   --- read the number of points and point locations ---
c
         read(iounit,ERR=7000) iseg, npoint
         if( npoint .NE. nptsrc ) goto 7002
         if( npoint .LE. 0 ) goto 9999
         read(iounit,ERR=7000,END=7001)  (idum, idum,
     &                                izcel(i), rdum, rdum,i=1,npoint)
c
c   --- read the emissions for this hour ---
c
         do 40 ispc=1,nspcpt(igroup)
            read(iounit,ERR=7000,END=7001) iseg, (iname(i),i=1,10),
     &                                            (emspts(i),i=1,nptsrc)
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c

             idx = idxpts(igroup,ispc)
             if( idx .LE. 0 ) goto 40
             if( .NOT. lvocsp(idx) .AND. .NOT. lnoxsp(idx) ) goto 40
c
c   --- if date and time does not match this hour, skip this record ---
c
             if( le1day ) then
                 if( bgtim .NE. ttime ) goto 40
             else
                 if( ndate .NE. ibgdat .OR. bgtim .NE. ttime ) goto 40
             endif
             lfound  = .TRUE.
c
c   --- convert to PPM (PPMC)
c
             do i=1,nptsrc
               emspts(i) = emspts(i)/(60.*dtems)
             enddo
             call fillpt(igroup,nptsrc,MXTEMF,MXPTSRC,lpig,lpigsa,
     &                   lvocsp(idx),lnoxsp(idx),crbnum(idx),emsvoc,
     &                   emsnox,emspts)

c
c   --- next species ---
c
  40     continue
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
             read(iounit)
             read(iounit)
             goto 111
         else
             goto 7001
         endif
c
c  --- get the next file ---
c
  10  continue
c
 9999 continue
      return
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDPTGRP:'
      write(iout,'(/,1X,A,I8.5,F8.1,2A)') 
     &      'Reading emissions after hour ',ibgdat, bgtim,' in file: ',
     &                                            fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDPTGRP:'
      write(iout,'(/,1X,3A)') 'Premature end-of-file reading ',
     &                  'emissions from file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDPTGRP:'
      write(iout,'(/,1X,A,I10,2A)') 'Number of points: ',npoint,
     &     ' is not consistent with regular emissions in file: ',
     &                                            fname(:istrln(fname))
      call camxerr()
c
      end
