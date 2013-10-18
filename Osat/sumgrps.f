      subroutine sumgrps(ndlast,ttlast,emstot,emslft,emsbas,emsoth,
     &                   emssum)
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
      include 'bndary.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, igroup, idx
      integer       iseg, npoint, idum, ispc, i, j, ndlast
      integer       igrid
      real          bgtim, edtim, emssum(MXSPEC,MXTRSP)
      real          emsgrd(MXCOLA,MXROWA), emspnt(MXPTSRC)
      real          emsbas(MXSPEC,MXTRSP), emsoth(MXSPEC,MXTRSP)
      real          emslft(MXCOLA,MXROWA), emstot(MXCOLA,MXROWA)
      real          ttlast
      logical       lpass
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- loop over all of the groups ----
c
      do 20 igroup=0,ngroup
c
c   --- only process if filename is supplied for group ---
c
         do 25 igrid=1,ngrid
            if( .NOT. ltemfl(igrid,igroup) .OR. .NOT. larsrc) goto 25
c
c   --- set the unit number for surface emissions file ---
c
            lpass = .FALSE.
            if( igroup .EQ. 0 ) then
                iounit = iarem(1,igrid)
                write(fname,'(A,I3)') 'EMISSIONS -- UNIT ',
     &                                                iarem(1,igrid)
            else
                iounit = IORTEM + (igrid-1)*ngroup + igroup
                fname = temfil(igrid,igroup)
            endif
c
c   --- read the date and time, again ---
c
  111     continue
          read(iounit,END=333) ibgdat, bgtim, iendat, edtim
          ichktm1 = NINT( 1000*(bgtim) )
          ichktm2 = NINT( 1000*(edtim) )
          if( NINT(edtim) .EQ. 0 ) ichktm2 = 24000
          ichkems = NINT( 1000*(dtems/60.) )
          if( (ichktm2 - ichktm1) .NE. ichkems ) then
            write(iout,'(//,a)')'ERROR in SUMGRPS:'
            write(iout,*) 'Time interval in surface emissions file does'
            write(iout,*) ' not match emissions update time interval.'
            write(iout,*) '  Beginning Date/Time (HHMM): ',ibgdat,bgtim
            write(iout,*) '  Ending Date/Time    (HHMM): ',iendat,edtim
            write(iout,*) '  Emiss Input interval (min): ',dtems
            call camxerr()
          endif
          if( INT(edtim) .EQ. 0 ) then
            edtim = 24.0
            iendat = iendat - 1
          endif
          if( le1day ) iendat = ndlast
c
c   --- read the emissions for this hour ---
c
            do 30 ispc=1,nspcem(igrid,igroup)
                read(iounit) iseg, (iname(i),i=1,10), 
     &                 ((emsgrd(i,j),i=1,ncol(igrid)),j=1,nrow(igrid))
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
                idx = idxems(igrid,igroup,ispc)
                if( idx .LE. 0 ) goto 30
                if( .NOT. lvocsp(idx) .AND. .NOT. lnoxsp(idx)) goto 30
c
                call sum1grd(igroup,igrid,idx,ibeg,iend,emssum,
     &                          emsgrd,emsbas,emsoth,emslft,emstot)
c
c   --- next species ---
c
  30        continue
c
c   --- check if end of simulation read ---
c
            if( iendat .LT. ndlast  .OR. (iendat .EQ. ndlast .AND. 
     &                          INT(edtim) .LT. INT(ttlast)) ) goto 111
            goto 25
c
c   --- if using 1 day emissions, we might need to keep going
c       to finish out the simulation ----
c
  333       continue
            if( le1day ) then
                if( lpass ) goto 7005
                lpass = .TRUE.
                rewind(iounit)
                read(iounit)
                read(iounit)
                read(iounit)
                read(iounit)
                goto 111
            endif
   25    continue 
c
c   --- only process if filename is supplied for group ---
c
         if( .NOT. ltptfl(igroup) .OR. .NOT. lptsrc ) goto 20
c
c   --- set the unit number for elevated points emissions file ---
c
         lpass = .FALSE.
         if( igroup .EQ. 0 ) then
             iounit = iptem
             write(fname,'(A,I3)') 'PTSOURCE -- UNIT ',iptem
         else
             iounit = IORTPT + igroup
             fname = tptfil(igroup)
         endif
c
c   --- read the date and time, again ---
c
  222    continue
         read(iounit,END=444) ibgdat, bgtim, iendat, edtim
         ichktm1 = NINT( 1000*(bgtim) )
         ichktm2 = NINT( 1000*(edtim) )
         if( NINT(edtim) .EQ. 0 ) ichktm2 = 24000
         ichkems = NINT( 1000*(dtems/60.) )
         if( (ichktm2 - ichktm1) .NE. ichkems ) then
            write(iout,'(//,a)')'ERROR in SUMGRPS:'
            write(iout,*) 'Time interval in surface emissions file does'
            write(iout,*) ' not match emissions update time interval.'
            write(iout,*) '  Beginning Date/Time (HHMM): ',ibgdat,bgtim
            write(iout,*) '  Ending Date/Time    (HHMM): ',iendat,edtim
            write(iout,*) '  Emiss Input interval (min): ',dtems
            call camxerr()
         endif
         if( INT(edtim) .EQ. 0 ) then
             edtim = 24.0
             iendat = iendat - 1
         endif
         if( le1day ) iendat = ndlast
c
c   --- read the emissions for this hour ---
c
         read(iounit,ERR=7000,END=7001) iseg, npoint
         if( npoint .GT. MXPTSRC ) goto 7002
         read(iounit,ERR=7000,END=7001) idum
         do 60 ispc=1,nspcpt(igroup)
            read(iounit) iseg, (iname(i),i=1,10), 
     &                                        (emspnt(i),i=1,npoint)
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
            idx = idxpts(igroup,ispc)
            if( idx .LE. 0 ) goto 60
            if( .NOT. lvocsp(idx) .AND. .NOT. lnoxsp(idx)) goto 60
c
c   --- sum up the emissions for each point ---
c
            call sum1pnt(igroup,idx,npoint,emsbas,emsoth,emslft,
     &                                          emstot,emspnt,emssum)
c
c   --- next species ---
c
  60     continue
         if( iendat .LT. ndlast .OR. (iendat .EQ. ndlast .AND. 
     &                          INT(edtim) .LT. INT(ttlast)) ) goto 222
c
         goto 20
c
c   --- if using 1 day emissions, we might need to keep going
c       to finish out the simulation ----
c
 444     continue
         if( le1day ) then
             if( lpass ) goto 7005
             lpass = .TRUE.
             rewind(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             read(iounit)
             goto 222
         endif
c
c  --- get the next file ---
c
  20  continue
c
      return
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,A,I10.5,F10.1,2A)') 
     &      'Reading emissions after hour ',ibgdat, bgtim,
     &      ' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,3A)') 'Premature end-of-file reading ',
     &                    'emissions from file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,A,I10,2A)') 'Number of points: ',npoint,
     &                   ' exceeds max in file: ',fname(:istrln(fname))
      write(iout,'(1X,A)') 'Increase the parameter MXPTSRC.'
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in SUMGRPS:'
      write(iout,'(/,1X,2A,/,10X,2A)') 'Emissions file does cover',
     &    ' entire simulation period.','File = ',fname(:istrln(fname))
      call camxerr()
c
      end
