c***** RDARDDM.F
c
      subroutine rdarddm(igrid,ndate,ttime,nox,noy,nddmspc,emisddm)
c
c-----CAMx v4.02 030709
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions for the DDM process
c   and fills the approproate arrays.  The emissions file for one grid
c   but each emissions groups is read.
c    Argument descriptions:
c     Outputs:
c       emisddm   R    array of emissions for DDM speices
c     Inputs:
c       igrid     I    grid number
c       ndate     I    julian day of current hour
c       ttime     R    current hour
c       nox       I    number of columns in grid
c       noy       I    number of rows in grid
c       nddmspc   I    number of DDM species
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c       10/24/01  Removed BSWAP and converted integer strings to character*4
c       07/19/02  Added seperate source area map for each grids.
c       05/01/03  Time span of emissions must now match emiss update interval
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
c       
c-----------------------------------------------------------------------
c   Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      integer ndate
      real    ttime
      integer nox
      integer noy
      integer nddmspc
      real    emisddm(nox,noy,nddmspc)
c
c-----------------------------------------------------------------------
c   Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  cname
      character*4   iname(10)
      integer       ibgdat, iendat, iounit, igroup, idx, iseg
      integer       i, j, nedge, iptr, ispc, iddm, ibegcl, iendcl
      integer       imap, imod
      real          emstmp, emsgrd(MXCOLA,MXROWA), bgtim, edtim
      logical       lfound, luse, lpass
c
c-----------------------------------------------------------------------
c   Entry point:
c-----------------------------------------------------------------------
c
c
c   --- initialize emissions to zero ---
c
      do i=1,nox
         do j=1,noy
            do l=1,nddmspc
               emisddm(i,j,l) = 0.
            enddo
         enddo
      enddo
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
      else
        nedge = 1
      endif
c
c   --- loop over all of the groups ----
c
      do 10 igroup=1,ngroup
c
c   --- skip if filename not supplied ---
c
        if( .NOT. ltemfl(igrid,igroup) .OR. .NOT. larsrc ) goto 10
c
c   --- set the unit number for surface emissions file ---
c
        iounit = IORTEM + (igrid-1)*ngroup + igroup
        fname = temfil(igrid,igroup)
c
c   --- read the date and time, again ---
c
        lfound = .FALSE.
        lpass = .FALSE.
  111   continue
        read(iounit,END=222) ibgdat, bgtim, iendat, edtim
        ichktm1 = NINT( 1000*(bgtim) )
        ichktm2 = NINT( 1000*(edtim) )
        if( NINT(edtim) .EQ. 0 ) ichktm2 = 24000
        ichkems = NINT( 1000*(dtems/60.) )
        if( (ichktm2 - ichktm1) .NE. ichkems ) then
            write(iout,'(//,a)')'ERROR in READARDDM:'
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
        bgtim = bgtim*100.0
        edtim = edtim*100.0
c
c   --- read the emissions for this hour ---
c
         do 20 ispc=1,nspcem(igrid,igroup)
             read(iounit,ERR=7000) iseg, (iname(i),i=1,10), 
     &                         ((emsgrd(i,j),i=1,nox),j=1,noy)
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
c   --- if the species is a not modelled or not used for DDM ----
c
             idx = idxems(igrid,igroup,ispc)
             if( idx .LE. 0 ) goto 20
c
c  --- find this species in the modeled species list ---
c
             luse = .FALSE.
             do imod=1,nspec
                if( cname .EQ. spname(imod) ) then
                   do iddm=1,nemddm
                      luse = .FALSE.
                      if( emddmsp(iddm) .EQ. cname ) luse = .TRUE.
                      if( emddmsp(iddm) .EQ. NAMALL  ) luse = .TRUE.
                      if( emddmsp(iddm) .EQ. NAMVOC
     &                               .AND. lvocsp(imod) ) luse = .TRUE.
                      if( emddmsp(iddm) .EQ. NAMNOX
     &                               .AND. lnoxsp(imod) ) luse = .TRUE.
c
c  --- if this DDM species matches this modeled species, load it ---
c
                      if( luse ) then
                         do 30 j=2,noy-1
                            ibegcl = 2
                            iendcl = nox - 1
                            if( igrid .EQ. 1 ) then
                               if( ibeg(j) .EQ. -999 ) goto 30
                               ibegcl = ibeg(j)
                               iendcl = iend(j)
                            endif
                            do 40 i=ibegcl,iendcl
c
                              imap = igrmap(igrid,i,j)
                              if( imap .LE. 0 .OR. imap .GT. nregin ) 
     &                                                          goto 40
c
c  --- convert to emissions time and put into array ---
c
                              emstmp = emsgrd(i,j)/(60.*dtems)
                              iptr = iptddm(imod) + nicddm + 
     &                          nbcddm*nedge + (iddm-1)*ngroup*nregin +
     &                                      (imap-1)*ngroup + igroup - 1
                              emisddm(i,j,iptr) = emstmp
c
c  --- next cell ---
c
   40                       continue
   30                    continue
                      endif
c
c  --- next DDM emssions species ---
c
                   enddo
c
c  --- next modeled species ---
c
                endif
             enddo
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
            do i=1,nspcem(igrid,igroup) + 1
                backspace(iounit)
            enddo
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
             if( lpass ) goto 7001
             if( .NOT. lpass ) lpass = .TRUE.
             goto 111
         else
             goto 7001
         endif
c
c  --- get the next file ---
c
  10   continue
       goto 9999
c
c----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in RDARDDM:' 
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'ERROR: Reading emissions ',
     &                  'after hour ',ibgdat, bgtim,' in file: ',fname
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in RDARDDM:' 
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 
     &              'ERROR: Premature end-of-file',
     &              ' in emissions file after hour ',ibgdat, bgtim,
     &                                               ' in file: ',fname
      call camxerr()
c
c----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
