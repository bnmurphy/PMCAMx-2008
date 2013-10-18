c***** RDPTDDM.F
c
      subroutine rdptddm(ndate,ttime)
c
c-----CAMx v4.02 030709
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions for the DDM process
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
c     10/24/01  Removed BSWAP and converted integer strings to character*4
c     07/19/02  Added seperate source area map for each grids.
c     05/01/03  Time span of emissions must now match emiss update interval
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
      integer       ibgdat, iendat, iounit, igroup, idx, iseg
      integer       i, j, nedge, iptr, ispc, iddm, igrd
      integer       icel, jcel, imap, npoint, idum, imod
      real          emspnt(MXPTSRC), bgtim, edtim, emstmp
      real          xloccrs, yloccrs, xloctmp, yloctmp
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
        if( .NOT. ltptfl(igroup) .OR. .NOT. lptsrc ) goto 10
c
c   --- set the unit number for surface emissions file ---
c
        iounit = IORTPT + igroup
        fname = tptfil(igroup)
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
            write(iout,'(//,a)')'ERROR in READPTDDM:'
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
        read(iounit,ERR=7000,END=7001) iseg, npoint
        if( npoint .GT. MXPTSRC ) goto 7002
        read(iounit,ERR=7000,END=7001) idum
        do 20 ispc=1,nspcpt(igroup)
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
c   --- if the species is a not modelled or not used for DDM ----
c
             idx = idxpts(igroup,ispc)
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
                         do 30 i=1,npoint
                            if( llatlon ) then
                               xloccrs = xlocpt(i) - xorg
                               yloccrs = ylocpt(i) - yorg
                            else
                               xloccrs = xlocpt(i)/1000. - xorg
                               yloccrs = ylocpt(i)/1000. - yorg
                            endif
                            icel = 1 + INT( xloccrs/delx )
                            jcel = 1 + INT( yloccrs/dely )
                            if( icel .LE. 0 .OR. icel .GT. ncol(1) ) 
     &                                                          goto 30
                            if( jcel .LE. 0 .OR. jcel .GT. nrow(1) ) 
     &                                                          goto 30
                            if( ibeg(jcel) .EQ. -999 ) goto 30
                            if( icel .LT. ibeg(jcel) .OR. 
     &                                   icel .GT. iend(jcel) ) goto 30
c
c   --- find out if a nest contains this source  ---
c
                            igrd = 1
                            do ig = 2,ngrid
                               xloctmp = xloccrs - (inst1(ig)-1)*delx
                               yloctmp = yloccrs - (jnst1(ig)-1)*dely
                               if( xloctmp .GT. 0. .AND. 
     &                                         yloctmp .GT. 0. ) then
                                    ii = 2 + INT( xloctmp/delx *
     &                                            FLOAT( meshold(ig) ) )
                                    jj = 2 + INT( yloctmp/dely *
     &                                            FLOAT( meshold(ig) ) )
                                    if( ii .GT. 1 .AND. jj .GT. 1 .AND. 
     &                                   ii .LT. ncol(ig) .AND. 
     &                                           jj .LT. nrow(ig) ) then
                                       igrd = ig
                                       icel = ii
                                       jcel = jj
                                    endif
                               endif
                            enddo
c
c  --- get the region for this cell from mapping array ----
c
                            imap = igrmap(igrd,icel,jcel)
                            if( imap .LE. 0 .OR. imap .GT. nregin ) 
     &                                                         goto 30
c
c  --- convert to emissions time and put into array ---
c
                            emstmp = emspnt(i)/(60.*dtems)
                            iptr = iptddm(imod) + nicddm + 
     &                          nbcddm*nedge + (iddm-1)*ngroup*nregin +
     &                                      (imap-1)*ngroup + igroup - 1
                              if( emstmp .GT. bdnl(imod) )
     &                                        sapnts(i,iptr) = emstmp
c
c  --- next point ---
c
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
            do i=1,nspcpt(igroup) + 3
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
      write(iout,'(//,a)')'ERROR in READPTDDM:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)') 'ERROR: Reading emissions ',
     &    'after hour ',ibgdat, bgtim,' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)')'ERROR in READPTDDM:'
      write(iout,'(/,1X,2A,I8.5,F8.1,2A)')
     &              'ERROR: Premature end-of-file',
     &              ' in point source file after hour ',ibgdat, bgtim,
     &                                ' in file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)')'ERROR in READPTDDM:'
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
