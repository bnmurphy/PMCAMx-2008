c*** RDFGSA 
c
      subroutine rdfgsa(idate,btim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine reads the instantaneous file for the tracer
c     species. 
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c    Argument description:
c     Outputs:
c     Inputs:
c        idate   I   begining date of simulation (YYJJJ)
c        btim    R   begining time of simulation
c
c----------------------------------------------------------------------- 
c   LOG: 
c----------------------------------------------------------------------- 
c 
c     11/06/01   --cemery--     Input dates are now Julian
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*60 msg
      character*10 cname, cspec(MXTRSP)
      integer      numfin, nspcs, ndate, indate
      integer      ixfb, jyfb, ixfe, jyfe, nhf, nvf
      integer      nfx(MXGRID), nfy(MXGRID), nfz(MXGRID)
      integer      ifgptr, ifglvl, idx, i, j, k, l, igrd, ifin
      integer      ibgdat, irdgrd(MXGRID)
      real         ttime, timein, begtim, cnctmp(MXCOLA,MXROWA)
      logical      lmatch, lread(MXGRID)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- initialize all of the tracer concs to zero ----
c
      do igrd=2,ngrid
         do l=1,nsaspc
            do k=1,nlay(igrd)
               do j=1,nrow(igrd)
                  do i=1,ncol(igrd)
                     idx = i + ncol(igrd)*(j-1) + 
     &                      ncol(igrd)*nrow(igrd)*(k-1) +
     &                           ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                     ptconc(ipsa3d(igrd)-1+idx) = 0.
                  enddo
               enddo
            enddo
         enddo
      enddo
c
c --- initialize the "read-in" flag ---
c
      do i=1,ngrid
        lread(i) = .FALSE.
      enddo
c
c --- if the file was not supplied then all of the grids are
c     to be assigned ----
c
      if( inifil(IDXFIN) .EQ. ' ' ) goto 222
c
c   --- open the fine grids instantaneous file ---
c
      open(unit=IORINI,file=inifil(IDXFIN),
     &                  form='UNFORMATTED',ERR=7000,status='UNKNOWN')
c
c   --- convert the date and time ---
c
      ndate = idate
      ttime = btim
      if( ttime .EQ. 2400.0 ) then
          ttime = 0.
          ndate = ndate + 1
          if( MOD(ndate,1000) .GT. 365 ) then
             if( MOD(INT(ndate/1000),4) .EQ. 0 ) then
                if( MOD(ndate,1000) .EQ. 367 )
     &                     ndate = (INT(ndate/1000)+1)*1000 + 1
             else
                ndate = (INT(ndate/1000)+1)*1000 + 1
             endif
          endif
      endif
c
c   --- file description header ---
c
      read(IORINI,ERR=7000) msg
      read(IORINI,ERR=7000) numfin, nspcs
c
c  ---- species list ---
c
      read(IORINI,ERR=7000) (cspec(i),i=1,nspcs)
c
c  --- make sure species list is compatable, up to timing tracers ---
c
      do 10 l=1,nspcs
         cname = cspec(l)
         if( cname .NE. ptname(l) ) goto 7004
   10 continue
c
c  --- read the fine grid descriptions, and check for 
c      consistancy ---
c
      do 20 ifin=1,numfin
          read(IORINI,ERR=7000) ixfb, jyfb, ixfe, jyfe, nhf, nvf, 
     &               nfx(ifin), nfy(ifin), nfz(ifin), ifgptr, ifglvl
          do igrd=2,ngrid
             lmatch = .TRUE.
             if( ixfb .NE. inst1(igrd) ) lmatch = .FALSE.
             if( jyfb .NE. jnst1(igrd) ) lmatch = .FALSE.
             if( ixfe .NE. inst2(igrd) ) lmatch = .FALSE.
             if( jyfe .NE. jnst2(igrd) ) lmatch = .FALSE.
             if( nhf .NE. meshold(igrd) ) lmatch = .FALSE.
             if( nvf .NE. meshold(igrd) ) lmatch = .FALSE.
             if( nfx(ifin) .NE. ncol(igrd) ) lmatch = .FALSE.
             if( nfy(ifin) .NE. nrow(igrd) ) lmatch = .FALSE.
             if( nfz(ifin) .NE. nlay(igrd) ) lmatch = .FALSE.
             if( lmatch ) then 
               lread(mapgrd(igrd)) = .TRUE.
               irdgrd(mapgrd(igrd)) = ifin
             endif
          enddo
   20  continue
c
c  --- make sure time span is correct ---
c
      read(IORINI,ERR=7000) timein, indate
      ibgdat = indate
      begtim = timein
      if( begtim .EQ. 2400.0 ) then
          begtim = 0.
          ibgdat = ibgdat + 1
          if( MOD(ibgdat,1000) .GT. 365 ) then
             if( MOD(INT(ibgdat/1000),4) .EQ. 0 ) then
                if( MOD(ibgdat,1000) .EQ. 367 )
     &                     ibgdat = (INT(ibgdat/1000)+1)*1000 + 1
             else
                ibgdat = (INT(ibgdat/1000)+1)*1000 + 1
             endif
          endif
      endif
      if( ibgdat .NE. ndate .OR. 
     &                        ABS(begtim-ttime) .GT. 0.00001 ) goto 7002 
c
c   --- read the data for this hour ----
c
      do 50 ifin=1,numfin
         do 60 l=1,nspcs
            do 70 k=1,nfz(ifin)
               read(IORINI) ((cnctmp(i,j),i=1,nfx(ifin)),
     &                                          j=1,nfy(ifin))
c
c-----Find the grid in this setup ----
c
               do igrd = 2,ngrid
                  if( irdgrd(igrd) .EQ. ifin ) then
                     if( l .EQ. 1 .and. k .EQ. 1 )
     &                     write(iout,'(a40,f5.0,i8.5,a,i3)')
     &                        'Read fine grid restart file at ',
     &                                    timein, ibgdat,' grid',igrd
                     do 80 j=1,nrow(igrd)
                        do 90 i=1,ncol(igrd)
                           idx =  i + ncol(igrd)*(j-1) + 
     &                         ncol(igrd)*nrow(igrd)*(k-1) +
     &                           ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                           ptconc(ipsa3d(igrd)-1+idx) = cnctmp(i,j)
   90                   continue
   80                continue
                  endif
               enddo
   70       continue
   60    continue
   50 continue
c
c  --- close file and return to calling routine ---
c
      close(IORINI)
c
c-----Load data for the grids not found in restart file----
c
  222 continue
      do ip = 1,ngrid
         do ic = 1,nchdrn(ip)
           ig = idchdrn(ic,ip)
           if( .NOT. lread(ig) ) then
              write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Assigning SA concs from parent grid at ',
     &                             btim, idate,' grid',ig
              call rassgn4d(ncol(ip),nrow(ip),nlay(ip),nsaspc,i1(ig),
     &                      j1(ig),nmesh(ig),nmshv(1,ig),ncol(ig),
     &                      nrow(ig),nlay(ig),ptconc(ipsa3d(ip)),
     &                      ptconc(ipsa3d(ig)) )
           endif
         enddo
      enddo
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDFGSA:'
      write(iout,9000,ERR=9999)'Reading the header of the ',
     &   'tracer concentration initialization file for fine grids: ',
     &                          inifil(IDXFIN)(:istrln(inifil(IDXFIN)))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDFGSA:'
      write(iout,9000,ERR=9999) 'Tracer concentration ',
     &   'initialization file for fine grids is not for correct ',
     &        'time period; ',inifil(IDXFIN)(:istrln(inifil(IDXFIN)))
      write(iout,9001,ERR=9999) 'Simulation Start: ',idate,btim
      write(iout,9001,ERR=9999) 'Initialization File: ',ibgdat,timein
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in RDFGSA:'
      write(iout,9000,ERR=9999) 'Species list for tracer ',
     &            'concentration initialization file for fine grids ',
     &         'is inconsistent with user parameters ',
     &                       inifil(IDXFIN)(:istrln(inifil(IDXFIN)))
      write(iout,9004,ERR=9999) 'Species number: ',l,' should be ',
     &                      ptname(l),' and ',cname,' was read instead.'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,1X,6A)
 9001 format(1X,A,I10.5,F10.1)
 9004 format(1X,A,I5,5A)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
