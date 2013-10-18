      subroutine rdfgcon(idate,btim)
c
c-----CAMx v4.02 030709
c
c     RDFGCON reads the instantaneous file for the fine grids
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        idate               beginning date of simulation (YYJJJ)
c        btim                beginning time of simulation (HHMM)
c
c     Output arguments:
c        none
c
c     Routines called:
c        none
c
c     Called by:
c        READCNC
c
      include 'camx.prm'
      include 'filunit.com'
      include 'camxfld.com'
      include 'grid.com'
      include 'chmstry.com'
c
      integer   idate
      real      btim
c
      character*60 msg
      character*10 cname, cspec(MXSPEC)
      integer      numfin, nspcs, ndate, indate
      integer      ixfb, jyfb, ixfe, jyfe, nhf, nvf
      integer      nfx(MXGRID), nfy(MXGRID), nfz(MXGRID)
      integer      ifgptr, ifglvl, idx, i, j, k, l, igrd
      integer      ibgdat, irdgrd(MXGRID)
      real         ttime, timein, begtim, cnctmp(MXCOLA,MXROWA)
      logical      lmatch, lread(MXGRID)
c
c-----Entry point:
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
      if(  irstf .LE. 0 ) goto 222
c
c --- otherwise get whatever grids are available ---
c
      ndate = idate
      ttime = btim
      if (ttime.eq.2400.) then
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
c-----File description header
c
      read(irstf,ERR=7000) msg
      read(irstf,ERR=7000) numfin, nspcs
c
c-----Species list
c
      read(irstf,ERR=7000) (cspec(i),i=1,nspcs)
c
c-----Make sure species list is compatable
c
      do 10 l = 1,nspec
        do 20 lic = 1,nspcs
          cname = cspec(lic)
          if (cname.eq.'HNO2      ') cname = 'HONO      '
          if (cname.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        cname = 'FORM      '
          if (cname.eq.spname(l)) then
            if (lic.ne.licmap(l,1)) goto 7004
            goto 10
          endif
   20   continue
   10 continue
c
c-----Read the fine grid descriptions, and check for consistancy
c
      do 30 ifin=1,numfin
        read(irstf,ERR=7000) ixfb, jyfb, ixfe, jyfe, nhf, nvf, 
     &                           nfx(ifin), nfy(ifin), nfz(ifin), 
     &                           ifgptr, ifglvl
        do igrd = 2,ngrid
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
   30 continue
c
c-----Make sure time span is correct
c
  111 continue
      read(irstf,ERR=7000,END=7002) timein, indate
      ibgdat = indate
      begtim = timein
      if (begtim.eq.2400.) then
        begtim = 0.
        ibgdat = ibgdat + 1
        if( MOD(ibgdat,1000) .GT. 365 ) then
           if( MOD(INT(ibgdat/1000),4) .EQ. 0 ) then
              if( MOD(ibgdat,1000) .EQ. 367 )
     &                       ibgdat = (INT(ibgdat/1000)+1)*1000 + 1
           else
              ibgdat = (INT(ibgdat/1000)+1)*1000 + 1
           endif
        endif
      endif
c
c-----Skip if not the correct hour
c
      if (ibgdat.ne.ndate .or. abs(begtim-ttime).gt.0.00001) then
        do 40 ifin=1,numfin
          do 50 l=1,nspcs
            do 60 k=1,nfz(ifin)
   60       continue
   50     continue
   40   continue
        goto 111
      endif
c
c-----Read the data for this hour
c
      do 70 ifin=1,numfin
        igrd = 0
        do n = 2,ngrid
          if( irdgrd(n) .EQ. ifin ) igrd = n
        enddo
        do 80 l=1,nspcs
          do 90 k=1,nfz(ifin)
            read(irstf) ((cnctmp(i,j),i=1,nfx(ifin)),
     &                                          j=1,nfy(ifin))
c
c-----Find the grid in this setup and load into gridded array---
c
            if( igrd .NE. 0 ) then
               do 11 j=1,nrow(igrd)
                 do 21 i=1,ncol(igrd)
                   idx = i + ncol(igrd)*(j-1) + 
     &                   ncol(igrd)*nrow(igrd)*(k-1) +
     &                   ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                   conc(iptr4d(igrd)-1+idx) = cnctmp(i,j)
   21            continue
   11          continue
            endif
   90     continue
   80   continue
        if( igrd .NE. 0 ) write(iout,'(a40,f7.0,i8.5,a,i3)')
     &     'Read fine grid restart file at ',timein, indate,' grid',igrd
   70 continue
c
c-----Load data for the grids not found in restart file----
c
  222 continue
      do ip = 1,ngrid
         do ic = 1,nchdrn(ip)
           ig = idchdrn(ic,ip)
           if( .NOT. lread(ig) ) then
              write(iout,'(a40,f7.0,i8.5,a,i3)')
     &        'Assigning concs from parent grid at ',
     &         btim, idate,' grid',ig
              call rassgn4d(ncol(ip),nrow(ip),nlay(ip),nspec,i1(ig),
     &                      j1(ig),nmesh(ig),nmshv(1,ig),ncol(ig),
     &                      nrow(ig),nlay(ig),conc(iptr4d(ip)),
     &                      conc(iptr4d(ig)) )
           endif
         enddo
      enddo
c
      goto 9999
c
c-----Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDFGCON:'
      write(iout,9000,ERR=9999)'Reading the header of the ',
     &     'concentration initialization file for fine grids. '
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDFGCON:'
      write(iout,9000,ERR=9999) 'Concentration ',
     &   'initialization file for fine grids is not for correct ',
     &                               'time period; '
      write(iout,9001,ERR=9999) 'Premature EOF reached.'
      write(iout,9001,ERR=9999) 'Simulation Start: ',idate,btim
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in RDFGCON:'
      write(iout,9000,ERR=9999) 'Species list for ',
     &            'concentration initialization file for fine grids ',
     &          'is inconsistent with user parameters.'
      write(iout,9004,ERR=9999) 'Species number: ',l,' should be ',
     &                      spname(l),' and ',cname,' was read instead.'
      call camxerr()
c
c-----Format statements
c
 9000 format(/,1X,6A)
 9001 format(1X,A,I10.5,F10.1)
 9004 format(1X,A,I5,5A)
c
 9999 continue
      return
      end
