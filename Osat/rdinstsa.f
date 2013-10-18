c*** RDINSTSA
c
      subroutine rdinstsa(idum,idate,btim,nox,noy,noz,nspsa,saconc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine reads the instantaneous file for the tracer
c     species.  A file is read for each grid nest and the data is
c     stored in the tracer concentration array.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c    Argument description:
c     Outputs:
c        saconc  R   intitial concentrations
c     Inputs: 
c        idate   I   begining date of simulation (YYJJJ)
c        btim    R   begining time of simulation
c        nox     I   number of X cells in grid
c        noy     I   number of Y cells in grid
c        noz     I   number of layers in grid
c        nspsa   I   number of species in conc array
c
c-----------------------------------------------------------------------
c   LOG:
c     1/20/99   Grid cell size from file should be meters for all
c               cartesian projections (UTM, LCP, PSP)
c     10/24/01  Removed BSWAP and converted integer strings to character*4
c     11/06/01  Input dates are now Julian
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'flags.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
      integer   nox
      integer   noy
      integer   noz
      integer   nspsa
      real      saconc(nox,noy,noz,nspsa)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 cname
      character*4  fname(10), fnote(60), ispec(10,MXTRSP)
      integer      nsegs, nspcs, ibgdat, iendat
      integer      nzonin, noxgin, noygin, nzclin, nvloin, nvupin, idum
      integer      ntoday, isegmt, ndate
      integer      i, j, k, l
      real         begtim, endtim, xutmin, yutmin, xorgin, yorgin
      real         delxin, delyin, dzsfin, dzmlin, dzmuin, ttime
      real         cnctmp(MXCOLA,MXROWA)
      logical      lmatch
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ndate = idate
      ttime = btim / 100. 
c
c   --- open the coarse grid instantaneous file ---
c
      open(unit=IORINI,file=inifil(IDXCRS),
     &                    form='UNFORMATTED',ERR=7006,status='UNKNOWN')
c
c   --- file description header ---
c
      read(IORINI,ERR=7000,END=7000) fname, fnote, nsegs, nspcs, ibgdat, 
     &                                        begtim, iendat, endtim
c
c  --- make sure this is the right file ---
c
      write(cname,'(10A1)') (fname(i),i=1,10) 
      if( cname .NE. 'AIRQUALITY' ) goto 7001
c
c  --- make sure time span is correct ---
c
      if( begtim .EQ. 24.0 ) then
          begtim = 0.
          ibgdat = ibgdat + 1
      endif
      if( MOD(ibgdat,1000) .GT. 365 ) then
          if( MOD(INT(ibgdat/1000),4) .EQ. 0 ) then
             if( MOD(ibgdat,1000) .EQ. 367 )
     &                     ibgdat = (INT(ibgdat/1000)+1)*1000 + 1
          else
             ibgdat = (INT(ibgdat/1000)+1)*1000 + 1
          endif
      endif
      if( ibgdat .NE. ndate .OR. 
     &                        ABS(begtim-ttime) .GT. 0.00001 ) goto 7002 
c
c   --- region description header ----
c
      read(IORINI,ERR=7000) xutmin, yutmin, nzonin, xorgin, yorgin, 
     &                delxin, delyin, noxgin, noygin, nzclin, nvloin,
     &                                  nvupin, dzsfin, dzmlin, dzmuin
      lmatch = .TRUE.
      if( nzonin .NE. iuzon ) lmatch = .FALSE.
c
c   --- set scaling factor for coordinates ---
c
      if( .NOT.llatlon) then
         factr = 1000.0
      else
         factr = 1.0
      endif
      xorgin = xorgin/factr
      yorgin = yorgin/factr
      delxin = delxin/factr
      delyin = delyin/factr
      if( ABS(xorgin-xorg) .GT. 0.0001 ) lmatch = .FALSE.
      if( ABS(yorgin-yorg) .GT. 0.0001 ) lmatch = .FALSE.
      if( ABS(delxin-delx) .GT. 0.0001 ) lmatch = .FALSE.
      if( ABS(delyin-dely) .GT. 0.0001 ) lmatch = .FALSE.
      if( noxgin .NE. nox ) lmatch = .FALSE.
      if( noygin .NE. noy ) lmatch = .FALSE.
      if( nzclin .NE. noz ) lmatch = .FALSE.
      if( .NOT. lmatch ) goto 7003
c
c  --- segment description header ---
c
      read(IORINI,ERR=7000) 
c
c  ---- species list ---
c
      read(IORINI,ERR=7000) ((ispec(i,j),i=1,10),j=1,nspcs)
c
c  --- make sure species list is compatable, up to timing tracers ---
c
      do 10 l=1,ipttim-1
         write(cname,'(10A1)') (ispec(i,l),i=1,10)
         if( cname .NE. ptname(l) ) goto 7004
   10 continue
c
c  --- shift the species list around to get the previous timing tracers
c      in the list as well as the timing tracers for this run ---
c
      if( ntrtim .GT. 0 ) then
         ntoday = ntotsp - ipttim + 1
         ntotsp = nspcs + ntoday
         npttim = npttim + ntoday
         nreles = (nspcs - ipttim +1) / (2*nregin) 
         do 20 l=ntoday,1,-1
            ptname(nspcs+l) = ptname(ipttim-1+l)
   20    continue
         do 30 l=ipttim,nspcs
            write(cname,'(10A1)') (ispec(i,l),i=1,10)
            ptname(l) = cname
   30    continue
         nsaspc = nspcs
      endif
c
c   --- check for array overflow ---
c
      if( ntotsp .GT. MXTRSP ) goto 7005
c
c   --- read the data for this hour ----
c
      read(IORINI,ERR=7000) ibgdat, begtim, iendat, endtim
      do 40 l=1,nspcs
         do 50 k=1,noz
            read(IORINI) isegmt, (ispec(i,l),i=1,10), 
     &                         ((cnctmp(i,j),i=1,nox),j=1,noy)
            do 60 j=1,noy
               do 70 i=1,nox
                  saconc(i,j,k,l) = cnctmp(i,j)
   70          continue
   60       continue
   50    continue
   40 continue
c
c  --- initialize all of the tracers concs to zero to start off ---
c
      do 80 l=nspcs+1,nsaspc
         do 90 k=1,noz
             do 11 j=1,noy
                do 21 i=1,nox
                   saconc(i,j,k,l) = 0.
   21           continue
   11        continue
   90    continue
         do 31 i=1,MXRECP
            conrcp(l,i) = 0.
   31    continue
   80 continue
c
c  --- close file and return to calling routine ---
c
      close(IORINI)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in RDINSTSA:'
      write(iout,9000,ERR=9999)'Reading the header of the ',
     &    'tracer concentration initialization file for coarse grid: ',
     &                          inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in RDINSTSA:'
      write(iout,9000,ERR=9999) 'Tracer concentration ',
     &    'initilization file for coarse grid appears to be type: ',
     &                                                        cname
      write(iout,9000,ERR=9999) 'It should be type: AIRQUALITY'
      write(iout,9000,ERR=9999) 'Filename: ',
     &                          inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in RDINSTSA:'
      write(iout,9000,ERR=9999) 'Tracer concentration ',
     &   'initialization file for coarse grid is not for correct ',
     &                               'time period; ',
     &                          inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
      write(iout,9001,ERR=9999) 'Simulation Start: ',ndate,ttime
      write(iout,9001,ERR=9999) 'Initialization File: ',ibgdat,begtim
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in RDINSTSA:'
      write(iout,9000,ERR=9999) 'Region definition from ',
     &  ' tracer concentration file for coarse grid does not match ',
     &                                           'simulation control.'
      write(iout,9000,ERR=9999) '      Initialization file: ',
     &                          inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
      write(iout,9003,ERR=9999) 'UTM zone',nzonin
      write(iout,9002,ERR=9999) 'Region origin',xorgin,yorgin
      write(iout,9002,ERR=9999) 'Cell spacing',delxin, delyin
      write(iout,9004,ERR=9999) 'Number of cells',noxgin,noygin
      write(iout,9003,ERR=9999) 'Vertical cells',nzclin
      write(iout,9000,ERR=9999) 
      write(iout,9000,ERR=9999) '      Simulation Control'
      write(iout,9003,ERR=9999) 'UTM zone',iuzon
      write(iout,9002,ERR=9999) 'Region origin',xorg,yorg
      write(iout,9002,ERR=9999) 'Cell spacing',delx, dely
      write(iout,9004,ERR=9999) 'Number of cells',nox,noy
      write(iout,9003,ERR=9999) 'Vertical cells',noz
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in RDINSTSA:'
      write(iout,9000,ERR=9999) 'Species list for tracer ',
     &            'concentration initialization file for coarse grid',
     &          ' is inconsistent with user parameters ',
     &                          inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
      write(iout,9006,ERR=9999) 'Species number: ',l,' should be ',
     &                      ptname(l),' and ',cname,' was read instead.'
      call camxerr()
c
 7005 continue 
      write(iout,'(//,a)') 'ERROR in RDINSTSA:'
      write(iout,'(/,1X,A,I10)') 
     &          'Number of tracer species exceeds max: ',MXTRSP
      write(iout,'(1X,A)') 'Increase parameter: MXTRSP'
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in RDINSTSA:'
      write(iout,9000,ERR=9999) 'Opening the ',
     &    'tracer concentration initialization file for coarse grid: ',
     &                          inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
      write(iout,9005,ERR=9999) 'Make sure the name for output files ',
     &                           'reflects the correct simulation day.'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,1X,6A)
 9001 format(1X,A,I10.5,F10.1)
 9002 format(1X,A,T20,': (',F10.2,',',F10.2,')')
 9003 format(1X,A,T20,': ',I10)
 9004 format(1X,A,T20,': (',I10,',',I10,')')
 9005 format(10X,2A)
 9006 format(1X,A,I5,5A)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
