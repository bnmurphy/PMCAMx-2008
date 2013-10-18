      subroutine cncprep(endtim,enddate)
c 
c-----CAMx v4.02 030709
c 
c     CNCPREP reads the AIRQUALITY or INSTANT file header (depending on
c     whether this is a restart or not) and maps file species names to
c     the internal CAMx species list. The routine then writes headers
c     to new AVERAGE files
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        1/20/99   Grid cell size on file should be meters for all cartesian
c                  projections (UTM, LCP, PSP)
c        10/24/01  Removed BSWAP and converted integer strings to character*4
c        10/31/01  Added logic to ensure that a cold start reads an AIRQUALITY
c                  file and a restart reads an INSTANT file
c        11/06/01  Input dates are now Julian
c 
c     Input arguments: 
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        STARTUP 
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
c
      character*4 ifile(10),note(60),icspec(10,MXSPEC)
      integer enddate
      character*10 aqfil,infil,cnfil,avfil,icspc
      integer ifgptr(MXGRID),ifglvl(MXGRID)
c
      data aqfil /'AIRQUALITY'/
      data cnfil /'INSTANT   '/
      data avfil /'AVERAGE   '/
      data nseg,izero /1,0/
      data zero /0./
c
c-----Entry point
c
      lairqul = .false.
      iunit = iic
      if (lrstrt) iunit = irstc
c
c-----Read 1st IC header record and check inputs
c
      read(iunit,ERR=7000,END=7000) ifile,note,iseg,nicspc,idat1,
     &                                             tim1,idat2,tim2
      if (INT(tim2) .EQ. 24 ) then
        idat2 = idat2 + 1
        tim2 = 0.
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      write(infil,'(10a1)') (ifile(n),n=1,10)
      if (.not.lrstrt .and. infil.ne.aqfil) then
        write(iout,'(//,a)') 'ERROR in CNCPREP:'
        write(iout,*)'This is a cold start from Initial Conditions'
        write(iout,*)'IC input file is not labelled AIRQUALITY'
        call camxerr()
      endif
      if (lrstrt .and. infil.ne.cnfil) then
        write(iout,'(//,a)') 'ERROR in CNCPREP:'
        write(iout,*)'This is a restart from a Restart File'
        write(iout,*)'IC input file is not labelled INSTANT'
        call camxerr()
      endif
      if (infil.eq.aqfil) lairqul = .true.
      if (nicspc.gt.MXSPEC) then
        write(iout,'(//,a)') 'ERROR in CNCPREP:'
        write(iout,*)'Number of species on IC file > MXSPEC'
        write(iout,*) nicspc,MXSPEC
        call camxerr()
      endif
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if (idat1.gt.begdate) then
        write(iout,'(//,a)') 'ERROR in CNCPREP:'
        write(iout,*)'IC start date > simulation start date'
        write(iout,'(a,i10.5,a,i10.5)')
     &        'IC file: ',idat1,' Sim start: ',begdate
        call camxerr()
      elseif (idat1.eq.begdate .and. tim1.gt.begtim) then
        write(iout,'(//,a)') 'ERROR in CNCPREP:'
        write(iout,*)'IC start time > simulation start time'
        write(iout,*)'IC file: ',tim1,' Sim start: ',begtim
        call camxerr()
      endif
c
c-----Read 2nd IC header record and check inputs
c
      read(iunit,ERR=7001) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz
      if (.NOT.llatlon) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if (dx.ne.delx .or. dy.ne.dely) then
        write(iout,'(//,a)') 'WARNING in CNCPREP:'
        write(iout,*)'IC cell size not equal to model cell size'
        write(iout,*)'IC file: ',dx,dy,' model: ',delx,dely
      elseif (nx.ne.ncol(1) .or. ny.ne.nrow(1) 
     &                          .or. nz.ne.nlay(1)) then
        write(iout,'(//,a)') 'ERROR in CNCPREP:'
        write(iout,*)'IC grid size not equal to model grid size'
        write(iout,*)'IC file: ',nx,ny,nz,
     &               ' model: ',ncol(1),nrow(1),nlay(1)
        call camxerr()
      endif 
c
c-----Read 3rd & 4th IC header 
c
      read(iunit,ERR=7001) (idum,idum,idum,idum,n=1,iseg)
      read(iunit,ERR=7001) ((icspec(n,l),n=1,10),l=1,nicspc)
c
c-----Map IC species to model species
c
      do 20 l = 1,nspec
        licmap(l,1) = 0
        do 15 lic = 1,nicspc
          write(icspc,'(10a1)') (icspec(n,lic),n=1,10)
          if (icspc.eq.'HNO2      ') icspc = 'HONO      '
          if (icspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        icspc = 'FORM      '
          if (icspc.eq.spname(l)) then
            licmap(l,1) = lic
            write(idiag,'(2(a,i5,2x,a))')
     &                   'Initial species ',lic,icspc,
     &                   ' mapped to model species ',l,spname(l)
            goto 20
          endif
 15     continue
        write(idiag,*)'Did not find species: ',spname(l),' on IC file'
        if (.not.lairqul) then
          write(iout,'(//,a)') 'ERROR in CNCPREP:'
          write(iout,*)'The INSTANT file must contain the same ',
     &                 'species as specified in the CHEMPARAM file'
          call camxerr()
        endif
 20   continue
      write(idiag,*)
c 
c-----Write average output concentration file headers
c
      idat1 = begdate
      idat2 = enddate
      tim1 = begtim/100.
      tim2 = endtim/100.
      if (.NOT.llatlon) then
        orgx = 1000.*xorg
        orgy = 1000.*yorg
        dx = 1000.*delx
        dy = 1000.*dely
        izone = 0
        if (lutm) izone = iuzon
      else
        orgx = xorg
        orgy = yorg
        dx = delx
        dy = dely
        izone = 0
      endif
      read(runmsg(1:60),'(60a1)') (note(n),n=1,60)
      read(avfil,'(10a1)') (ifile(n),n=1,10)
      do l = 1,navspc
        read(spname(lavmap(l)),'(10a1)') (icspec(n,l),n=1,10)
      enddo
      if (l3davg) then
        nlayer = nlay(1)
      else
        nlayer = 1 
      endif
c
c-----Coarse grid average header
c
      rewind(iavg)
      write(iavg) ifile,note,nseg,navspc,idat1,tim1,idat2,tim2
      write(iavg) zero,zero,izone,orgx,orgy,dx,dy,nx,ny,nlayer,
     &            izero,izero,zero,zero,zero
      write(iavg) izero,izero,nx,ny
      write(iavg) ((icspec(n,l),n=1,10),l=1,navspc)
      if( ngrid .EQ. 1 ) goto 9999
c
c-----Fine grid average header
c
      do i = 1,ngrid 
        ifglvl(i) = 0 
      enddo 
      do i = 1,ngrid 
        do j = 1,nchdrn(i) 
          idch = idchdrn(j,i) 
          ifgptr(idch) = i - 1 
          ifglvl(idch) = ifglvl(idch) + 1 
          do k = 1,nchdrn(idch) 
            ifglvl(idchdrn(k,idch)) = ifglvl(idchdrn(k,idch)) + 1 
          enddo 
        enddo 
      enddo
      rewind(ifavg)
      write(ifavg) runmsg 
      write(ifavg) ngrid-1, navspc 
      write(ifavg) (spname(lavmap(i)),i=1,navspc) 
      do igrd = 2,ngrid
        if (l3davg) then
          nlayer = nlay(igrd)
        else
          nlayer = 1 
        endif
        write(ifavg) inst1(igrd),jnst1(igrd),inst2(igrd), 
     &               jnst2(igrd),meshold(igrd),meshold(igrd), 
     &               ncol(igrd),nrow(igrd),nlayer, 
     &               ifgptr(igrd),ifglvl(igrd) 
      enddo
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in CNCPREP:'
      write(iout,'(A)',ERR=9999)'Reading restart file for coarse grid.'      
      write(iout,'(2A)',ERR=9999)'Make sure the filename is specified ',
     &          'correctly and the previous day finished correctly.'
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in CNCPREP:'
      write(iout,'(2A)',ERR=9999)'Reading the header of restart file ',
     &                                  'for coarse grid.'      
      call camxerr()
c
 9999 continue
      return
      end
