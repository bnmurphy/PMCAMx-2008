      program CAMx
c
c***********************************************************************
c
c                  CCCCCCC      AA      MM      MM                       
c                 CC          AA  AA    MMM    MMM   xx   xx  
c                 CC         AA    AA   MM MMMM MM    xx xx  
c                 CC         AAAAAAAA   MM  MM  MM     xxx  
c                 CC         AA    AA   MM      MM    xx xx
c                  CCCCCCC   AA    AA   MM      MM   xx   xx
c
c        C O M P R E H E N S I V E   A I R   Q U A L I T Y   M O D E L
c        -                           -                       - 
c                         with   E X T E N S I O N S
c                                  -
c
c                          VERSION  4.02  03-07-09
c
c           Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c                     ENVIRON International Corporation 
c                        101 Rowland Way, Suite 220 
c                          Novato, CA 94945-5010
c                             (415) 899 - 0700
c                                www.camx.com
c
c     Revision History:
c       07/09/03   Version 4.02 released (see Release_notes_4.02)
c       06/17/03   Version 4.01 released (see Release_notes_4.01)
c       05/31/03   Version 4.00 released (see Release_notes_4.00)
c       11/15/01   Version 3.10 released (see Release_notes_3.10)
c        4/25/01   Version 3.01 released (see Release_notes_3.01)
c       12/30/00   Version 3.00 released (see Release_notes_3.00) 
c       12/30/99   Version 2.03 released (see Release_notes_2.03) 
c        9/13/99   Version 2.02 released (see Release_notes_2.02) 
c        5/10/99   Version 2.01 released (see Release_notes_2.01) 
c       12/30/98   Version 2.00 released, original development
c
c     Routines Called: 
c        STARTUP,  READINP,  INTRPDAT, READCNC,  WRTMASS,  MASSUM,
c        TIMESTEP, KHORZ,    DRYDEP,   WRTDEP,   READPT,   READAR,
c        READBND,  READAHO,  IASSGN2D, AVERAGE,  UPTIME,   PIGINIT,
c        RADDRIVR, EMISTRNS, NESTING,  CHEMRXN,  PIGWALK,  PIGMSCAL,
c        WRTCON,   WRFGCON,  WRTPIG,   PIGDIAG
c
c***********************************************************************
c  
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'grid.com'
      include "ptemiss.com"
      include 'bndary.com'
      include 'chmstry.com'
      include 'flags.com'
      include 'ahomap.com'
      include 'filunit.com'
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
      include 'rtracchm.com'
c
c========================= Source Apportion End ========================
c
c
c========================= Process Analysis Begin ==============================
c
      include 'procan.com'
c
c     linit   --  logical, if true, initialize   cipr(IPR_INIT ,*,*,*)
c                          otherwise, initialize cipr(IPR_FINAL,*,*,*)
c
      logical linit
c
c========================== Process Analysis End ===============================
c
      integer inpdate,emsdate,hazdate,ozndate,bnddate,wrtdate,enddate 
      integer ldumap(MXSPEC)
      real inptim,emstim,haztim,ozntim,bndtim,wrttim,endtim 
      character*20 version
      character*10 name
c
      data version /'CAMx v4.02, 03-07-09'/
c
c-----Entry point
c
      nsteps = 0
      do i=1,MXSPEC
         ldumap(i) = i
      enddo
c
c-----Start model simulation
c
      write(*,'(a20,$)') 'startup ......'
      call flush(6)
      call startup(version,inptim,inpdate,emstim,emsdate,haztim,
     &             hazdate,ozntim,ozndate,bndtim,bnddate,wrttim,
     &             wrtdate,endtim,enddate) 
c
c-----Set average concentrations to zero
c
      do igrd = 1,ngrid
        nodes=ncol(igrd)*nrow(igrd)*nlay(igrd)*navspc
        call zeros(avcnc(iptr4d(igrd)),nodes)
        nodes = ncol(igrd)*nrow(igrd)*3*navspc
        call zeros(depfld(iptrdp(igrd)),nodes)
      enddo
c
c======================== Source Apportion Begin =======================
c
c   --- call routine to initialize the running averages ---
c
      if( ltrace .OR. lddm ) then
         do igrd=1,ngrid
            nodes=ncol(igrd)*nrow(igrd)*ntotsp
            call zeros(ptavrg(ipsa2d(igrd)),nodes)
         enddo
         nodes=ncol(1)*nrow(1)*MXSPEC
         call zeros(modavg,nodes)
         nodes = MXRTCEL* MXSPEC 
         call zeros(rcpdcy,nodes)
      endif
c
c========================= Source Apportion End ========================
c
c
      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      tcpu = etime(tarray)
      write(iout,*)
      write(iout,'(a,f10.3)') 
     &       'Elapsed time (s) for model startup: ',tcpu
      write(iout,*)
      write(idiag,*)
      write(iout,'(a)')'BEGIN SIMULATION'
      write(idiag,'(a)')'BEGIN SIMULATION'
      write(iout,'(a)') 'Starting main integration loop...'
      write(idiag,*)
c
c--------------------  Main time-integration loop  ---------------------
c
  100 continue
      nsteps = nsteps + 1
      xyordr = mod(nsteps,2)
      write(iout,*)
      write(iout,'(a,i8.5,f8.2)') 'Date/time: ',date,time
      write(iout,*)
      write(*,*)
      write(*,'(a,i8.5,f8.2)') 'Date/time: ',date,time
      write(*,*)
c
c-----Check if time-varying grid/met data are to be read
c
      if (date.eq.inpdate .and. abs(time-inptim).lt.0.01) then
        write(*,'(a20,$)') 'readinp ......'
        tcpu = dtime(tarray2)
        do 10 igrd = 1,ngrid
          call readinp(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                 height(iptr3d(igrd)),phpt(iptr3d(igrd)),
     &                 hnxt(iptr3d(igrd)),
     &                 press(iptr3d(igrd)),pppt(iptr3d(igrd)),
     &                 pnxt(iptr3d(igrd)),
     &                 windu(iptr3d(igrd)),pupt(iptr3d(igrd)),
     &                 unxt(iptr3d(igrd)),
     &                 windv(iptr3d(igrd)),pvpt(iptr3d(igrd)),
     &                 vnxt(iptr3d(igrd)),
     &                 tsurf(iptr2d(igrd)),pspt(iptr2d(igrd)),
     &                 tsnxt(iptr2d(igrd)),
     &                 tempk(iptr3d(igrd)),ptpt(iptr3d(igrd)),
     &                 tnxt(iptr3d(igrd)),
     &                 water(iptr3d(igrd)),fcloud(iptr3d(igrd)),
     &                 rkv(iptr3d(igrd)),cwc(iptr3d(igrd)),
     &                 pwc(iptr3d(igrd)),cod(iptr3d(igrd)),
     &                 cldtrns(iptr3d(igrd)),pi0(iptr3d(igrd)) ) 
 10     continue
c
c------Estimate those fields that were not read
c
        call intrpdat
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        whr = aint(inptim/100.)
        wmn = amod(inptim,100.)
        inptim = 100.*(whr + aint((wmn + dtinp)/60.)) + 
     &             amod((wmn + dtinp),60.)
        if (inptim.ge.2400.) then
          inptim = inptim - 2400.
          inpdate = inpdate + 1
          if( MOD(inpdate,1000) .GT. 365 ) then
            if( MOD(INT(inpdate/1000),4) .EQ. 0 ) then
               if( MOD(inpdate,1000) .EQ. 367 )
     &                     inpdate = (INT(inpdate/1000)+1)*1000 + 1
            else
               inpdate = (INT(inpdate/1000)+1)*1000 + 1
            endif
         endif
        endif
c
c-----Initialize concentrations from an AIRQUALITY file or RESTART files
c
        if (nsteps.eq.1) then
          write(*,'(a20,$)') 'readcnc ......'
          call readcnc
          do igrd = 1,ngrid
            call wrtmass(igrd,date,time,0)
            call massum(igrd,nspec,ncol(igrd),nrow(igrd),nlay(igrd),
     &                  deltax(1,igrd),deltay(igrd),depth(iptr3d(igrd)),
     &                  conc(iptr4d(igrd)),xmass0(1,igrd))
            do l = 1,nspec
              xmsold(l,igrd) = xmass0(l,igrd)
            enddo
c
c --- initialize the mass summary array for PiG ---
c
            if( ipigflg .EQ. GRESPIG ) then
               do i=1,4
                  pigdump(i,igrd) = 0.0
               enddo
            else if( ipigflg .EQ. IRONPIG ) then
               do i=1,nspec
                  pigdump(i,igrd) = 0.0
               enddo
            endif
          enddo
c
c======================== Source Apportion Begin =======================
c
c   --- if this is not a restart, call routine to fill the initial
c       conditions arrays for the tracers ---
c
           if( ltrace .AND. .NOT. lrstrt ) then
               if( tectyp .EQ. RTRAC ) then
                  call rdicrt(ncol(1),nrow(1),nlay(1),ntotsp,ptconc(1))
c
c   --- interpolate initial conditions to nests ---
c
                  if(ngrid.gt.1) then
                     do ip = 1,ngrid
                       do ic = 1,nchdrn(ip)
                         ig = idchdrn(ic,ip)
                         call intrpcnc(ntotsp,ncol(ip),nrow(ip),
     &                               nlay(ip),i1(ig),j1(ig),nmesh(ig),
     &                               nmshv(1,ig),ncol(ig),nrow(ig),
     &                               nlay(ig),ptconc(ipsa3d(ip)),
     &                               ptconc(ipsa3d(ig)) )
                       enddo
                     enddo
                  endif
c
c --- call routine to convert IC to ug/m ----
c
                  do igrd=1,ngrid
                       call cvticrt(igrd,ncol(igrd),nrow(igrd),
     &                          nlay(igrd),ntotsp,ptconc(ipsa3d(igrd)),
     &                          tempk(iptr3d(igrd)),press(iptr3d(igrd)))
                  enddo
               else
                  do 20 igrd=1,ngrid
                     call filaqsa(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                          nspec,ntotsp,conc(iptr4d(igrd)),
     &                          ptconc(ipsa3d(igrd)) )
   20             continue
               endif
           endif
c
c========================= Source Apportion End ========================
c
c
c======================== DDM Begin =======================
c
c   --- call routine to read the IC file for DDM ---
c
           if( lddm .AND. .NOT. lrstrt .AND. nicddm .GT. 0) then
              call rdicddm(ncol(1),nrow(1),nlay(1),ntotsp,ptconc(1))
c
c   --- interpolate initial conditions to nests ---
c 
              if(ngrid.gt.1) then
                 do ip = 1,ngrid
                   do ic = 1,nchdrn(ip)
                     ig = idchdrn(ic,ip)
                     call intrpcnc(ntotsp,ncol(ip),nrow(ip),nlay(ip),
     &                               i1(ig),j1(ig),nmesh(ig),
     &                               nmshv(1,ig),ncol(ig),nrow(ig),
     &                               nlay(ig),ptconc(ipsa3d(ip)),
     &                               ptconc(ipsa3d(ig)) )
                   enddo
                 enddo
              endif
c
c --- call routine to convert IC to ug/m ----
c
              do igrd=1,ngrid
                   call cvticddm(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                        ntotsp, ptconc(ipsa3d(igrd)),
     &                          tempk(iptr3d(igrd)),press(iptr3d(igrd)))
              enddo
           endif
c
c========================= DDM End ========================
c
c
c========================= Process Analysis Begin ==============================
c
          if( lipr .OR. lirr ) then
             call pazero()
             if( lipr ) then
               linit = .TRUE.
               do igrd = 1,ngrid
                  call initipr(linit,igrd,nspec,ncol(igrd),nrow(igrd),
     &                                   nlay(igrd),conc(iptr4d(igrd)))
               enddo
             endif
          endif
c
c========================== Process Analysis End ===============================
c
          tcpu = dtime(tarray2)
          write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        endif
c
c-----Calculate timestep
c
        write(*,'(a20,$)') 'timestep ......'
        call timestep()
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
c
c-----Calculate horizontal diffusion coefficients
c
        write(*,'(a20,$)') 'khorz ......'
        do igrd = 1,ngrid
          call khorz(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &               deltax(1,igrd),deltay(igrd),deltat(igrd),
     &               windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &               idfin(iptr2d(igrd)),rkx(iptr3d(igrd)),
     &               rky(iptr3d(igrd)) )
        enddo
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
c
c-----Calculate dry deposition rates
c
        if (ldry) then
          write(*,'(a20,$)') 'drydep ......'
          do igrd = 1,ngrid
            call drydep(igrd,ncol(igrd),nrow(igrd),nlay(igrd),itzon,
     &                  tsurf(iptr2d(igrd)),cellat(iptr2d(igrd)),
     &                  cellon(iptr2d(igrd)),pwc(iptr3d(igrd)),
     &                  cwc(iptr3d(igrd)),height(iptr3d(igrd)),
     &                  press(iptr3d(igrd)),
     &                  windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &                  fcloud(iptr3d(igrd)),cldtrns(iptr3d(igrd)),
     &                  water(iptr3d(igrd)),fsurf(iptrlu(igrd)),
     &                  tempk(iptr3d(igrd)),
cbk     &                  vdep(iptrem(igrd)) )
     &                  vdep(iptrem(igrd)),conc(iptr4d(igrd)) )
            call depsmry(igrd,ncol(igrd),nrow(igrd),nspec,
     &                   vdep(iptrem(igrd)) )
c
c======================== Source Apportion Begin =======================
c
            if( ltrace .AND. tectyp .EQ. RTRAC ) then
                call drydeprt(igrd,ncol(igrd),nrow(igrd),nlay(igrd),
     &                  nrtrac,itzon,                 
     &                  tsurf(iptr2d(igrd)),cellat(iptr2d(igrd)),
     &                  cellon(iptr2d(igrd)),pwc(iptr3d(igrd)),
     &                  cwc(iptr3d(igrd)),height(iptr3d(igrd)),
     &                  press(iptr3d(igrd)),
     &                  windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &                  fcloud(iptr3d(igrd)),
     &                  water(iptr3d(igrd)),fsurf(iptrlu(igrd)),
     &                  tempk(iptr3d(igrd)),
     &                  vdeprt(ipsa2d(igrd)))
            endif
c
c======================== Source Apportion End =======================
c
          enddo
          tcpu = dtime(tarray2)
          write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        else
          do igrd = 1,ngrid
            call zeros(vdep(iptrem(igrd)),ncol(igrd)*nrow(igrd)*nspec)
          enddo
        endif
      endif
c
c-----Check if emissions data are to be read
c
      if (date.eq.emsdate .and. abs(time-emstim).lt.0.01) then
        if (lptsrc) then
          write(*,'(a20,$)') 'readpt ......'
          call readpt
c
c======================== Source Apportion Begin =======================
c
c   --- call routine to read the points source emissions files
c       and load the tracer emission arrays ---
c
          if( ltrace ) then
            if( tectyp .EQ. RTRAC ) then
                call rdptrt(date,time)
            else
                call readptsa(date,time)
            endif
          endif
c
c========================= Source Apportion End ========================
c
c
c======================== DDM Begin =======================
c
c   --- call routine to read the points source emissions files
c       and load the tracer emission arrays ---
c
          if( lddm ) then
            call rdptddm(date,time)
          endif
c
c========================= DDM End ========================
c
          tcpu = dtime(tarray2)
          write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        endif
c
        if( larsrc ) then
c
c------Read area emissions if data is available----
c
           write(*,'(a20,$)') 'readar ......'
           do igrd = 1,ngrid
              if( iarem(1,igrd) .GT. 0 ) then
                 call readar(igrd,ncol(igrd),nrow(igrd),
     &                            iarem(1,igrd),iout,
     &                            aremis(iptrem(igrd)),
     &                            deltax(1,igrd),deltay(igrd))
              else
c
c------Otherwise assign values from the parent---
c
                 do ip=1,ngrid
                   do ic = 1,nchdrn(ip)
                     if( igrd .EQ. idchdrn(ic,ip) ) then
                        write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                         'Assigning emissions from parent grid',
     &                                    time, date,' grid',igrd
                        call emassign(ncol(ip),nrow(ip),
     &                    i1(igrd),j1(igrd),nmesh(igrd),ncol(igrd),
     &                        nrow(igrd),narspc,aremis(iptrem(ip)),
     &                                         aremis(iptrem(igrd)) )
                     endif
                   enddo
                 enddo
              endif
           enddo
c
c======================== Source Apportion Begin =======================
c
c   --- call routine to read the emissions files and load
c       the tracer emission arrays ---
c
           if( ltrace ) then
              do igrd=1,ngrid
                if( tectyp .EQ. RTRAC ) then
                    call rdarrt(igrd,date,time,ncol(igrd),nrow(igrd),
     &                                    ntotsp,saemis(ipsa2d(igrd)) )
                else
                    call readarsa(igrd,date,time,
     &                        ncol(igrd),nrow(igrd),ntotsp,
     &                        deltax(1,igrd),deltay(igrd),
     &                        saemis(ipsa2d(igrd)) )
                endif 
              enddo
           endif
c
c========================= Source Apportion End ========================
c
c
c======================== DDM Begin =======================
c
c   --- call routine to read the emissions files and load
c       the DDM tracer emission arrays ---
c
          if( lddm .AND. nemddm .GT. 0) then
             do igrd=1,ngrid
               call rdarddm(igrd,date,time,ncol(igrd),nrow(igrd),ntotsp,
     &                                            saemis(ipsa2d(igrd)) )
             enddo
          endif
c
c========================= DDM End ========================
c
          tcpu = dtime(tarray2)
          write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        endif
c
        whr = aint(emstim/100.)
        wmn = amod(emstim,100.)
        emstim = 100.*(whr + aint((wmn + dtems)/60.)) +
     &             amod((wmn + dtems),60.)
        if (emstim.ge.2400.) then
          emstim = emstim - 2400.
          emsdate = emsdate + 1
          if( MOD(emsdate,1000) .GT. 365 ) then
             if( MOD(INT(emsdate/1000),4) .EQ. 0 ) then
                if( MOD(emsdate,1000) .EQ. 367 )
     &                     emsdate = (INT(emsdate/1000)+1)*1000 + 1
             else
                emsdate = (INT(emsdate/1000)+1)*1000 + 1
             endif
          endif
        endif
      endif
c
c-----Check if coarse grid boundary data are to be read
c
      if (date.eq.bnddate .and. abs(time-bndtim).lt.0.01) then
        write(*,'(a20,$)') 'readbnd ......'
        call readbnd(bndtim,bnddate)
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
c
c======================== Source Apportion Begin =======================
c
        if( ltrace .AND. tectyp .EQ. RTRAC ) then
        endif
c
c======================== Source Apportion End =======================
c
      endif
c
c======================== Source Apportion Begin =======================
c
c   --- call routine to clear the old boundary cells and then
c       call routine to fill with new boundary concentrations ----
c       For the Coarse Grid ---
c
      if( ltrace ) then
         if( tectyp .EQ. RTRAC ) then
            call rdbcrt(nsteps,ncol(1),nrow(1),nlay(1),ntotsp,ptconc(1),
     &                                                tempk(1),press(1))
         else
            call clrbdysa(1,ncol(1),nrow(1),nlay(1),ntotsp,ptconc(1))
            call filbdysa(1,ncol(1),nrow(1),nlay(1),nspec,ntotsp,
     &                                       conc(1),ptconc(1))
c
c   --- make sure tracer concs are greater than lower bound ---
c
            do 30 i=1,MXSA3D
               ptconc(i) = AMAX1(ptconc(i),BNDLPT)
   30       continue
         endif
      endif
c
c========================= Source Apportion End ========================
c
c
c======================== DDM Begin =======================
c
c   --- call routine to clear the old boundary cells and then
c       call routine to fill with new boundary concentrations ----
c       For the Coarse Grid ---
c
      if( lddm .AND. nbcddm .GT. 0) then
        call clrbdyddm(ncol(1),nrow(1),nlay(1),ntotsp,ptconc(1))
        call rdbcddm(ncol(1),nrow(1),nlay(1),ntotsp,ptconc(1),
     &                                            tempk(1),press(1))
      endif
c
c========================= DDM End ========================
c
c-----Check if haze data are to be read
c
      if (lchem) then
        if (date.eq.hazdate .and. abs(time-haztim).lt.0.01) then
          name = 'HAZE      '
          write(*,'(a20,$)') 'readaho (haze)'
          call readaho(ncol(1),nrow(1),time,date,ly2k,name,haztim,
     &                                                hazdate,icdhaz)
          tcpu = dtime(tarray2)
          write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        endif
c
c-----Check if ozone column data are to be read
c
        if (date.eq.ozndate .and. abs(time-ozntim).lt.0.01) then
          name = 'OZONE COL '
          write(*,'(a20,$)') 'readaho (o3)..'
          call readaho(ncol(1),nrow(1),time,date,ly2k,name,ozntim,
     &                                               ozndate,icdozn)
          tcpu = dtime(tarray2)
          write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        endif
c
c-----Assign haze and ozone column values for fine grids
c
        do ip = 1,ngrid
          do ic = 1,nchdrn(ip)
            igrd = idchdrn(ic,ip)
            call iassgn2d(ncol(ip),nrow(ip),i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    icdhaz(iptr2d(ip)),icdhaz(iptr2d(igrd)))
            call iassgn2d(ncol(ip),nrow(ip),i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    icdozn(iptr2d(ip)),icdozn(iptr2d(igrd)))
          enddo
        enddo
      endif
c
c----- Update coarse grid average concentrations
c
      call average(.FALSE.,1,deltat(1)/2.0,ncol(1),nrow(1),
     &             nlay(1),nlay(1),navspc,nspec,lavmap,tempk(1),
     &             press(1),conc(1),avcnc(1),ipacl_3d(1))
c
c----- Add PiG masses to average
c
      if( ipigflg .EQ. IRONPIG ) then
         call avepig(1,deltat(1)/2.0,ncol(1),nrow(1),nlay(1),
     &            nlay(1),deltax(1,1),deltay(1),mapscl(1),
     &            height(1),navspc,nspec,lavmap,tempk(1),press(1),
     &            conc(1),avcnc(1))
      endif
c
c======================== Source Apportion Begin =======================
c======================== DDM Begin =======================
c
c   --- call routine to update the running averages ---
c
      if( ltrace .OR. lddm ) then
         call average(.TRUE.,1,deltat(1)/2.0,ncol(1),nrow(1),
     &               nlay(1),1,ntotsp,ntotsp,lsamap,tempk(1),
     &                         press(1),ptconc(1),ptavrg(1),ipacl_3d(1))
         call average(.TRUE.,1,deltat(1)/2.0,ncol(1),nrow(1),
     &                  nlay(1),1,nspec,nspec,ldumap,tempk(1),
     &                              press(1),conc(1),modavg,ipacl_3d(1))
c
c   --- if WALL OF CELLS receptors exist, add averages ---
c
         if( lwalls ) then
           do igrd=1,ngrid
            call avgwal(igrd,ncol(igrd),nrow(igrd),nlay(igrd),nspec,
     &                  ntotsp,deltat(igrd)/2.0,deltax(1,igrd),
     &                  deltay(igrd),depth(iptr3d(igrd)),
     &                  tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &                  ptconc(ipsa3d(igrd)),conc(iptr4d(igrd)))
          enddo
        endif
      endif
c
c========================= Source Apportion End ========================
c======================== DDM End =======================
c
c-----Update date/time for this timestep
c
      call uptime(time,date,deltat(1))
c
c-----Initialize new pig puffs
c
      if( ipigflg .EQ. IRONPIG .OR. ipigflg .EQ. GRESPIG ) then
        write(*,'(a20,$)') 'piginit ......'
        write(iout,'(/,a20,$)') 'piginit ......'
        call piginit(deltat(1),ptemis,height,windu,windv,tempk,press)
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
      endif
c
c-----Update PiG locations 
c
      if( ipigflg .EQ. IRONPIG .OR. ipigflg .EQ. GRESPIG ) then
        write(*,'(a20,$)') 'pigwalk ......'
        write(iout,'(a20,$)') 'pigwalk ......'
        call pigwalk(deltat(1))
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
      endif
c
c-----Initialize radical concentrations
c
      if (lchem.and.nsteps.eq.1) then
        write(*,'(a20,$)') 'raddrivr ......'
        write(iout,'(a20,$)') 'raddrivr ......'
        call raddrivr
        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
      endif
c
c-----Coarse grid emission and transport
c
      call emistrns(1)
c
c-----Computation for nested grids
c
      if (ngrid.gt.1) call nesting
c     
c-----Coarse grid chemical reactions
c
      call chemrxn(1)
c
c-----Compute PiG mass
c
      if( ipigflg .EQ. IRONPIG .OR. ipigflg .EQ. GRESPIG ) then
        if( ipigflg .EQ. GRESPIG ) then
           call gresmscl(ngrid,time,date,idiag,pigdump)
        else if( ipigflg .EQ. IRONPIG ) then
           call ironmscl(ngrid,time,date,idiag,pigdump)
        endif
      endif
c
c-----Update average concentrations
c
      call average(.FALSE.,1,deltat(1)/2.0,ncol(1),nrow(1),
     &             nlay(1),nlay(1),navspc,nspec,lavmap,tempk(1),
     &             press(1),conc(1),avcnc(1),ipacl_3d(1))
c
c----- Add PiG masses to average
c
      if( ipigflg .EQ. IRONPIG ) then
         call avepig(1,deltat(1)/2.0,ncol(1),nrow(1),nlay(1),
     &            nlay(1),deltax(1,1),deltay(1), mapscl(1),
     &            height(1),navspc,nspec,lavmap,tempk(1),press(1),
     &            conc(1),avcnc(1))
      endif
c
c======================== Source Apportion Begin =======================
c
c   --- call routine to update the running averages ---
c
      if( ltrace .OR. lddm ) then
         call average(.TRUE.,1,deltat(1)/2.0,ncol(1),nrow(1),
     &                      nlay(1),1,ntotsp,ntotsp,lsamap,tempk(1),
     &                       press(1),ptconc(1),ptavrg(1),ipacl_3d(1))
         call average(.TRUE.,1,deltat(1)/2.0,ncol(1),nrow(1),
     &                      nlay(1),1,nspec,nspec,ldumap,tempk(1),
     &                            press(1),conc(1),modavg,ipacl_3d(1))
c
c   --- if WALL OF CELLS receptors exist, add averages ---
c
         do igrd=1,ngrid
            call avgwal(igrd,ncol(igrd),nrow(igrd),nlay(igrd),nspec,
     &                  ntotsp,deltat(igrd)/2.0,deltax(1,igrd),
     &                  deltay(igrd),depth(iptr3d(igrd)),
     &                  tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &                  ptconc(ipsa3d(igrd)),conc(iptr4d(igrd)))
         enddo
      endif
c
c========================= Source Apportion End ========================
c
c-----Check if concentration fields are to be written
c
      if (date.eq.wrtdate .and. abs(time-wrttim).lt.0.01) then
c
c======================== Source Apportion Begin =======================
c======================== DDM Begin =======================
c
        if( ltrace .OR. lddm ) then
            if( tectyp .EQ. RTRAC ) then
               call wrrcprt(date,time)
            else
c
c   --- call routine to find the peak cell in the coarse grid ---
c
               call pekrcp(avcnc(1),ncol(1),nrow(1),nlay(1),navspc)
c
c   --- call routine to get the averages at the receptors and
c       write the receptor average file ----
c
               do igrd=1,ngrid
                  call addrcp(igrd,ncol(igrd),nrow(igrd),
     &                    ntotsp,nspec,ptavrg(ipsa2d(igrd)),modavg)
               enddo
               if( ltrace ) call avgrcp(date,time,ncol(1),nrow(1),
     &                                 ntotsp,nspec,ptavrg(1),modavg)
               if( lddm ) call avgrcpddm(date,time,ncol(1),nrow(1),
     &                                                  nspec,modavg)
            endif
        endif
c
c   --- call routine to write the tracer instantaneous files and the
c       gridded surface concentrations file ---
c
        if( ltrace .OR. lddm ) then
            call instsa(1,date,time,ncol(1),nrow(1),nlay(1),ntotsp,
     &                                         ptconc(1),ptavrg(1))
            if( ngrid .GT. 1 ) call wrfgsa(date,time)
c
c   --- call routine to re-initialize the running averages ---
c
            do igrd=1,ngrid
              nodes=ncol(igrd)*nrow(igrd)*ntotsp
              call zeros(ptavrg(ipsa2d(igrd)),nodes)
            enddo
            nodes=ncol(1)*nrow(1)*MXSPEC
            call zeros(modavg,nodes)
        endif
c
c========================= Source Apportion End ========================
c======================== DDM End =======================
c
c-----Coarse grid files
c
        write(*,'(a20,$)') 'wrtcon ......'
        write(iout,'(a20,$)') 'wrtcon ......'
        call wrtcon(0,time,date,iavg,ncol(1),nrow(1), 
     &              nlay(1),navspc,avcnc(1))
        if (mod(int(time/100.),2).eq.1) then
          iunit = iconc(1)
        else
          iunit = iconc(2)
        endif
        call wrtcon(1,time,date,iunit,ncol(1),nrow(1),
     &              nlay(1),nspec,conc(1))
        if (ldry) then
          call wrtdep(time,date,idep,ncol(1),nrow(1),3*navspc,nspec,
     &                vdep(1),depfld(1))
        endif
c
c-----Fine grid files, and zero out average concentration array
c
        if( ngrid .GT. 1) then
          call wrfgcon(date,time)
          if (ldry) call wrfgdep(date,time)
        endif
        do igrd = 1,ngrid
          nodes = ncol(igrd)*nrow(igrd)*nlay(igrd)*navspc
          call zeros(avcnc(iptr4d(igrd)),nodes)
          nodes = ncol(igrd)*nrow(igrd)*3*navspc
          call zeros(depfld(iptrdp(igrd)),nodes)
        enddo

        tcpu = dtime(tarray2)
        write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
        write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
c
c-----Write PiG restart file and diagnostics
c
        if( ipigflg .EQ. IRONPIG .OR. ipigflg .EQ. GRESPIG ) then
          write(*,'(a20,$)') 'wrtpig ......'
          write(iout,'(a20,$)') 'wrtpig ......'
          call wrtpig(date,time,begdate,begtim)
          tcpu = dtime(tarray2)
          write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
          write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
        endif
c
c========================= Process Analysis Begin ==============================
c
c-----Get final concentration
c
        if (lipr .OR. lirr) then
           if( lipr ) then
             linit = .FALSE.
             do igrd = 1,ngrid
                call initipr(linit,igrd,nspec,ncol(igrd),nrow(igrd),
     &                                 nlay(igrd),conc(iptr4d(igrd)))
             enddo
           endif
c
c-----Write PA results
c
           if (lipr) call wrtipr(date,time)
           if (lirr) call wrtirr(date,time)
c
c   --- call routine to zero out all Process Analysis data structures ---
c
           call pazero()
c
c-----Get initial concentration for next loop
c
           if( lipr ) then
              linit = .TRUE.
              do igrd = 1,ngrid
                 call initipr(linit,igrd,nspec,ncol(igrd),nrow(igrd),
     &                                 nlay(igrd),conc(iptr4d(igrd)))
              enddo
           endif
        endif
c
c========================== Process Analysis End ==============================
c
c-----Write model mass
c
        do igrd = 1,ngrid
          call wrtmass(igrd,date,time,1)
        enddo
c
c-----Flush file units
c
        call flush(iavg)
        call flush(ifavg)
        call flush(iconc(1))
        call flush(iconc(2))
        call flush(ifconc(1))
        call flush(ifconc(2))
        call flush(iout)
        call flush(idiag)
        call flush(imass)
c
        whr = aint(wrttim/100.)
        wmn = amod(wrttim,100.)
        wrttim = 100.*(whr + aint((wmn + dtout)/60.)) + 
     &             amod((wmn + dtout),60.)
        if (wrttim.ge.2400.) then
          wrttim = wrttim - 2400.
          wrtdate = wrtdate + 1
          if( MOD(wrtdate,1000) .GT. 365 ) then
              if( MOD(INT(wrtdate/1000),4) .EQ. 0 ) then
                 if( MOD(wrtdate,1000) .EQ. 367 )
     &                     wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
              else
                 wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
              endif
           endif
        endif
      endif
c
c-----Check for end of simulation
c
      tcpu = etime(tarray)
      write(iout,'(a,f10.3)') 'Accumulative CPU time:  ',tcpu 
      write(*,'(a,f10.3)') 'Accumulative CPU time:  ',tcpu
      if (date.lt.enddate) goto 100
      if (date.eq.enddate .and. time.lt.endtim - 0.01) goto 100
c
c------------------  End main time-integration loop  -------------------
c
      call closefl(iavg)
      call closefl(ncpig)
      write(iout,'(/,a,i8.5,f8.2,/)') 'time/Date: ',date,time
      write(*,'(/,a,i8.5,f8.0,/)')'Date/time: ',date,time
      write(iout,'(a)')'END SIMULATION'
      write(*,'(a)')'END SIMULATION'
      write(iout,'(a,i10)') 'TOTAL COARSE GRID TIME STEPS: ', nsteps
      write(*,'(a,i10)') 'TOTAL COARSE GRID TIME STEPS: ', nsteps
c
      end
