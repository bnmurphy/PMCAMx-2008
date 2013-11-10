      subroutine readinp(igrd,ngcol,ngrow,nglay,hght,phptim,hnext,
     &                   presure,ppptim,pnext,wndu,puptim,unext,wndv,
     &                   pvptim,vnext,tsrf,psptim,tsnext,temper,ptptim,
     &                   tnext,vapor,cloud,rkvgrd,cldwtr,pcpwtr,
     &                   cldod,energy,rdum)
c 
c-----CAMx v4.03 031205
c 
c     READINP cycles through and reads all time-variant environmental 
c     input files.  The following files are read:
c     - Height/pressure fields at next update date/time, all grids (optional)
c     - Wind fields            at next update date/time, all grids (optional)
c     - Temperature fields     at next update date/time, all grids (optional)
c     - Water vapor fields     at current date/time,     all grids (optional)
c     - Cloud/rain fields      at current date/time,     all grids (optional)
c     - Kv fields              at current date/time,     all grids (optional)
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        1/13/99   Added checks to convert input hours at 2400 to 0 on the
c                  following day
c       11/06/01   CAMx now assumes that all input file dates are in
c                  Julian format (YYJJJ) if the simulation year is 2000
c                  or greater.  Added call to DATEERR for time mismatches.
c        8/30/02   Cloud and rain files combined into one file (CAMx formatted
c                  file only), and now cloud/rain and water vapor files can be 
c                  read in for each nest
c        1/23/03   Cloud file contains new parameters:
c                  -cloud water content (cwc)
c                  -precip water content (pwc)
c                  -optical depth (cod)
c        4/2/03    Removed option for UAM-V type cloud file
c 
c     Input arguments: 
c        igrd                grid index
c        ngcol               number of columns
c        ngrow               number of rows
c        nglay               number of layers
c        hght                current layer interface field (m)
c        presure             current pressure field (mb)
c        wndu                current u-component wind field (m/s)
c        wndv                current v-component wind field (m/s)
c        tsrf                current surface temperature field (K)
c        temper              current temperature field (K)
c 
c     Output arguments: 
c        phptim              time-rate change of layer interface height (m/s)
c        hnext                next layer interface height field (m)
c        ppptim              time-rate change of pressure (mb/s)
c        pnext                next pressure field (mb)
c        puptim              time-rate change of u-component wind (m/s2)
c        unext                next u-component wind field (m/s)
c        pvptim              time-rate change of v-component wind (m/s2)
c        vnext                next v-component wind field (m/s)
c        psptim              time-rate change of surface temperature (K/s)
c        tsnext               next surface temperature field (m/s)
c        ptptim              time-rate change of 3-D temperature (K/s)
c        tnext                next 3-D temperature field (m/s)
c        vapor               water vapor field (ppm)
c        cloud               fractional cloud cover field (decimal fraction)
c        rkvgrd              vertical exchange coefficient (m2/s)
c        cldwtr              cloud water content (g/m3)
c        pcpwtr              precipitation water content (g/m3)
c        cldod               cloud optical depth
c        energy              energy transmission coefficient
c 
c     Routines Called: 
c        JULDATE
c        TIMRATES
c        CVTWIND
c        ZEROS
c        RASSGN3D
c 
c     Called by: 
c        CAMx 
c 
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'flags.com'
      include 'grid.com'
      include 'camxfld.com'
c
      integer hdate
      dimension hnext(ngcol,ngrow,nglay),pnext(ngcol,ngrow,nglay),
     &          unext(ngcol,ngrow,nglay),vnext(ngcol,ngrow,nglay),
     &          tnext(ngcol,ngrow,nglay)
      dimension hght(ngcol,ngrow,nglay),phptim(ngcol,ngrow,nglay),
     &          presure(ngcol,ngrow,nglay),ppptim(ngcol,ngrow,nglay),
     &          wndu(ngcol,ngrow,nglay),puptim(ngcol,ngrow,nglay),
     &          wndv(ngcol,ngrow,nglay),pvptim(ngcol,ngrow,nglay),
     &          temper(ngcol,ngrow,nglay),ptptim(ngcol,ngrow,nglay),
     &          vapor(ngcol,ngrow,nglay),rkvgrd(ngcol,ngrow,nglay)
      real cloud(ngcol,ngrow,nglay),cldwtr(ngcol,ngrow,nglay),
     &     pcpwtr(ngcol,ngrow,nglay),cldod(ngcol,ngrow,nglay),
     &     energy(ngcol,ngrow,nglay)
      dimension tsrf(ngcol,ngrow),psptim(ngcol,ngrow),
     &          tsnext(ngcol,ngrow)
      character*60 string
c
c-----Entry point
c
      whr = aint(time/100.)
      wmn = amod(time,100.)
      htim = 100.*(whr + aint((wmn + dtinp)/60.)) +
     &             amod((wmn + dtinp),60.)
      hdate = date 
      if (htim.ge.2400.) then
        htim = htim - 2400.
        hdate = hdate + 1
      endif
      if( MOD(hdate,1000) .GT. 365 ) then
         if( MOD(INT(hdate/1000),4) .EQ. 0 ) then
            if( MOD(hdate,1000) .EQ. 367 )
     &                     hdate = (INT(hdate/1000)+1)*1000 + 1
         else
            hdate = (INT(hdate/1000)+1)*1000 + 1
         endif
      endif
c
c-----Read height/pressure file for coarse grid and optionally for any fine
c     grids
c
      iunit = ihtp(igrd)
      if( iunit .GT. 0 ) then
 101     continue
         do k = 1,nglay
           read(iunit,end=7000) hr,idt,
     &                        ((hnext(i,j,k),i=1,ngcol),j=1,ngrow) 
           read(iunit,end=7000) hr,idt,
     &                        ((pnext(i,j,k),i=1,ngcol),j=1,ngrow) 
         enddo   
c
         if (.not.ly2k .and. idt.gt.100000) call juldate(idt)
         if (hr.ge.2400.) then
           hr = hr - 2400.
           idt = idt + 1
           if( MOD(idt,1000) .GT. 365 ) then
              if( MOD(INT(idt/1000),4) .EQ. 0 ) then
                 if( MOD(idt,1000) .EQ. 367 )
     &                     idt = (INT(idt/1000)+1)*1000 + 1
              else
                 idt = (INT(idt/1000)+1)*1000 + 1
              endif
           endif
         endif
         write(iout,'(a40,f7.0,i8.5,a,i3)') 
     &         'Read height/pressure file at ',hr,idt,' grid',igrd 
c             
         if (idt.lt.hdate .or. (idt.eq.hdate .and. hr.lt.htim)) goto 101 
         if (idt.gt.hdate .or. (idt.eq.hdate .and. hr.gt.htim)) then 
           write(iout,'(//,a)') 'ERROR in READINP:'
           write(iout,'(a,f10.1,i10.5)')
     &           'Past current time/date ',time,date 
           if( igrd .EQ. 1 ) then
              string = 'Reading height/pressure file for coarse grid.'
           else
              write(string,'(A,I5)')
     &                   'Reading height/pressure file for grid:',igrd
           endif
           call dateerr(ly2k,string)
         endif 
c
c-----Interpolate/assign the values from the parent, if no data available
c
      else
         do ip=1,ngrid
           do ic = 1,nchdrn(ip)
             if( igrd .EQ. idchdrn(ic,ip) ) then
               write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Assigning heights from parent grid',
     &                             time, date,' grid',igrd
               call rassgn3d(ncol(ip),nrow(ip),nlay(ip),
     &           i1(igrd),j1(igrd),nmesh(igrd),ncol(igrd),nrow(igrd),
     &                            hnxt(iptr3d(ip)),hnxt(iptr3d(igrd)) )
               write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating pressure from parent grid',
     &                             time, date,' grid',igrd
               call interp2d(ncol(ip),nrow(ip),nlay(ip),
     &             i1(igrd),j1(igrd),nmesh(igrd),ncol(igrd),nrow(igrd),
     &                            pnxt(iptr3d(ip)),pnxt(iptr3d(igrd)) )
               call interpv(ncol(ip),nrow(ip),nlay(ip),ncol(igrd),
     &                   nrow(igrd),nlay(igrd),nmesh(igrd),
     &                   nmshv(1,igrd),i1(igrd),j1(igrd),
     &                   hnxt(iptr3d(ip)),
     &                   hnxt(iptr3d(igrd)),pnxt(iptr3d(igrd)) )
             endif
           enddo
         enddo
      endif
c
c-----Calculate time rates of change for layer interface heights and pressure
c
      call timrates(ngcol,ngrow,nglay,hght,hnext,phptim)
      call timrates(ngcol,ngrow,nglay,presure,pnext,ppptim)
c
c-----Read wind file for coarse grid and optionally for any fine grids
c
      iunit = iwind(igrd)
      if (iunit.gt.0) then
 201    continue
        read(iunit,end=7000) hr,idt
        do k = 1,nglay
          read(iunit,end=7000) ((unext(i,j,k),i=1,ngcol),j=1,ngrow)
          read(iunit,end=7000) ((vnext(i,j,k),i=1,ngcol),j=1,ngrow)
        enddo
        read(iunit,end=7000) dum
c
        if (.not.ly2k .and. idt.gt.100000) call juldate(idt)
        if (hr.ge.2400.) then
          hr = hr - 2400.
          idt = idt + 1
          if( MOD(idt,1000) .GT. 365 ) then
             if( MOD(INT(idt/1000),4) .EQ. 0 ) then
                if( MOD(idt,1000) .EQ. 367 )
     &                     idt = (INT(idt/1000)+1)*1000 + 1
             else
                idt = (INT(idt/1000)+1)*1000 + 1
             endif
          endif
        endif
        write(iout,'(a40,f7.0,i8.5,a,i3)') 
     &        'Read wind file at ',hr,idt,' grid',igrd
c             
        if (idt.lt.hdate .or. (idt.eq.hdate .and. hr.lt.htim)) goto 201 
        if (idt.gt.hdate .or. (idt.eq.hdate .and. hr.gt.htim)) then 
          write(iout,'(//,a)') 'ERROR in READINP:'
          write(iout,'(a,f10.1,i10.5)')
     &          'Past current time/date ',time,date 
           if( igrd .EQ. 1 ) then
              string = 'Reading wind file for coarse grid.'
           else
              write(string,'(A,I5)') 'Reading wind file for grid:',igrd
           endif
          call dateerr(ly2k,string)
        endif 
c
c-----Convert wind vectors if they are not staggered
c
        if (.not.lstagw) call cvtwind(ngcol,ngrow,nglay,unext,vnext)
c
c-----Calculate time rates of change for u and v winds on coarse grid
c
        if (igrd.eq.1) then
          call timrates(ngcol,ngrow,nglay,wndu,unext,puptim)
          call timrates(ngcol,ngrow,nglay,wndv,vnext,pvptim)
        endif
      endif
c
c-----Read temperature file for coarse grid and optionally any fine grids
c
      iunit = itemp(igrd)
      if (iunit.gt.0) then
 301    read(iunit,end=7000) hr,idt,
     &                   ((tsnext(i,j),i=1,ngcol),j=1,ngrow)
        do k = 1,nglay
          read(iunit,end=7000) hr,idt,
     &                 ((tnext(i,j,k),i=1,ngcol),j=1,ngrow)
        enddo
c
        if (.not.ly2k .and. idt.gt.100000) call juldate(idt)
        if (hr.ge.2400.) then
          hr = hr - 2400.
          idt = idt + 1
          if( MOD(idt,1000) .GT. 365 ) then
             if( MOD(INT(idt/1000),4) .EQ. 0 ) then
                if( MOD(idt,1000) .EQ. 367 )
     &                     idt = (INT(idt/1000)+1)*1000 + 1
             else
                idt = (INT(idt/1000)+1)*1000 + 1
             endif
          endif
        endif
        write(iout,'(a40,f7.0,i8.5,a,i3)')
     &        'Read temperature file at ',hr,idt,' grid',igrd
c             
        if (idt.lt.hdate .or. (idt.eq.hdate .and. hr.lt.htim)) goto 301
        if (idt.gt.hdate .or. (idt.eq.hdate .and. hr.gt.htim)) then 
          write(iout,'(//,a)') 'ERROR in READINP:'
          write(iout,'(a,f10.1,i10.5)')
     &          'Past current time/date ',time,date 
           if( igrd .EQ. 1 ) then
              string = 'Reading temperature file for coarse grid.'
           else
              write(string,'(A,I5)') 
     &                     'Reading temperature file for grid:',igrd
           endif
          call dateerr(ly2k,string)
        endif 
c
c-----Calculate time rates of change for surface and 3-D temperature
c
        call timrates(ngcol,ngrow,1,tsrf,tsnext,psptim)
        call timrates(ngcol,ngrow,nglay,temper,tnext,ptptim)
      endif
c
c-----Read water vapor file for coarse grid and optionally any fine grids
c
      if (lchem) then
        iunit = ih2o(igrd)
        if (iunit.gt.0) then
 401      do k = 1,nglay
            read(iunit,end=7000) hr,idt,
     &               ((vapor(i,j,k),i=1,ngcol),j=1,ngrow)
          enddo
c
          if (.not.ly2k .and. idt.gt.100000) call juldate(idt)
          if (hr.ge.2400.) then
            hr = hr - 2400.
            idt = idt + 1
            if( MOD(idt,1000) .GT. 365 ) then
               if( MOD(INT(idt/1000),4) .EQ. 0 ) then
                  if( MOD(idt,1000) .EQ. 367 )
     &                     idt = (INT(idt/1000)+1)*1000 + 1
               else
                  idt = (INT(idt/1000)+1)*1000 + 1
               endif
            endif
          endif
          write(iout,'(a40,f7.0,i8.5,a,i3)')
     &          'Read water vapor file at ',hr,idt,' grid',igrd
c
          if (idt.lt.date .or. (idt.eq.date .and. hr.lt.time)) goto 401 
          if (idt.gt.date .or. (idt.eq.date .and. hr.gt.time)) then 
            write(iout,'(//,a)') 'ERROR in READINP:'
            write(iout,'(a,f10.1,i10.5)')
     &            'Past current time/date ',time,date 
           if( igrd .EQ. 1 ) then
              string = 'Reading water vapor file for coarse grid.'
           else
              write(string,'(A,I5)') 
     &                  'Reading water vapor file for grid:',igrd
           endif
            call dateerr(ly2k,string)
          endif 
        endif
      endif
c
c-----Read cloud file for coarse grid and optionally any fine grids
c
      if (lchem .or. lwet .or. ldry) then
        iunit = icld(igrd)
        if (igrd.eq.1 .and. icld(igrd).eq.0) then
          call zeros(cloud,ngcol*ngrow*nglay)
          call zeros(cldwtr,ngcol*ngrow*nglay)
          call zeros(pcpwtr,ngcol*ngrow*nglay)
          call zeros(cldod, ngcol*ngrow*nglay)
          call zeros(energy,ngcol*ngrow*nglay)
          if (lwet) then
            write(iout,'(//,a)')'ERROR in READINP:'
            write(iout,*)'Wet Deposition is selected but ',
     &                   'cloud/rain file is not given'
            call camxerr()
          endif
        elseif (iunit.gt.0) then
c
c-----Read CAMx cloud/rain file
c
 501      read(iunit,end=7000,err=7001) hr,idt
          do k = 1,nglay
            read(iunit,end=7000,err=7001)
     &           ((cldwtr(i,j,k),i=1,ngcol),j=1,ngrow)
            read(iunit,end=7000,err=7001)
     &           ((pcpwtr(i,j,k),i=1,ngcol),j=1,ngrow) 
            read(iunit,end=7000,err=7001)
     &           ((cldod(i,j,k),i=1,ngcol),j=1,ngrow) 
          enddo 
c
          if (.not.ly2k .and. idt.gt.100000) call juldate(idt)
          if (hr.ge.2400.) then
            hr = hr - 2400.
            idt = idt + 1
            if( MOD(idt,1000) .GT. 365 ) then
               if( MOD(INT(idt/1000),4) .EQ. 0 ) then
                  if( MOD(idt,1000) .EQ. 367 )
     &               idt = (INT(idt/1000)+1)*1000 + 1
               else
                  idt = (INT(idt/1000)+1)*1000 + 1
               endif
            endif
          endif
          write(iout,'(a40,f7.0,i8.5,a,i3)')
     &          'Read cloud/rain file at ',hr,idt,' grid',igrd
c               
          if (idt.lt.date .or.(idt.eq.date .and. hr.lt.time)) goto 501
          if (idt.gt.date .or. (idt.eq.date .and. hr.gt.time)) then 
            write(iout,'(//,a)') 'ERROR in READINP:'
            write(iout,'(a,f10.1,i10.5)')
     &            'Past current time/date ',time,date 
            if( igrd .EQ. 1 ) then
              string = 'Reading cloud/rain file for coarse grid.'
            else
              write(string,'(A,I5)') 
     &                 'Reading cloud/rain file for grid:',igrd
            endif
            call dateerr(ly2k,string)
          endif 
c
c-----Calculate energy transmission coefficient if using RADM cloud adjustment
c
          do k = 1,nglay
            do j = 1,ngrow
              do i = 1,ngcol
                tau = cldod(i,j,k)
                if (tau.lt.5.) then
                  energy(i,j,k) = 1.
                  cloud(i,j,k) = 0.
                else
                  energy(i,j,k) = (5. - exp(-tau))/(4. + 0.42*tau)
                  cloud(i,j,k) = 1.
                endif
              enddo
            enddo
          enddo

        endif
      endif
c
c-----Read vertical diffusion coefficient for coarse grid and optionally
c     for any fine grids
c
      iunit = ikv(igrd)
      if (iunit.gt.0) then
 701    do k = 1,nglay
          read(iunit,end=7000,err=7002) hr,idt,
     &               ((rkvgrd(i,j,k),i=1,ngcol),j=1,ngrow)
        enddo
c
        if (.not.ly2k .and. idt.gt.100000) call juldate(idt)
        if (hr.ge.2400.) then
          hr = hr - 2400.
          idt = idt + 1
          if( MOD(idt,1000) .GT. 365 ) then
             if( MOD(INT(idt/1000),4) .EQ. 0 ) then
                if( MOD(idt,1000) .EQ. 367 )
     &                     idt = (INT(idt/1000)+1)*1000 + 1
             else
                idt = (INT(idt/1000)+1)*1000 + 1
             endif
          endif
        endif
        write(iout,'(a40,f7.0,i8.5,a,i3)')
     &        'Read KV file at ',hr,idt,' grid',igrd
c             
        if (idt.lt.date .or. (idt.eq.date .and. hr.lt.time)) goto 701 
        if (idt.gt.date .or. (idt.eq.date .and. hr.gt.time)) then 
          write(iout,'(//,a)') 'ERROR in READINP:'
          write(iout,'(a,f10.1,i10.5)')
     &          'Past current time/date ',time,date 
           if( igrd .EQ. 1 ) then
              string = 'Reading KV file for coarse grid.'
           else
              write(string,'(A,I5)')'Reading KV file for grid:',igrd
           endif
          call dateerr(ly2k,string)
        endif 
      endif
      goto 9999
c
 7000 continue
      write(iout,'(//,a)')'ERROR in READINP:'
      write(iout,*)'End of input file reached.  Make sure the file '
      write(iout,*)'is for the correct day and contains all hours.'
      write(iout,*)
      if (ly2k) then
        write(iout,*)'You are modeling year 2000 or later --'
        write(iout,*)'CAMx assumes that all input file dates',
     &               ' are in Julian format (YYJJJ)'
      endif
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)')'ERROR in READINP:'
      write(iout,*)'There was an error reading a record of the'
      write(iout,*)'cloud/rain file.  This file has a new format'
      write(iout,*)'so make sure it contains a rain record.'
      write(iout,*)'Also make sure you have provided cloud/rain file'
      write(iout,*)'records in CAMx.in for all fine grids;'
      write(iout,*)'if you do not have fine grid cloud/rain files,'
      write(iout,*)'the respective CAMx.in records should be blank.'
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)')'ERROR in READINP:'
      write(iout,*)'There was an error reading a record of the'
      write(iout,*)'diffusivity file -- it may be trying to read'
      write(iout,*)'an old rain file. Check to make sure you do'
      write(iout,*)'not have an old rain file listed in CAMx.in;'
      write(iout,*)'rain info is now contained in the new cloud/rain'
      write(iout,*)'file. See the list of I/O files opened by CAMx in'
      write(iout,*)'the .out file'
      call camxerr()

 9999 continue
      return
      end
