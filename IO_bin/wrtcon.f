      subroutine wrtcon(iflag,tim2,idat2,iunit,nox,noy,noz,
     &                  nsptmp,cncfld)
c
c-----CAMx v4.03 031205
c 
c     WRTCON writes average and instantaneous concentration fields; for
c     average files, optionally writes only layer 1; for instantaneous files,
c     writes all layers.
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        1/20/99   Grid cell size on file should be meters for all cartesian
c                  projections (UTM, LCP, PSP)
c        10/24/01  Removed BSWAP and converted integer strings to character*4
c 
c     Input arguments:
c        iflag               output type flag (0=average/1=instantaneous)
c        tim2                output time (HHMM)
c        idat2               output date (YYJJJ)
c        iunit               output unit
c        nox                 number of cells in x-direction
c        noy                 number of cells in y-direction
c        noz                 number of layers
c        nsptmp              number of species
c        cncfld              concentration field to output (ppm or umol/m3)
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        CAMx
c 
      include 'camx.prm'
      include 'camx.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
c
      character*4 ispec(10,MXSPEC), ifile(10), note(60)
      real cncfld(nox,noy,noz,nsptmp)
      character*10 cnfil
c
c-----Data statements
c
      data cnfil /'INSTANT   '/
      data nseg,izero /1,0/
      data zero /0./
c
c-----Entry point
c
c-----Determine time/date range
c
      idat1 = idat2 
      etim = ANINT(tim2)/100.
      if( dtout .GE. 60.0 ) then
          btim = ANINT( 1000*(etim - ANINT(dtout)/60.) )/1000.
      else
          btim = ANINT( 1000*(etim - ANINT(dtout)/100.) )/1000.
      endif
      if (btim.lt.0.) then 
        btim = btim + 24. 
        idat1 = idat1 - 1 
      endif 
      idat3 = idat2
      etim3 = etim + 0.1
      if (etim3.gt.24.) then 
        etim3 = 24. - etim3 
        idat3 = idat3 + 1 
        if( MOD(idat3,1000) .GT. 365 ) then
            if( MOD(INT(idat3/1000),4) .EQ. 0 ) then
               if( MOD(idat3,1000) .EQ. 367 )
     &                     idat3 = (INT(idat3/1000)+1)*1000 + 1
            else
               idat3 = (INT(idat3/1000)+1)*1000 + 1
            endif
         endif
      endif 
c
c-----For instantaneous files, rewind file and write the header
c
      nlayer = noz
      if (iflag.eq.1) then
        read(cnfil,'(10a1)') (ifile(n),n=1,10)
        read(runmsg(1:60),'(60a1)') (note(n),n=1,60)
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
        do l=1,nsptmp
          read(spname(l),'(10a1)') (ispec(n,l),n=1,10) 
        enddo
c
        rewind(iunit)
        write(iunit) ifile,note,nseg,nsptmp,idat2,etim,idat3,etim3
        write(iunit) zero,zero,izone,orgx,orgy,dx,dy,
     &               nox,noy,noz,izero,izero,zero,zero,zero
        write(iunit) izero,izero,nox,noy
        write(iunit) ((ispec(n,l),n=1,10),l=1,nsptmp)
        write(iunit) idat2,etim,idat3,etim3
      else
        do l = 1,nsptmp
          read(spname(lavmap(l)),'(10a1)') (ispec(n,l),n=1,10)
        enddo
        if (.not.l3davg) nlayer = 1 
        write(iunit) idat1,btim,idat2,etim
      endif
c
c-----Write gridded concentration field
c
      do l = 1,nsptmp
        do k = 1,nlayer
          write(iunit) nseg,(ispec(n,l),n=1,10),
     &                 ((cncfld(i,j,k,l),i=1,nox),j=1,noy)
        enddo
      enddo
c
      return
      end
