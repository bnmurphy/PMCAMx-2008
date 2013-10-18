      subroutine parntchd(ngrid,ncol,nrow,inst1,inst2,jnst1,jnst2,
     &                    meshold,mapgrd,i1,i2,j1,j2,nmesh,nchdrn,
     &                    idchdrn)
c
c-----CAMx v4.02 030709
c
c     PARNTCHD figures out the parent-children relationship and relative
c     mesh numbers between the grids
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        ngrid               number of grids
c        ncol                number of columns in each grid
c        nrow                number of rows in each grid
c        inst1/inst2         starting/ending i-indices relative to CG     
c        jnst1/jnst2         starting/ending j-indices relative to CG
c        meshold             meshing factor relative to CG
c
c     Output arguments:
c        mapgrd              grid mapping indices
c        i1/i2               starting/ending i-indices relative to parent 
c        j1/j2               starting/ending j-indices relative to parent
c        nmesh               meshing factor relative to parent
c        nchdrn              number of children per parent grid
c        idchdrn             indices of children per parent grid
c
c     Routines called:
c        none
c
c     Called by:
c        NSTPREP
c
      include "camx.prm"
      include "filunit.com"
c
      dimension ncol(ngrid),nrow(ngrid),inst1(ngrid),inst2(ngrid),
     &          jnst1(ngrid),jnst2(ngrid),meshold(ngrid),
     &          i1(ngrid),i2(ngrid),j1(ngrid),j2(ngrid),
     &          nmesh(ngrid),nchdrn(ngrid),idchdrn(MXCHDRN,MXGRID)
      dimension igen(MXGRID),mapgrd(ngrid),
     &          itmp1(MXGRID),itmp2(MXGRID),jtmp1(MXGRID),jtmp2(MXGRID)
c
c-----Entry point
c
c-----Check for overlap
c
      do l1 = 2,ngrid
        do l2 = 2,ngrid
          if (l1.eq.l2) then
            if (inst1(l2).ge.inst1(l1) .and. 
     &          jnst1(l2).ge.jnst1(l1) .and. 
     &          inst1(l2).le.inst2(l1) .and. 
     &          jnst1(l2).le.jnst2(l1)) then
              if (inst2(l2).gt.inst2(l1) .or. 
     &            jnst2(l2).gt.jnst2(l1)) then
                write(iout,'(//,a)') 'ERROR in PARNTCHD:'
                write(iout,*) 'Two grids overlap: ',l1,l2
                call camxerr()
              endif
            endif
            if (inst2(l2).ge.inst1(l1) .and. 
     &          jnst1(l2).ge.jnst1(l1) .and. 
     &          inst2(l2).le.inst2(l1) .and. 
     &          jnst1(l2).le.jnst2(l1)) then
              if (inst1(l2).gt.inst1(l1) .or. 
     &            jnst2(l2).gt.jnst2(l1)) then
                write(iout,'(//,a)') 'ERROR in PARNTCHD:'
                write(iout,*) 'Two grids overlap: ',l1,l2
                call camxerr()
              endif
            endif
          endif
        enddo
      enddo
c
c-----Set up grids
c
      igen(1) = 1
      do igrd = 2,ngrid
        igen(igrd) = 2
      enddo
c
c-----Determine the generations by sorting
c
      do iter = 1,ngrid*ngrid
        indx = 0
        do l1 = 2,ngrid
          do l2 = 2,ngrid
            if (l1.ne.l2) then
              if (inst1(l2).ge.inst1(l1) .and. jnst1(l2).ge.jnst1(l1) 
     &                                   .and. 
     &            inst2(l2).le.inst2(l1) .and. jnst2(l2).le.jnst2(l1))
     &                                                              then
                if (igen(l2).le.igen(l1)) then
                  igen(l2) = igen(l1) + 1
                  indx = 1
                endif
              endif
            endif
          enddo
        enddo
        if (indx.eq.0) goto 50
      enddo
  50  continue
c
c-----Assign grid ID
c
      do itmp = 1,ngrid
        mapgrd(itmp) = itmp 
      enddo
c
      do iter = 1,ngrid*ngrid
        indx = 0
        do igrd = 2,ngrid-1
          if (igen(igrd).gt.igen(igrd+1)) then
            igntmp = igen(igrd)
            igen(igrd) = igen(igrd+1)
            igen(igrd+1) = igntmp
            igntmp = mapgrd(igrd)
            mapgrd(igrd) = mapgrd(igrd+1)
            mapgrd(igrd+1) = igntmp
            indx = 1
          endif
        enddo
        if (indx.eq.0) goto 60
      enddo
  60  continue
c
c-----Assign initial relative cell index range
c
      do igrd = 2,ngrid
        i1(igrd) = inst1(mapgrd(igrd)) 
        j1(igrd) = jnst1(mapgrd(igrd)) 
        i2(igrd) = inst2(mapgrd(igrd)) 
        j2(igrd) = jnst2(mapgrd(igrd)) 
      enddo
c
c-----Determine number of children and associated grid indices
c
      i1(1) = 1
      j1(1) = 1
      i2(1) = ncol(1)
      j2(1) = nrow(1)
      do ip = 1,ngrid
        nchdrn(ip) = 0
        do igrd = 2,ngrid
          if (igen(igrd).eq.igen(ip)+1) then
            if (i1(igrd).ge.i1(ip) .and. j1(igrd).ge.j1(ip) .and. 
     &          i2(igrd).le.i2(ip) .and. j2(igrd).le.j2(ip)) then
              nchdrn(ip) = nchdrn(ip) + 1
              idchdrn(nchdrn(ip),ip) = igrd
            endif
          endif
        enddo
      enddo
c
c-----Calculate relative cell index range
c
      mapgrd(1) = 1
      meshold(1) = 1
      nmesh(1) = 1
      i1(1) = 2
      j1(1) = 2 
      do ip=1,ngrid
        meshp = meshold(mapgrd(ip))
        do ic=1,nchdrn(ip)
          igrd = idchdrn(ic,ip)
          meshc= meshold(mapgrd(igrd))
          if (mod(meshc,meshp).ne.0) then
            write(iout,'(//,a)') 'ERROR in PARNTCHD:'
            write(iout,*) 'Inappropriate mesh number for grid ',
     &                     mapgrd(igrd), ' mesh # = ', meshc 
            call camxerr()
          endif
          nmesh(igrd) = meshold(mapgrd(igrd))/meshold(mapgrd(ip))
          itmp1(igrd) = (i1(igrd) - i1(ip))*meshold(mapgrd(ip)) + 2
          jtmp1(igrd) = (j1(igrd) - j1(ip))*meshold(mapgrd(ip)) + 2
          itmp2(igrd) = (i2(igrd) - i1(ip) + 1)*meshold(mapgrd(ip)) + 1
          jtmp2(igrd) = (j2(igrd) - j1(ip) + 1)*meshold(mapgrd(ip)) + 1
        enddo
      enddo
c
      i1(1) = 0
      j1(1) = 0 
      do igrd = 2,ngrid
        i1(igrd) = itmp1(igrd)
        j1(igrd) = jtmp1(igrd)
        i2(igrd) = itmp2(igrd)
        j2(igrd) = jtmp2(igrd)
      enddo
c
      return
      end
