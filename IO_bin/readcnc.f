      subroutine readcnc
c 
c     
c-----CAMx v4.02 030709
c
c     READCNC operates in two modes:
c        1) when LRSTRT = F, reads and cycles through the AIRQUALITY file
c           to current time/date, initializes coarse grid concentrations,
c           and maps concentrations to any nested grids
c        2) when LRSTRT = T, reads and cycles through the coarse and
c           fine grid INSTANT files (from previous run) to current
c           time/date, and initializes all grid concentrations directly
c     Finally, the routine reads the TOPCON file and initializes CALOFT
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002
c     ENVIRON International Corporation
c           
c     Modifications: 
c        1/7/99    Added conditional to call rdfgcon for only if there are nests
c      10/31/01    Improved assignment of coarse grid IC's based on whether
c                  this is a restart or not
c      04/17/03    Changed the logic so it no longer such a large local array.
c                  It was causing problems with stack size on the latest PGI
c                  compiler.
c  
c     Input arguments: 
c        none
c             
c     Output arguments: 
c        none
c             
c     Routines Called: 
c        INTRPCNC
c        RDFGCON
c             
c     Called by: 
c        CAMx
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'bndary.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
c
      character*4 icspec(10)
      dimension cinit(MXCOLA,MXROWA)
      character*10 tpspc
c
c-----Entry point
c
      iunit = iic
      if (lrstrt) iunit = irstc
      nx = ncol(1)
      ny = nrow(1)
      nz = nlay(1)
c
c-----Read through coarse grid concentration records until current time/date
c
      read(iunit,end=900) idat1,tim1,idat2,tim2
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      write(iout,'(a40,2(f7.0,i8.5))') 
     &      'Read initial condition file at ',tim1,idat1,tim2,idat2
      do lread = 1,nicspc
        do k = 1,nz
          read(iunit) idum,(icspec(n),n=1,10), 
     &                ((cinit(i,j),i=1,nx),j=1,ny) 
          if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time)) 
     &                 .and. (idat2.gt.date .or. (idat2.eq.date .and. 
     &                                            tim2.gt.time))) then
c
             do 90 lmod = 1,nspec
                lic = licmap(lmod,1)
                if( lic .NE. lread ) goto 90
                do j = 1,ny
                  do i = 1,nx
                    n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                    n4d = n3d + nx*ny*nz*(lmod - 1)
                    if (lrstrt) then
                      conc(n4d) = cinit(i,j)
                    else
                      conc(n4d) = bdnl(lmod)
                      conc(n4d) = amax1(conc(n4d),cinit(i,j))
                    endif
c-----Add if statements to make IC zero----------------------
c                    if (lread.ge.36.and.lread.le.41) then
c                      conc(n4d) = 0.001
c                    endif
c                    if (lread.ge.18.and.lread.le.23) then
c                      conc(n4d) = 0.001
c                    endif
c------------------------------------------------------------

c                    if (i.eq.5.and.j.eq.5.and.k.eq.1) then
c                      write(6,*)lic,lread,conc(n4d)
c                    endif
                  enddo
                enddo
  90         continue
          endif
        enddo
      enddo
c
c-----If this is not a restart, interpolate coarse grid concentrations
c     to all fine grids
c
      if (.not.lrstrt) then
        if (ngrid.gt.1) then
          do ip = 1,ngrid
            do ic = 1,nchdrn(ip)
              ig = idchdrn(ic,ip)
              call intrpcnc(nspec,ncol(ip),nrow(ip),nlay(ip),i1(ig),
     &                      j1(ig),nmesh(ig),nmshv(1,ig),ncol(ig),
     &                      nrow(ig),nlay(ig),conc(iptr4d(ip)),
     &                      conc(iptr4d(ig)) )
            enddo
          enddo
        endif
c
c-----Convert from ppm to umol/m3
c
        do igrd = 1,ngrid
          nx = ncol(igrd)
          ny = nrow(igrd)
          nz = nlay(igrd)
          do l = 1,nspec
            do k = 1,nz
              do j = 1,ny
                do i = 1,nx
                  n3d = i + nx*(j - 1) + nx*ny*(k - 1)
                  n4d = n3d + nx*ny*nz*(l - 1)
                  if (l.le.ngas) then
                    convfac = densfac*273./tempk(iptr3d(igrd)-1+n3d)*
     &                        press(iptr3d(igrd)-1+n3d)/1013.
                  else
                    convfac = 1.
                  endif
cbk                  conc(iptr4d(igrd)-1+n4d) = 
cbk     &                AMAX1( bdnl(l),convfac*conc(iptr4d(igrd)-1+n4d) )
                  conc(iptr4d(igrd)-1+n4d) = convfac*
     &                     AMAX1( bdnl(l), conc(iptr4d(igrd)-1+n4d) )
                enddo
              enddo
            enddo
          enddo
        enddo
      else
c
c-----Otherwise, read fine grid concentrations from fine grid restart file
c
        if( ngrid .GT. 1 ) call rdfgcon(idat1,tim1)
      endif
c
c-----Read TOPCON file and initialize top concentrations
c
      do l = 1,nspec
        caloft(l) = 0.
      enddo
c
      write(idiag,*)
 200  read(itopc,'(a10,f10.0)',end=901) tpspc,ctin
        if (tpspc.eq.'HNO2      ') tpspc = 'HONO      '
        if (tpspc.eq.'HCHO      ' .and. kHCHO.eq.nspec+1)
     &                                        tpspc = 'FORM      '
      do l = 1,nspec
        if (tpspc.eq.spname(l)) then
          write(idiag,'(a,a10,a)') 
     &          'Topcon species ',tpspc,' read from file'
          caloft(l) = ctin
          goto 200
        endif
      enddo
      write(idiag,'(a,a10)') 'Topcon species not on internal list',
     &                        tpspc
      goto 200
c
c-----Use BDNL value for species not on CALOFT file
c
 901  do l = 1,nspec
        caloft(l) = amax1(caloft(l), bdnl(l))
      enddo
      goto 999
c
c-----End of IC file reached
c
 900  write(iout,'(//,a)') 'ERROR in READCNC:'
      write(iout,*)'End of IC file'
      write(iout,*)'Make sure initial condition file contains the ',
     &                                   'simulation beginning hour.'
      call camxerr()
c
 999  return
      end
