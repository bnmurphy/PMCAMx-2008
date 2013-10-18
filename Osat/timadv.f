      subroutine timadv(igrd,xyordr,ncol,nrow,nlay,nspc,deltat,dx,dy,
     &                   windu,windv,depth,
     &                   mapscl,saconc)
c
c-----CAMx v4.02 030709
c
c     TIMADV drives 2-D advection of timing tracers
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        02/05/03    --gwilson--  Removed SMOLAR advections solver.
c        none
c
c     Input arguments:
c        igrd              grid index
c        xyordr            order of x & y advection
c        ncol              number of columns
c        nrow              number of rows
c        nlay              number of layers
c        nspc              number of species
c        deltat            time step (s)
c        dx                cell size in x-direction (m)
c        dy                cell size in y-direction (m)
c        windu             wind speed in x-direction (m/s)
c        windv             wind speed in y-direction (m/s)
c        depth             layer depth (m)
c        mapscl            map-scale factor at cell centroids
c        saconc            species concentrations (umol/m3)
c
c     Output arguments:
c        saconc             species concentrations (umol/m3)
c
c     Routines Called:
c        HADVSMO
c        HADVBOT
c
c     Called by:
c        EMISTRNS
c
      include "camx.prm"
      include "bndary.com"
      include "chmstry.com"
      include "filunit.com"
      include "flags.com"
      include "tracer.com"
c
      integer xyordr
      dimension saconc(ncol,nrow,nlay,nspc)
      real windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &     depth(ncol,nrow,nlay),
     &     mapscl(ncol,nrow),
     &     dx(nrow)
      real c1d(MX1D),v1d(MX1D),flxarr(MX1D),m1d(MX1D),
     &     saflux(MX1D),fpc(MX1D),fmc(MX1D)
      real*8 flux1,flux2
      dimension tarray2(2)
c
c======================== Process Analysis Begin ====================================
c
      real fc1(MX1D), fc2(MX1D)
c
c========================= Process Analysis End =====================================
c
c
c-----Entry point
c
      if (xyordr.eq.0) goto 200
c
c-----Advection in x-direction
c
 100  write(*,'(a20,$)') '  SA x advect ......'  
      write(iout,'(a20,$)') '  SA x advect ......'  
      do 20 k = 1,nlay
        do 10 j = 2,nrow-1
          do 21 ispc = ipttim,nspc
            i1 = 1
            i2 = ncol
            if (igrd.eq.1) then
              if (ibeg(j).eq.-999) goto 10
              i1 = ibeg(j) - 1
              i2 = iend(j) + 1
            endif
c
            l = 0
            do i = i1,i2
              l = l + 1
              v1d(l) = windu(i,j,k)
              if (i.lt.i2)
     &                 v1d(l) = 2.*v1d(l)/(mapscl(i+1,j) + mapscl(i,j))
              c1d(l) = saconc(i,j,k,ispc)*dy*depth(i,j,k)
              m1d(l) = mapscl(i,j)*mapscl(i,j)
            enddo
            nn = i2 - i1 + 1
            if (iadvct.eq.2) then
              call hadvbot(nn,deltat,dx(j),c1d,v1d,m1d,flxarr,flux1,
     &                                flux2,saflux,fpc,fmc,fc1,fc2)
            elseif (iadvct.eq.3) then
              call hadvppm(nn,deltat,dx(j),c1d,v1d,m1d,flxarr,flux1,
     &                                       flux2,saflux,fc1,fc2)
            endif
c
            l = 1
            do i = i1+1,i2-1
              l = l + 1
              saconc(i,j,k,ispc) = c1d(l)/dy/depth(i,j,k)
            enddo
  21      continue
c
  10    continue
  20  continue
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1) 
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
      if (xyordr.eq.0) goto 300
c
c-----Advection in y-direction
c
 200  write(*,'(a20,$)') '  SA y advect ......'  
      write(iout,'(a20,$)') '  SA y advect ......'  
      do 40 k = 1,nlay
        do 30 i = 2,ncol-1
         do 41 ispc = 1,nspc
            j1 = 1
            j2 = nrow
            if (igrd.eq.1) then
              if (jbeg(i).eq.-999) goto 30
              j1 = jbeg(i) - 1
              j2 = jend(i) + 1
            endif
c
            l = 0
            do j = j1,j2
              l = l + 1
              v1d(l) = windv(i,j,k)
              if (j.lt.j2)
     &          v1d(l) = 2.*v1d(l)/(mapscl(i,j+1) + mapscl(i,j))
              c1d(l) = saconc(i,j,k,ispc)*dx(j)*depth(i,j,k)
              m1d(l) = mapscl(i,j)*mapscl(i,j)
            enddo
            nn = j2 - j1 + 1
c
            if (iadvct.eq.2) then
              call hadvbot(nn,deltat,dy,c1d,v1d,m1d,flxarr,flux1,
     &                                    flux2,saflux,fpc,fmc,fc1,fc2)
            elseif (iadvct.eq.3) then
              call hadvppm(nn,deltat,dy,c1d,v1d,m1d,flxarr,flux1,
     &                                           flux2,saflux,fc1,fc2)
            endif
c
            l = 1
            do j = j1+1,j2-1
              l = l+1
              saconc(i,j,k,ispc) = c1d(l)/dx(j)/depth(i,j,k)
            enddo
  41      continue
c
  30    continue
  40  continue
      tcpu = dtime(tarray2) 
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1) 
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
      if (xyordr.eq.0) goto 100
c
 300  continue
      return
      end
