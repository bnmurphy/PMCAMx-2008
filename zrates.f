      subroutine zrates(igrid,xyordr,ncol,nrow,nlay,nadv,deltat,dx,dy,
     &                  depth,phpt,pppt,ptpt,windu,windv,tempk,press,
     &                  mapscl,dilut,entrn)
c 
c-----CAMx v4.02 030709
c  
c     ZRATES calculates new vertical velocity and dilution/entrainment
c     rates resulting from the time-varying vertical grid.
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c 
c     Modifications:
c        4/26/99   Added Piecewise Parabolic Method for horizontal advection
c        12/3/99   A dummy atmospheric density is added above the top of
c                  the model based upon an extrapolation to an overlying
c                  layer of the same thickness as the top model layer
c        9/26/01   Revised the calculation of vertical velocity (W) to be
c                  consistent with the way the implicit solver uses W
c        4/10/03   X/Y advection now uses layer-dependent timestep
c                            
c     Input arguments:  
c        igrid               grid index 
c        xyordr              x/y advection order
c        ncol                number of columns 
c        nrow                number of rows 
c        nlay                number of layers
c        nadv                number of sub-steps per timestep
c        deltat              timestep (s)
c        dx                  cell size in x direction (m)
c        dy                  cell size in y direction (m)
c        depth               layer depth (m) 
c        phpt                time-rate change of layer interface height (m/s)
c        pppt                time-rate change of pressure (mb/s)
c        ptpt                time-rate change of temperature (K/s)
c        windu               u-component windspeed (m/s)
c        windv               v-component windspeed (m/s)
c        tempk               temperature (k)
c        press               pressure (mb)
c        mapscl              map-scale factor at cell centroid
c             
c     Output arguments:  
c        dilut               dilution rate (m/s) 
c        entrn               entrainment rate (m/s) 
c             
c     Routines Called:  
c        none 
c             
c     Called by:  
c        EMISTRNS 
c 
      include 'camx.prm'
      include 'bndary.com'
      include 'flags.com'
c
      integer xyordr
      integer nadv(nlay)
      real phpt(ncol,nrow,nlay), 
     &     windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &     depth(ncol,nrow,nlay),tempk(ncol,nrow,nlay),
     &     press(ncol,nrow,nlay),dilut(ncol,nrow,nlay),
     &     entrn(ncol,nrow,nlay),pppt(ncol,nrow,nlay),
     &     ptpt(ncol,nrow,nlay),mapscl(ncol,nrow),dx(nrow)
      real rho(MXCOLA,MXROWA,MXLAYA+1),rhon(MXCOLA,MXROWA,MXLAYA),
     &     fluxx(MXCOLA,MXROWA,MXLAYA),fluxy(MXCOLA,MXROWA,MXLAYA)
      real flxarr(MX1D),v1d(MX1D),c1d(MX1D),m1d(MX1D),saflux(MX1D)
      real ronxt(MXLAYA+1)
      real fpc(MX1D),fmc(MX1D)
      real*8 flux1,flux2
c
c========================= Process Analysis Begin ==============================
c-----These are not used here, but are needed to give the correct
c     number of arguments in the call to the horizontal advection solvers
c
      real fc1(MX1D),fc2(MX1D)
c
c========================= Process Analysis End ================================
c
c
c-----Entry point 
c
      do 30 k = 1,nlay
c
c-----Load current atmospheric density (rho) and advected density (rhon)
c
        do j = 1,nrow
          do i = 1,ncol
            rho(i,j,k) = press(i,j,k)/tempk(i,j,k)
            rhon(i,j,k) = rho(i,j,k)
            if (k.eq.nlay) then
              drodz = 2.*(rho(i,j,k) - rho(i,j,k-1))/
     &                (depth(i,j,k) + depth(i,j,k-1))
              rho(i,j,nlay+1) = rho(i,j,k) + drodz*depth(i,j,k)
            endif
          enddo
        enddo
c
c-----Determine atmospheric density flux in x-direction
c
        if (xyordr.eq.0) goto 200
 100    continue
        do 10 j = 2,nrow-1
          i1 = 1
          i2 = ncol
          if (igrid.eq.1) then
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
     &        v1d(l) = 2.*v1d(l)/(mapscl(i+1,j) + mapscl(i,j)) 
            c1d(l) = rhon(i,j,k)*dy*depth(i,j,k)
            m1d(l) = mapscl(i,j)*mapscl(i,j)
          enddo
          nn = i2 - i1 + 1
c
          do i = i1,i2-1
            fluxx(i,j,k) = 0.
          enddo
          dtuse = deltat/nadv(k)
          do istep = 1,nadv(k)
            if (iadvct.eq.2) then
              call hadvbot(nn,dtuse,dx(j),c1d,v1d,m1d,flxarr,flux1,
     &                                flux2,saflux,fpc,fmc,fc1,fc2)
            elseif( iadvct .eq. 3) then
              call hadvppm(nn,dtuse,dx(j),c1d,v1d,m1d,flxarr,flux1,
     &                                        flux2,saflux,fc1,fc2)
            endif
c
            l = 0
            do i = i1,i2-1
              l = l + 1
              fluxx(i,j,k) = fluxx(i,j,k) + flxarr(l)/nadv(k)
              if (i.gt.i1) rhon(i,j,k) = c1d(l)/dy/depth(i,j,k)
            enddo
          enddo
 10     continue
        if (xyordr.eq.0) goto 30
c
c-----Determine atmospheric density flux in y-direction
c
 200    continue
        do 20 i = 2,ncol-1
          j1 = 1
          j2 = nrow
          if (igrid.eq.1) then
            if (jbeg(i).eq.-999) goto 20
            j1 = jbeg(i) - 1
            j2 = jend(i) + 1
          endif
c
          l = 0
          do j = j1,j2
            l = l + 1
            v1d(l) = windv(i,j,k)
            if (j.lt.j2) 
     &        v1d(l) = 2.*v1d(l)/(mapscl(i,j+1) + mapscl(i,j)) 
            c1d(l) = rhon(i,j,k)*dx(j)*depth(i,j,k)
            m1d(l) = mapscl(i,j)*mapscl(i,j)
          enddo
          nn = j2 - j1 + 1
c
          do j = j1,j2-1
            fluxy(i,j,k) = 0.
          enddo
          dtuse = deltat/nadv(k) 
          do istep = 1,nadv(k)
            if (iadvct.eq.2) then
              call hadvbot(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1,
     &                               flux2,saflux,fpc,fmc,fc1,fc2)
            elseif( iadvct .eq. 3) then
              call hadvppm(nn,dtuse,dy,c1d,v1d,m1d,flxarr,flux1,
     &                                     flux2,saflux,fc1,fc2)
            endif
c
            l = 0
            do j = j1,j2-1
              l = l + 1
              fluxy(i,j,k) = fluxy(i,j,k) + flxarr(l)/nadv(k)
              if (j.gt.j1) rhon(i,j,k) = c1d(l)/dx(j)/depth(i,j,k)
            enddo
          enddo
 20     continue
        if (xyordr.eq.0) goto 100
 30   continue
c
c-----Loop over grid to calculate vertical velocity, dilution and
c     entrainment rates
c             
      do 60 j = 2,nrow-1 
        i1 = 2 
        i2 = ncol - 1 
        if (igrid.eq.1) then 
          if(ibeg(j).eq.-999) goto 60 
          i1 = ibeg(j) 
          i2 = iend(j) 
        endif 
        do 50 i = i1,i2
c
c-----Diagnose density at end of timestep
c
          rhow = 0.
          do k = 1,nlay
            pnxt = press(i,j,k) + deltat*pppt(i,j,k)
            tnxt = tempk(i,j,k) + deltat*ptpt(i,j,k)
            ronxt(k) = pnxt/tnxt
          enddo
c
c-----Calculate horizontal divergence
c
          do 40 k = 1,nlay 
            drhou = (fluxx(i,j,k) - fluxx(i-1,j,k))*mapscl(i,j)**2
            drhov = (fluxy(i,j,k) - fluxy(i,j-1,k))*mapscl(i,j)**2
            div = (drhou + drhov)/dy/dx(j)/depth(i,j,k)
c
c-----Calculate actual density tendency
c
            drhodt = (ronxt(k) - rho(i,j,k))/deltat
            totrat = drhodt + div
            rhow = rhow - depth(i,j,k)*totrat
c
            if (rhow.ge.0.) then
              windw = rhow/ronxt(k) 
            else
              if (k.eq.nlay) then
                ronxt(k+1) = rho(i,j,k+1) + deltat*rhow/depth(i,j,nlay)
              endif
              windw = rhow/ronxt(k+1)
            endif
c             
c-----Calculate entrainment and dilution rates
c
            dilut(i,j,k) = phpt(i,j,k)
            if (k.gt.1) dilut(i,j,k) = phpt(i,j,k) - phpt(i,j,k-1) 
            entrn(i,j,k) = phpt(i,j,k) - windw
 40       continue 
 50     continue
 60   continue 
c
      return
      end
