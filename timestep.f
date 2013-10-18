      subroutine timestep()
c
c-----CAMx v4.02 030709
c
c     TIMESTEP computes the time step size for each grid 
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        11/30/99  Removed timestep calculation for horizontal diffusion
c        1/25/02   Revised I/O frequencies and max time step to minutes,
c                  added code to check winds at next update time to ensure
c                  CFL stability throughout period
c        4/10/03   Now calculating max timestep by layer, AND
c                  DTMAX moved from control file to parameter defined here
c
c     Input arguments:
c        none
c
c     Output Arguments:
c        none
c
c     Routines called:
c        GETDELT
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "camx.com"
      include "camxfld.com"
      include "bndary.com"
      include "grid.com"
c
      real dtmx1(MXLAYA),dtmx2(MXLAYA),dtlay(MXLAYA)
c
c-----Entry point
c
c-----Coarse grid: maximum time step size according to stability criterion
c
      igrd = 1
      call getdelt(ncol(1),nrow(1),nlay(1),deltax(1,1),deltay(1),
     &             windu(1),windv(1),height(1),dtmx1,kmax)
      call getdelt(ncol(1),nrow(1),nlay(1),deltax(1,1),deltay(1),
     &             unxt(1),vnxt(1),hnxt(1),dtmx2,kmax)
c
c-----Make sure an integer number of coarse grid time steps are to be
c     completed each hour and between I/O times
c
      dtio = 60.*amin1(60.,dtinp,dtems,dtout)
      dtmx = dtio
      do k = 1,nlay(1)
        dtlay(k) = amin1(dtmx1(k),dtmx2(k),dtio,60.*dtmax)
        if (k.le.kmax) dtmx = amin1(dtmx,dtlay(k))
      enddo
      nsteps = INT( 0.999*dtio/dtmx ) + 1
      deltat(1) = dtio/FLOAT( nsteps )
c
c-----Number of advection steps per layer in coarse grid
c
      do k = 1,nlay(1)
        nadv(k,1) = INT( 0.999*deltat(1)/dtlay(k) ) + 1
      enddo
c
c-----Fine grids: make sure an integer number of fine grid time steps 
c     are to be completed between coarser grid steps
c
      do ip = 1,ngrid
        do ic = 1,nchdrn(ip)
          igrd = idchdrn(ic,ip)
          call getdelt(ncol(igrd),nrow(igrd),nlay(igrd),
     &                 deltax(1,igrd),deltay(igrd),windu(iptr3d(igrd)),
     &                 windv(iptr3d(igrd)),height(iptr3d(igrd)),dtmx1,
     &                 kmax)
          call getdelt(ncol(igrd),nrow(igrd),nlay(igrd),
     &                 deltax(1,igrd),deltay(igrd),unxt(iptr3d(igrd)),
     &                 vnxt(iptr3d(igrd)),hnxt(iptr3d(igrd)),dtmx2,
     &                 kmax)
          dtmx = deltat(ip)
          do k = 1,nlay(igrd)
            dtlay(k) = amin1(dtmx1(k),dtmx2(k),deltat(ip))
            if (k.le.kmax) dtmx = amin1(dtmx,dtlay(k))
          enddo
          ntim(igrd) = INT( 0.999*deltat(ip)/dtmx ) + 1
          deltat(igrd) = deltat(ip)/FLOAT( ntim(igrd) )
c
c-----Number of advection steps per layer in this fine grid
c
          do k = 1,nlay(igrd)
            nadv(k,igrd) = INT( 0.999*deltat(igrd)/dtlay(k) ) + 1
          enddo
c
        enddo
      enddo
c
      return
      end
