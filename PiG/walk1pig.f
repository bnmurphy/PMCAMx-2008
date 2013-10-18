      subroutine walk1pig(kount,kpuff,dt,igrd,igrdo,n,nx,ny,nz,windu,
     &                    windv,height)
c
c-----CAMx v4.02 030709
c
c     WALK1PIG calculates new pig location by integrating dx/u(x) = dt
c     where u(x) = u(i-1) + (u(i) - u(i-1))*(x/deltax).  Each puff is
c     transported the entire dt or to the nearest cell interface,
c     whichever occurs first.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        kount               counter
c        kpuff               Puff layer index (kount > 1)
c        dt                  remaining timestep (s)
c        igrd                grid index
c        igrdo               previous grid index
c        n                   puff index
c        nx                  number of cells in x-direction
c        ny                  number of cells in y-direction
c        nz                  number of layers
c        windu               wind speed in x-direction (m/s)
c        windv               wind speed in y-direction (m/s)
c        height              layer height (m)
c
c     Output arguments:
c        kpuff               Puff layer index (kount = 1)
c
c     Subroutines Called:
c        none
c
c     Called by:
c        PIGWALK
c
      include "camx.prm"
      include "grid.com"
      include "pigsty.com"
      include "bndary.com"
      include "filunit.com"
c
      dimension windu(nx,ny,nz),windv(nx,ny,nz),height(nx,ny,nz)
c 
c-----Entry point
c
      iedge = 0
c
c-----Calculate current location
c
      i = iipig(n)
      j = jjpig(n)
      if (kount.gt.1 .and. igrd.eq.igrdo) goto 10
      igrdo = igrd
      do kpuff = 1,nlay(igrd)
        if (height(i,j,kpuff).gt.zpig(n)) goto 10
      enddo
      kpuff = nlay(igrd)
  10  k = kpuff
c
c-----Calculate puff distance from West/South cell edge
c
      if (igrd.eq.1) then
        wedge = (i - 1)*delx
        sedge = (j - 1)*dely
      else
        wedge = (inst1(igrd) - 1)*delx + (i - 2)*delx/meshold(igrd)
        sedge = (jnst1(igrd) - 1)*dely + (j - 2)*dely/meshold(igrd)
      endif
      xcell0 = (xpig(n) - wedge)/delx*deltax(j,igrd)*meshold(igrd)
      ycell0 = (ypig(n) - sedge)/dely*deltay(igrd)*meshold(igrd)
c
c-----Check if the wind is constant across the cell, which determines the
c     equations used for calculating the distance (or time)
c
      iu = INT( ABS( windu(i,j,k) - windu(i-1,j,k))*100. )
      iv = INT( ABS( windv(i,j,k) - windv(i,j-1,k))*100. )
c
c-----Compute the new position after the whole time step
c
      if (iu.eq.0) then
        xcell = xcell0 + windu(i,j,k)*dt
      else
        ax = windu(i-1,j,k)
        bx = (windu(i,j,k) - windu(i-1,j,k))/deltax(j,igrd)
        xcell = (-ax + (ax + bx*xcell0)*exp(bx*dt))/bx
      endif
      if (iv.eq.0) then
        ycell = ycell0 + windv(i,j,k)*dt
      else
        ay = windv(i,j-1,k)
        by = (windv(i,j,k) - windv(i,j-1,k))/deltay(igrd)
        ycell = (-ay + (ay + by*ycell0)*exp(by*dt))/by
      endif
c
c-----If the new position is beyond the boundaries of cell, compute
c     the time needed to reach the boundary
c
      dtx = dt
      dty = dt
      if (xcell.gt.deltax(j,igrd)) then
        iedge = 1
        if (iu.eq.0) then
          dtx = (deltax(j,igrd)*1.001 - xcell0)/windu(i,j,k)
        else
          dtx = alog((ax + bx*deltax(j,igrd)*1.001)
     &              /(ax + bx*xcell0))/bx
        endif
      elseif (xcell.lt.0.) then
        iedge = 1
        if (iu.eq.0) then
          dtx = (-deltax(j,igrd)*0.001 - xcell0)/windu(i,j,k)
        else
          dtx = alog((ax - bx*deltax(j,igrd)*0.001)
     &              /(ax + bx*xcell0))/bx
        endif
      endif
      if (ycell.gt.deltay(igrd)) then
        iedge = 1
        if (iv.eq.0) then
          dty = (deltay(igrd)*1.001 - ycell0)/windv(i,j,k)
        else
          dty = alog((ay + by*deltay(igrd)*1.001)
     &              /(ay + by*ycell0))/by
        endif
      elseif (ycell.lt.0.) then
        iedge = 1
        if (iv.eq.0) then
          dty = (-deltay(igrd)*0.001 - ycell0)/windv(i,j,k)
        else
          dty = alog((ay - by*deltay(igrd)*0.001)
     &              /(ay + by*ycell0))/by
        endif
      endif
c
c-----If the puff reaches a boundary, compute the remaining time and new 
c     position
c
      if (iedge.eq.1) then
        dttmp = amin1(dtx,dty,dt)
        if (iu.eq.0) then
          xcell = xcell0 + windu(i,j,k)*dttmp
        else
          xcell = (-ax + (ax + bx*xcell0)*exp(bx*dttmp))/bx
        endif
        if (iv.eq.0) then
          ycell = ycell0 + windv(i,j,k)*dttmp
        else
          ycell = (-ay + (ay + by*ycell0)*exp(by*dttmp))/by
        endif
        dt = dt - dttmp
      else
        dt = 0.
      endif
      xpig(n) = wedge + xcell/deltax(j,igrd)*delx/meshold(igrd)
      ypig(n) = sedge + ycell/deltay(igrd)*dely/meshold(igrd)
c
c-----If a puff reaches a boundary, compute new ingrd, iipig, jjpig
c
      if (iedge.ne.0) then
        i = 1 + INT( xpig(n)/delx )
        j = 1 + INT( ypig(n)/dely )
c
        ingrd(n) = 1
        do ip = 1,ngrid
          do ic = 1,nchdrn(ip)
            igrd = idchdrn(ic,ip)
            ig = mapgrd(igrd)
            if (i.ge.inst1(ig) .and. i.le.inst2(ig) .and.             
     &          j.ge.jnst1(ig) .and. j.le.jnst2(ig)) 
     &      ingrd(n) = igrd
          enddo
        enddo
c
        igrd = ingrd(n)
        if (igrd.eq.1) then
          iipig(n) = i
          jjpig(n) = j
          if (jbeg(i).eq.-999 .or. j.lt.jbeg(i) .or. j.gt.jend(i) .or.
     &        ibeg(j).eq.-999 .or. i.lt.ibeg(j) .or. i.gt.iend(j)) then
             write(idiag,*)'Puff beyond Coarse Grid boundary ',n
             ingrd(n) = 0
             dt = 0.
          endif 
        else
          xtmp = xpig(n) - (inst1(igrd) - 1)*delx
          ytmp = ypig(n) - (jnst1(igrd) - 1)*dely
          iipig(n) = 2 + INT( xtmp/delx*FLOAT( meshold(igrd) ) )
          jjpig(n) = 2 + INT( ytmp/dely*FLOAT( meshold(igrd) ) )
        endif
      endif
c
      return
      end
