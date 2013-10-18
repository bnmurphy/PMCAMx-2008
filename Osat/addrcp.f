
      subroutine addrcp(igrd,nox,noy,nspsa,nspc,saconc,conc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine adds the concentrations for each receptor and to
c     the running sum of concentrations.  This data will then be
c     averaged over all time steps in an hour to get hourly
c     average concentration at the receptor locations. 
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c         igrd   I  grid number for this grid
c         nox    I  number of X cells in grid
c         noy    I  number of Y cells in grid
c         nspsa  I  number of tracer species
c         nspc   I  number of regular model species
c         saconc R  traxcer concentrations 
c         conc   R  concentrations from the regular model
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c         
c    04/02/03  --gwilson--   Added grid number to recptors defined by
c                            cell index
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      real      dxcell
      real      dycell
      integer   nox
      integer   noy
      integer   nspc
      integer   nspsa
      real      conc(nox,noy,nspc)
      real      saconc(nox,noy,nspsa)
      integer   igrd
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   ircp, i, j, icl, jcl
      integer   ixlow, ixhigh, iylow, iyhigh
      real      xnest, ynest
      real      xlow, xhigh, ylow, yhigh, xylow, xyhigh
      real      xcelwd, ycelwd, xcentr, ycentr
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- loop over all receptors, except the hourly peak (recptor 1) ---
c
      do 10 ircp=2,nrecep
c
c  --- skip receptor if it is not in this grid ---
c
         if( igrdrcp(ircp) .NE. igrd ) goto 10
c
c  --- if receptor is a single cell type just add the concentrations ---
c
          if( idrcp(ircp) .EQ. IDCEL ) then
              icl = irecep(ircp,1)
              jcl = jrecep(ircp,1)
              do 20 i=1,nsaspc
                 conrcp(i,ircp) = conrcp(i,ircp) + saconc(icl,jcl,i)
   20         continue
              do 30 i=1,nspc
                 if( lvocsp(i) ) then
                     rcpvoc(ircp) = rcpvoc(ircp) + 
     &                                     conc(icl,jcl,i) * crbnum(i)
                 else if( lnoxsp(i) ) then
                     rcpnox(ircp) = rcpnox(ircp) + conc(icl,jcl,i)
                 else if( lo3sp(i) ) then
                     rcpo3(ircp) = rcpo3(ircp) + conc(icl,jcl,i)
                 endif
   30         continue
c
c   --- if receptor is a cell average type calculate the average of 
c       each tracer concentration ---
c
          else if( idrcp(ircp) .EQ. IDAVG ) then
              do 40 j=1,nclrcp(ircp)
                 icl = irecep(ircp,j)
                 jcl = jrecep(ircp,j)
                 do 50 i=1,nsaspc
                     conrcp(i,ircp) = conrcp(i,ircp) + 
     &                               saconc(icl,jcl,i) / nclrcp(ircp)
   50            continue
                 do 60 i=1,nspc
                     if( lvocsp(i) ) then
                        rcpvoc(ircp) = rcpvoc(ircp) + 
     &                     conc(icl,jcl,i) * crbnum(i) / nclrcp(ircp)
                     else if( lnoxsp(i) ) then
                        rcpnox(ircp) = rcpnox(ircp) + 
     &                                conc(icl,jcl,i) / nclrcp(ircp)
                     else if( lo3sp(i) ) then
                        rcpo3(ircp) = rcpo3(ircp) + 
     &                                conc(icl,jcl,i) / nclrcp(ircp)
                     endif
   60            continue
   40         continue
c
c  --- else if receptor is a point type, perform bilinear interpolation
c      at the location and add to receptor total ----
c
          else if( idrcp(ircp) .EQ. IDPNT ) then
c
c   --- calculate the cell containing the receptor ---
c
              icl = irecep(ircp,1)
              jcl = jrecep(ircp,1)
c
c   --- if the cell is on the boundary, don't do interpolation ---
c       (why would you have a receptor in a bounday cell?)
c
              if( icl .LE. 1 .OR. icl .GE. nox .OR. 
     &                           jcl .LE. 1 .OR. jcl .GE. noy ) then
                 do 70 i=1,nsaspc
                   conrcp(i,ircp) = conrcp(i,ircp) + saconc(icl,jcl,i)
   70            continue
                 do 80 i=1,nspc
                    if( lvocsp(i) ) then
                        rcpvoc(ircp) = rcpvoc(ircp) + 
     &                                    conc(icl,jcl,i) * crbnum(i)
                    else if( lnoxsp(i) ) then
                        rcpnox(ircp) = rcpnox(ircp) + conc(icl,jcl,i)
                    else if( lo3sp(i) ) then
                        rcpo3(ircp) = rcpo3(ircp) + conc(icl,jcl,i)
                    endif
   80            continue
                 goto 10
              endif
c
c   --- find the locations of centers surrounding the receptor ---
c
              dxcell = delx / meshold(igrd)
              dycell = dely / meshold(igrd)
              xnest = xorg + FLOAT(inst1(igrd)-1)*delx - dxcell
              ynest = yorg + FLOAT(jnst1(igrd)-1)*dely - dycell
              xcentr = (FLOAT(icl)-1.5)*dxcell + xnest
              if( recepx(ircp) - xcentr .LE. dxcell ) then
                 ixlow = icl-1
                 xlow = xcentr
                 ixhigh = icl
                 xhigh =  (FLOAT(icl)-0.5)*dxcell + xnest
              else
                 ixlow = icl
                 xlow =  (FLOAT(icl)-0.5)*dxcell + xnest
                 ixhigh = icl+1
                 xhigh =  (FLOAT(icl)+0.5)*dxcell + xnest
              endif
              ycentr = (FLOAT(jcl)-1.5)*dycell + ynest
              if( recepy(ircp) - ycentr .LE. dycell ) then
                 iylow = jcl-1
                 ylow = ycentr
                 iyhigh = jcl
                 yhigh =  (FLOAT(jcl)-0.5)*dycell + ynest
              else
                 iylow = jcl
                 ylow =  (FLOAT(jcl)-0.5)*dycell + ynest
                 iyhigh = jcl+1
                 yhigh =  (FLOAT(jcl)+0.5)*dycell + ynest
              endif
              xcelwd = xhigh - xlow
              ycelwd = yhigh - ylow
              do 90 i=1,nsaspc
c
c  --- interpolate to the X position, at both Y positions ---
c
                  xylow = (recepx(ircp) - xlow) / xcelwd * 
     &                                           saconc(ixlow,iylow,i)
                  xylow = xylow + (xhigh - recepx(ircp)) / xcelwd * 
     &                                          saconc(ixhigh,iylow,i)
                  xyhigh = (recepx(ircp) - xlow) / dxcell * 
     &                                          saconc(ixlow,iyhigh,i)
                  xyhigh = xyhigh + (xhigh - recepx(ircp)) / 
     &                                dxcell * saconc(ixhigh,iyhigh,i)
c
c  --- interpolate the the exact position using the high/low X values ----
c
                 conrcp(i,ircp) = conrcp(i,ircp) + 
     &                     (recepy(ircp) - ylow)/ycelwd * xylow +
     &                            (yhigh - recepy(ircp))/ycelwd * xyhigh
   90         continue
              do 11 i=1,nspc
c
c  --- interpolate to the X position, at both Y positions ---
c
                  xylow = (recepx(ircp) - xlow) / xcelwd * 
     &                                           conc(ixlow,iylow,i)
                  xylow = xylow + (xhigh - recepx(ircp)) / xcelwd * 
     &                                           conc(ixhigh,iylow,i)
                  xyhigh = (recepx(ircp) - xlow) / dxcell * 
     &                                           conc(ixlow,iyhigh,i)
                  xyhigh = xyhigh + (xhigh - recepx(ircp)) /
     &                                 dxcell * conc(ixhigh,iyhigh,i)
c
c  --- interpolate the the exact position using the high/low X values ----
c
                 if( lvocsp(i) ) then
                     rcpvoc(ircp) = rcpvoc(ircp) + crbnum(i) *
     &                     ( (recepy(ircp) - ylow)/ycelwd * xylow +
     &                          (yhigh - recepy(ircp))/ycelwd * xyhigh )
             
                 else if( lnoxsp(i) ) then
                     rcpnox(ircp) = rcpnox(ircp) + 
     &                      (recepy(ircp) - ylow)/ycelwd * xylow +
     &                          (yhigh - recepy(ircp))/ycelwd * xyhigh 
                 else if( lo3sp(i) ) then
                     rcpo3(ircp) = rcpo3(ircp) + 
     &                      (recepy(ircp) - ylow)/ycelwd * xylow +
     &                          (yhigh - recepy(ircp))/ycelwd * xyhigh 
                 endif
   11         continue
          endif
c
c  --- next receptor ----
c
   10 continue
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
c
