      subroutine grdprep(nx,ny,cellon,cellat,mapscl)
c     
c---- CAMx v4.51 080522
c     
c     GRDPREP calculates coarse grid parameters for the following projections:
c     1. Universal Transverse Mercator (UTM)
c     2. Polar Stereographic (PSP)
c     3. Lambert Conformal (LCP)
c     4. lat/lon grid
c     
c     Copyright 1996-2008
c     ENVIRON International Corporation
c     
C     Modifications by Pavan (pavan_nandan_racherla@alumni.cmu.edu)
C     Updated Oct 13 2008
c     
c     Input arguments: 
c     nx                  number of columns
c     ny                  number of rows
c     
c     Output arguments: 
c     cellon              array of cell centroid longitude (deg)
c     cellat              array of cell centroid latitude (deg)
c     mapscl              array of map-scale factors centered on each cell
c     
c     Routines Called: 
c     UTMGEO
c     PSPGEO
c     LCPGEO
c     
c     Called by: 
c     STARTUP
c
 
      include 'camx.prm'
      include 'grid.com'
      include 'flags.com'

c     
      real cellon(nx,ny), cellat(nx,ny), mapscl(nx,ny)
c     
      data deg2km /111.1338/
      data deg2rad /0.01745329/

C     Pavan:
      real gs, cLat, cLon, eLon, eLat, wLon, wLat

c     
c-----Entry point
c     
c-----Cartesian coordinates are selected
c     Input cell sizes (delx and dely) are already in km
c     
      if( .NOT.llatlon .and. .NOT.lpolar ) then
c     
c-----Determine deltax and deltay (m) for coarse grid
c     
         deltay(1) = 1000.*dely
         do 10 j = 1,ny
            deltax(j,1) = 1000.*delx
c     
c-----Calculate lat/lon at cell centroids, and map-scale factor
c     
            yloc = yorg + dely*(float(j) - 0.5)
            do 15 i = 1,nx
               xloc = xorg + delx*(float(i) - 0.5)
               if (lutm) then
                  call utmgeo(1,iuzon,xloc,yloc,alon,alat)
                  xloce = xloc + delx/2. 
                  call utmgeo(1,iuzon,xloce,yloc,elon,elat)
                  xlocw = xloc - delx/2. 
                  call utmgeo(1,iuzon,xlocw,yloc,wlon,wlat)
               elseif (lambrt) then
                  call lcpgeo(1,ylatc,xlonc,tlat1,tlat2,xloc,yloc,alon,alat)
                  xloce = xloc + delx/2. 
                  call lcpgeo(1,ylatc,xlonc,tlat1,tlat2,xloce,yloc,elon,elat)
                  xlocw = xloc - delx/2. 
                  call lcpgeo(1,ylatc,xlonc,tlat1,tlat2,xlocw,yloc,wlon,wlat)
               endif
               dx = deg2km*(elon - wlon)*cos(deg2rad*alat)  
               dy = deg2km*(elat - wlat) 
               mapscl(i,j) = delx/sqrt(dx*dx + dy*dy) 
               cellon(i,j) = alon
               cellat(i,j) = alat
 15         continue
 10      continue

C     
c*****Polar Stereographic Projection
c*****Modified by Pavan
c*****
      else if ( LPOLAR ) then

         gs = 1.0e3 * delx
         deltay(1) = 1.0e3 * dely

         do J = 1, NY

            deltax(j, 1) = 1000. * delx

            do I = 1, NX
               
c--------------Lon/Lat of the cell centroid:
               call pspgeo(real(I), real(J), gs,
     $              psp_stdlon, psp_truelat1, psp_lon1, psp_lat1,
     $              cLat, cLon)               
               cellon(i,j) = cLon
               cellat(i,j) = cLat
               
c--------------Lon/Lat of the west face of the cell:
               call pspgeo(real(I), real(J+.5), gs,
     $              psp_stdlon, psp_truelat1, psp_lon1, psp_lat1,
     $              wLat, wLon)

c--------------Lon/Lat of the east face of the cell:
               call pspgeo(real(I+1), real(J+.5), gs,
     $              psp_stdlon, psp_truelat1, psp_lon1, psp_lat1,
     $              eLat, eLon)
               
c--------------Calculate the distortion factor for each cell:               
               dx = deg2km * (eLon - wLon) * cos(deg2rad * cLat)
               dy = deg2km * (eLat - wLat)               
               mapscl(I,J) = delx / sqrt(dx**2 + dy**2)
               
            end do            
         end do
c*****
c*****Pavan: end of my modifications.
         
c     
c-----Lat/lon coordinates are selected
c     
      else
c     
c-----Determine deltax and deltay (m) on coarse grid
c     
         deltay(1) = 1000.*deg2km*dely
         do 20 j = 1,ny
            alat = yorg + dely*(float(j) - 0.5)
            deltax(j,1) = 1000.*deg2km*delx*cos(alat*deg2rad)
c     
c-----Calculate lat/lon at cell centroids
c     
            do 25 i = 1,nx
               alon = xorg + delx*(float(i) - 0.5)
               cellon(i,j) = alon
               cellat(i,j) = alat
               mapscl(i,j) = 1.
 25         continue
 20      continue
      endif

      return
      end
      
