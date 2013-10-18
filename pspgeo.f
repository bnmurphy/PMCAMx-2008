      subroutine pspgeo(iway,polelon,polelat,xloc,yloc,qlon,qlat)
c
c-----CAMx v4.02 030709
c 
c     PSPGEO performs polar stereographic to geodetic (lat/lon) translation
c 
c     Portions Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c 
c     Modifications: 
c        none 
c 
c     Input arguments: 
c        iway                Conversion type
c                            0 = geodetic to Polar Stereographic
c                            1 = Polar Stereographic to geodetic
c        polelat             Latitude of projection pole (deg) 
c        polelon             Longitude of projection pole (deg,negative for W) 
c        xloc/yloc           Projection coordinates (km) 
c        qlon/qlat           latitude/longtitude (degrees) 
c 
c     Output arguments: 
c        xloc/yloc           Projection coordinates (km) 
c        qlon/qlat           latitude/longtitude (deg) 
c 
c     Routines called: 
c        none 
c 
c     Called by: 
c        GRDPREP 
c
      parameter (erad=6.367e6,erad2=1.2734e7,pi180=3.14159265/180.)
c
c-----Entry point
c
      xloc = 1000.*xloc
      yloc = 1000.*yloc
c
c-----Evaluate sine and cosine of latitude and longitude of pole point p.
c
      sinplat = sin(polelat * pi180)
      cosplat = cos(polelat * pi180)
      sinplon = sin(polelon * pi180)
      cosplon = cos(polelon * pi180)
c
c-----Compute (x3,y3,z3) coordinates of the pole point where the origin is the 
c     center of the earth, the z axis is through the north pole, the x axis is
c     through the equator and prime meridian, and the y axis is through the
c     equator and 90 E.
c
      x3p = erad * cosplat * cosplon
      y3p = erad * cosplat * sinplon
      z3p = erad * sinplat
c
c-----Polar stereographic to lat/lon conversion
c
      if (iway.eq.1) then
c
c----Compute distance d from given point R on the polar stereographic plane
c    to the pole point P:
c
        d = sqrt (xloc**2 + yloc**2)
c 
c-----Compute angle QCP where C is the center of the Earth.  This is twice 
c     angle QAP where A is the antipodal point.  Angle QAP is the same as
c     angle RAP:
c
        alpha = 2. * atan2(d,erad2)
c
c-----Compute zq, the height of Q relative to the polar stereographic plane:
c
        zq = erad * (cos(alpha) - 1.)
c
c-----Compute the parameter t which is the the distance ratio AQ:AR
c
        t = (erad2 + zq) / erad2
        xq = t * xloc
        yq = t * yloc
c
c-----Transform location of Q from (x,y,z) coordinates to (x3,y3,z3):
c
        x3q = x3p - xq * sinplon - yq * cosplon * sinplat
     +        + zq * cosplat * cosplon
        y3q = y3p + xq * cosplon - yq * sinplon * sinplat 
     +        + zq * cosplat * sinplon
        z3q = z3p + yq * cosplat + zq * sinplat
c
c-----Compute the latitude and longitude of Q:
c
        qlon = atan2(y3q,x3q) / pi180
        r3q = sqrt(x3q ** 2 + y3q ** 2)
        qlat = atan2(z3q,r3q) / pi180
c
c-----Lat/lon to Polar Stereographic conversion
c
      else
        sinqlat = sin(qlat * pi180)
        cosqlat = cos(qlat * pi180)
        sinqlon = sin(qlon * pi180)
        cosqlon = cos(qlon * pi180)
c
        z3q = erad * sinqlat
        x3q = erad * cosqlat * cosqlon
        y3q = erad * cosqlat * sinqlon
c
c-----Transform q point from (x3,y3,z3) coordinates in the above system to
c     polar stereographic coordinates (x,y,z):
c
        xq = - sinplon * (x3q-x3p) + cosplon * (y3q-y3p)  
        yq =   cosplat * (z3q-z3p)
     +       - sinplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )
        zq =   sinplat * (z3q-z3p)
     +       + cosplat * ( cosplon * (x3q-x3p) + sinplon * (y3q-y3p) )
c
c-----Parametric equation for line from antipodal point at (0,0,-2 erad) to
c     point q has the following parameter (t) value on the polar stereographic
c     plane:
c
        t = erad2 / (erad2 + zq)
        xloc = xq*t
        yloc = yq*t
      endif
c
      xloc = xloc/1000.
      yloc = yloc/1000.
c
      return
      end
