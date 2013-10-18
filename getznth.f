      subroutine getznth(lat,lon,time,date,itzon,zenith,ldark)
c
c-----CAMx v4.02 030709
c
c     GETZNTH calculates solar zenith and sets the darkness flag
c
c     Based on equations given by Paltridge and Platt [1976]
c        "Radiative Processes in Meteorology and Climatology",
c         Elsevier, pp. 62,63.
c     Fourier coefficients originally from:
c        Spencer, J.W., 1971, Fourier series representation of the position 
c        of the sun, Search, 2:172.
c     NOTE: this approximate program has no changes from year to year.
c
c     Portions Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        lat                 latitude (deg) 
c        lon                 longitude (deg)
c        time                Current simulation time (HHMM)
c        date                Current simulaiton date (YYJJJ)
c        itzon               Model time zone (5=EST,6=CST,7=MST,8=PST)
c
c     Output arguments: 
c        zenith              solar zenith angle (90 if zenith > 88)
c        ldark               darkness flag
c
c     Routines Called:
c        none
c
c     Called by:
c        CHEMDRIV
c        DRYDEP
c        RADDRIVR
c
      parameter (pi=3.14159265) 
      parameter (dr=pi/180.0)
c
      real lat,lon
      integer date
      real lbut,lzut
      real rlt
      real d, tz, rdecl, eqr, eqh, zpt
      real csz, zr
      real sintz, costz, sin2tz, cos2tz, sin3tz, cos3tz
      integer iyear, ijd
      logical ldark
c
c-----Entry point
c
c-----Convert to radians
c
      rlt = lat*dr
c
c-----Parse date and convert date and time to Greenwich
c
      iyear = date/1000
      ijd   = mod(date,1000)
      ut = aint(time/100.) + amod(time,100.)/60. + float(itzon)
      if (ut .ge. 24.) then
        ut = ut - 24.
        ijd = ijd + 1
        if ((ijd.gt.365 .and. mod(iyear,4).ne.0) .or.
     &      (ijd.gt.366 .and. mod(iyear,4).eq.0)) ijd = 1
      endif
c
c-----Calculate decimal Julian day from start of year:
c
      d = float(ijd-1) + ut/24.
c
c-----Equation 3.8 for "day-angle"
c
      tz = 2.*pi*d/365.
c
c-----Calculate sine and cosine from addition theoremes for 
c     better performance;  the computation of sin2tz,
c     sin3tz, cos2tz and cos3tz is about 5-6 times faster
c     than the evaluation of the intrinsic functions 
c
c     It is SIN(x+y) = SIN(x)*COS(y)+COS(x)*SIN(y)
c     and   COS(x+y) = COS(x)*COS(y)-SIN(x)*SIN(y)
c
c     sintz  = SIN(tz)      costz  = COS(tz)
c     sin2tz = SIN(2.*tz)   cos2tz = SIN(2.*tz)
c     sin3tz = SIN(3.*tz)   cos3tz = COS(3.*tz)
c
      sintz = SIN(tz)
      costz = COS(tz)
      sin2tz = 2.*sintz*costz
      cos2tz = costz*costz-sintz*sintz
      sin3tz = sintz*cos2tz + costz*sin2tz
      cos3tz = costz*cos2tz - sintz*sin2tz
c
c----Equation 3.7 for declination in radians
c
      rdecl = 0.006918 - 0.399912*costz  + 0.070257*sintz 
     &                 - 0.006758*cos2tz + 0.000907*sin2tz    
     &                 - 0.002697*cos3tz + 0.001480*sin3tz
c
c-----Equation 3.11 for Equation of time  in radians
c
      eqr   = 0.000075 + 0.001868*costz  - 0.032077*sintz
     $		       - 0.014615*cos2tz - 0.040849*sin2tz
c
c-----Convert equation of time to hours:
c
      eqh = eqr*24./(2.*pi) 
c
c-----Calculate local hour angle (hours):
c
      lbut = 12. - eqh - lon*24./360 
c
c-----Convert to angle from UT
c
      lzut = 15.*(ut - lbut)
      zpt = lzut*dr
c
c-----Equation 2.4 for cosine of zenith angle 
c
      csz = sin(rlt)*sin(rdecl) + cos(rlt)*cos(rdecl)*cos(zpt)
      zr = acos(csz)
      zenith = zr/dr
      ldark = .false.
      if (zenith.gt.88.) then
        zenith = 90.
        ldark = .true.
      endif
c
      return
      end
