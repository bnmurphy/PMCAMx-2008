      function areag(x1, x2, d, sigma)
      implicit none

c-----PMCAMx v3.01 020531
c
c --- calculate the fractional area of a Gaussian 
c     within the interval x1 to x2
c
c     Copyright 2002
c
c     corrected by bkoo (11/14/03)
c
      double precision areag, x1, x2, d, sigma
      double precision t, xerf, sqr2

      data sqr2 /1.414214/

c --- enforce x2 > x1

      if (x1 .gt. x2) then
        t = x1
        x1 = x2
        x2 = t
      elseif (x1 .eq. x2) then
        areag = 0.0
        return
      endif

c --- calculate area

      areag = 0.5*(xerf((x2-d)/sigma/sqr2) - xerf((x1-d)/sigma/sqr2))

      return
      end


      function xerf(x)
c
c     Returns the complementary error function erfc(x) with fractional error
c     everywhere less than 1.2e-7 (Numerical Recipes in Fortran, 2nd Ed., 1992)
c
c     modified by bkoo to return the error function erf(x) - (11/14/03)
c
      double precision xerf
      real erfcc
      double precision x
      real t,z
      z=abs(x)
      t=1./(1.+0.5*z)
      erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
     1     t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
     2     t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
      if (x.lt.0.) erfcc=2.-erfcc

      xerf = 1.-erfcc ! return erf(x)

      return
      end

