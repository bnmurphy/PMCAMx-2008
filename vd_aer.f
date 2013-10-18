      subroutine vd_aer(z0,deltaz,psih,ustar,diam,rhop,ts,vd)
c
c-----CAMx v4.02 030709
c
c     VD_AER calculates a deposition velocity for a specific aerosol size
c     bin, grid cell, and land use category.  The parallel resistance approach 
c     of Slinn and Slinn (1980) is used, as implemented in UAM-AERO 
c     (STI, 1996).
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c 
c     Modifications: 
c        3/21/03       Modified Ra calculation to use layer 1 midpoint height
c                      rather than default 10 m.
c        11/7/03       Removed small bug in final Vd equation
c
c     Input arguments: 
c        z0                  surface roughness length (m)
c        deltaz              Layer 1 midpoint height (m)
c        psih                similarity stability correction term 
c        ustar               friction velocity (m/s)
c        diam                log-mean sectional aerosol diameter (m)
c        rhop                aerosol density (g/m3)
c        ts                  surface temperature (K)
c      
c     Output arguments: 
c        vd                  deposition velocity (m/s)
c      
c     Routines called: 
c        none 
c      
c     Called by: 
c        DRYDEP
c 
      data vk/0.4/, rmin/1.0/, xmfp/6.5e-8/, g/9.8/, vabs/1.81e-2/
      data boltz/1.38e-20/, pi/3.1415927/, vair/1.5e-5/
c
c-----Entry point
c
c-----Speed correction factor and sedimendation velocity
c
      power = amin1(7.6,0.55*diam/xmfp)
      scf = 1. + (2.514 + 0.8*exp(-power))*xmfp/diam
      vsed = rhop*g*(diam*diam)*scf/(18.*vabs)
c
c-----Brownian diffusivity and Schmidt number
c
      difbrwn = boltz*ts*scf/(3.*pi*vabs*diam)
      schmidt = vair/difbrwn
c
c-----Stokes number
c
      stokes = vsed*(ustar*ustar)/(vair*g)
c
c-----Compute atmospheric resistance, RA
c
      ra = (alog(deltaz/z0) - psih)/(vk*ustar)
      ra = amax1(ra,rmin)
c
c-----Compute the deposition layer resistance, RD
c
      sc23 = schmidt**(-2./3.)
      power = -3./stokes
      if (power.lt.-37.) then
        xinert = 10.**(-37.)
      else
        xinert = 10.**(power)
      endif
      rd = 1./(ustar*(sc23 + xinert))
      rd = amax1(rd,rmin)
c
c-----Final deposition velocity for this cell, land use, and aerosol size
c
      vd = vsed + 1./(ra + rd + ra*rd*vsed)
c
      return
      end
