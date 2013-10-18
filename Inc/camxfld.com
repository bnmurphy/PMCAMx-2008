c-----CAMx v4.02 030709
c
c     CAMxFLD.COM contains all multidimensional fields for CAMx
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        1/10/03    Added deposition field array
c        1/23/03    New cloud parameters from cloud/rain file
c
c-----------------------------------------------------------------------
c     Variables for 2-dimensional fields:
c
c     cellon --  array of longitudes for the cell centroids (deg)
c     cellat --  array of lattitudes for the cell centroids (deg)
c     mapscl --  array of map-scale factors at cell centroids
c     tsurf  --  array for surface temperature field (K)
c     pspt   --  array for time-rate change of surface temperature (K/s) 
c-----------------------------------------------------------------------
c
      real   cellon(MXVEC2D)
      real   cellat(MXVEC2D)
      real   mapscl(MXVEC2D)
      real   tsurf(MXVEC2D)
      real   pspt(MXVEC2D)
c
      common /store2d/ cellon, cellat, mapscl, tsurf, pspt
c
c-----------------------------------------------------------------------
c     Variables for 3-dimensional fields:
c
c      windu  --  U-component of the wind field (m/s)
c      windv  --  V-component of the wind field (m/s)
c      pupt   --  time-rate change of U-component wind (m/s2)
c      pvpt   --  time-rate change of V-component wind (m/s2)
c      tempk  --  temperature field (K) 
c      ptpt   --  time-rate change of temperature (K/s)
c      press  --  pressure field (mb)
c      pppt   --  time-rate change of pressure (mb/s)
c      height --  layer interface height field (m)
c      phpt   --  time-rate change of layer interface height (m/s)
c      rkv    --  vertical diffusion coefficient field (m2/s)
c      water  --  water vapor concentration field (ppm)
c      fcloud --  vertically-accumulating cloud coverage field (fraction)
c      depth  --  layer depth field (m)
c      rkx    --  horizontal diffusion coefficient in X direction (m2/s)
c      rky    --  horizontal diffusion coefficient in Y direction (m2/s)
c      pi0    --  intitial pressure field (mb) -- Used for NETCDF I/O
c      aremis --  surface emissions (moles/hour or g/hour)
c      cwc    --  cloud water content (g/m3)
c      pwc    --  precipitation water content (g/m3)
c      cod    --  cloud optical depth
c      cldtrns -  Cloud energy transmission coefficient (fraction)
c-----------------------------------------------------------------------
c
      real   windu(MXVEC3D)
      real   windv(MXVEC3D)
      real   pupt(MXVEC3D)
      real   pvpt(MXVEC3D)
      real   tempk(MXVEC3D)
      real   ptpt(MXVEC3D)
      real   press(MXVEC3D)
      real   pppt(MXVEC3D)
      real   height(MXVEC3D)
      real   phpt(MXVEC3D)
      real   rkv(MXVEC3D)
      real   water(MXVEC3D)
      real   fcloud(MXVEC3D)
      real   depth(MXVEC3D)
      real   rkx(MXVEC3D)
      real   rky(MXVEC3D)
      real   pi0(MXVEC3D)
      real   aremis(MXVECEM)
      real   cwc(MXVEC3D)
      real   pwc(MXVEC3D)
      real   cod(MXVEC3D)
      real   cldtrns(MXVEC3D)
c
      common /store3d/ windu, pupt, windv, pvpt, tempk, ptpt, press,
     &                 pppt, height, phpt, rkv, water, fcloud, depth,
     &                 rkx, rky, pi0, aremis, cwc, pwc, cod, cldtrns
c
c-----------------------------------------------------------------------
c     Variables for 4-dimensional fields:
c
c      conc   --  species concentrations field (umol/m3)
c      avcnc  --  average species concentration (gas=ppm,other=ug/m3)
c      cncrad --  radical concentrations (ppm)
c-----------------------------------------------------------------------
c
      real   conc(MXVEC4D)
      real   avcnc(MXVEC4D)
      real   cncrad(MXVECRD)
c
      common /store4d/ conc, avcnc, cncrad
c
c-----------------------------------------------------------------------
c     Variables for calculating vertical transport:
c     NOTE:  These fields are over-written when processing each grid.
c
c     entrn  --  entrainment rate (m/s)
c     dilut  --  dilution rate (m/s)
c-----------------------------------------------------------------------
c   
      real   entrn(MXVEC3A)
      real   dilut(MXVEC3A)
c
      common /vrtrate/ entrn, dilut
c
c-----------------------------------------------------------------------
c     Variables for calculating depostion rates:
c
c     vdep   --  species-dependent deposition velocity field (m/s)
c     fsurf  --  fractional landuse cover field (fraction)
c     depfld --  2-D array containing dry, wet dep mass (mol/ha, g/ha) and
c                precip liquid concentration (mol/l, g/l)
c-----------------------------------------------------------------------
c
      real   vdep(MXVECEM)
      real   fsurf(MXVECLU)
      real   depfld(MXVECDP)
c
      common /depstn/  vdep, fsurf, depfld
c
c-----------------------------------------------------------------------
c     Variables for mass flux calculations:
c
c     xmass   -- current total grid mass by species (moles)
c     xmass0  -- initial total grid mass by species (moles)
c     armass  -- total area emissions mass by species (moles)
c     ptmass  -- total point emissions mass by species (moles)
c     fluxes  -- array of mass tranport by species (moles)
c     xmschem -- mass change due to chemistry by species (moles)
c     xmsold  -- total grid mass by species before last process (moles)
c     resid   -- residual (error) mass by species (moles)
c     xmsfin  -- current total fine grid mass by species (moles)
c     xmstmp  -- temporary total fine grid mass by species (moles)
c     pigdump -- total mass transferred from PiG to grid by species (moles)
c     pigmass -- total mass in PiG by species (moles)
c-----------------------------------------------------------------------
c
      real*8 xmass(MXSPEC,MXGRID)
      real*8 xmass0(MXSPEC,MXGRID)
      real*8 armass(MXSPEC,MXGRID)
      real*8 ptmass(MXSPEC,MXGRID)
      real*8 fluxes(MXSPEC*11,MXGRID)
      real*8 xmschem(MXSPEC,MXGRID)
      real*8 xmsold(MXSPEC,MXGRID)
      real*8 resid(MXSPEC,MXGRID)
      real*8 xmsfin(MXSPEC,MXGRID)
      real*8 xmstmp(MXSPEC,MXGRID)
      real*8 pigdump(MXSPEC,MXGRID)
      real*8 pigmass(MXSPEC,MXGRID)
c
      common /ms_flx/  xmass, xmass0, armass, ptmass, fluxes, xmschem,
     &                 xmsold, resid, xmsfin, xmstmp, pigdump, pigmass
c
c-----------------------------------------------------------------------
c     Variables for mass balance calculations for the extent of the
c     simulation:
c
c     tarmass  -- total area emissions mass by species (moles)
c     tptmass  -- total point emissions mass by species (moles)
c     tfluxes  -- array of mass tranport by species (moles)
c     txmschem -- mass change due to chemistry by species (moles)
c     tresid   -- residual (error) mass by species (moles)
c     txmsfin  -- current total fine grid mass by species (moles)
c-----------------------------------------------------------------------
c
      real*8 tarmass(MXSPEC,MXGRID)
      real*8 tptmass(MXSPEC,MXGRID)
      real*8 tfluxes(MXSPEC*12,MXGRID)
      real*8 tresid(MXSPEC,MXGRID)
      real*8 txmschem(MXSPEC,MXGRID)
      real*8 txmsfin(MXSPEC,MXGRID)
c
      common /flxacc/  tarmass, tptmass, tfluxes, tresid, txmschem, 
     &                 txmsfin
c
c-----------------------------------------------------------------------
c     Variables for storing the next hours meteorological fields:
c
c     hnxt   --  layer interface height at next hour (m)
c     pnxt   --  pressure field at next hour (mb)
c     unxt   --  U-component of the wind field at next hour (m/s)
c     vnxt   --  V-component of the wind field at next hour (m/s)
c     tnxt   --  temperature field at next hour (K) 
c     tsnxt  --  surface temperature field at next hour (K)
c-----------------------------------------------------------------------
c
      real   hnxt(MXVEC3D)
      real   pnxt(MXVEC3D)
      real   unxt(MXVEC3D)
      real   vnxt(MXVEC3D)
      real   tnxt(MXVEC3D)
      real   tsnxt(MXVEC2D)
c
      common /newdat/  hnxt, pnxt, unxt, vnxt, tnxt, tsnxt
