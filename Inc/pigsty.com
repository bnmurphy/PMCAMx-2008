c-----CAMx v4.02 030709
c 
c     PIGSTY.COM contains general Plume-in-Grid variables
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none
c
c-----------------------------------------------------------------------
c     Variables for PiG submodel:
c
c     npig     -- number of active PiG puffs
c     xlmax    -- maximum puff length (meters)
c     agemax   -- maximum puff age (hours)
c     ingrd    -- grid ID in which puff resides
c     idpig    -- stack ID from which puff is released
c     iipig    -- current i-location of puff
c     jjpig    -- current j-location of puff
c     xpig     -- current x-location of puff (km)
c     ypig     -- current y-location of puff (km)
c     zpig     -- current z-location of puff (meters)
c     xlength  -- puff length (meters)
c     axisy    -- width of puff (meters)
c     axisz    -- depth of puff (meters)
c     sigy     -- Gaussian standard deviation width (meters)
c     sigz     -- Gaussian standard deviation depth (meters)
c     lnewt    -- flag indicating newly released puff for transport
c     lnewg    -- flag indicating newly released puff for growth
c     puffmass -- puff pollutant mass (umol)
c     fmspig   -- puff volume parameter
c     agepig   -- puff age (seconds)
c     thetapig -- puff potential temperature (K)
c     npigon   -- number of active puffs in grid
c     idpigvar -- PiG variable identifier for NetCDF I/O
c     kntpig   -- puff counter for NetCDF I/O
c     pufftop  -- height of top of puff
c     puffbot  -- height of bottom of puff
c     cncradp  -- initial guess for radicals in puff reactors
c     nreactr  -- number of reactions in puff chemistry
c     pgmserr  -- error in mass due to aborted pig dumping(umol)
c     puffpol  -- puff pollutant mass (full chemistry - IRON PiG)
c     kopfail  -- counter for puff overlay failures
c
c-----------------------------------------------------------------------
c
      logical lnewt,lnewg
      common /pig1/ npig,xlmax,agemax,
     &              ingrd(MXPIG),idpig(MXPIG),iipig(MXPIG),jjpig(MXPIG),
     &              xpig(MXPIG),ypig(MXPIG),zpig(MXPIG),
     &              xlength(MXPIG),axisy(MXPIG),axisz(MXPIG),
     &              sigy(MXPIG),sigz(MXPIG),lnewt(MXPIG),lnewg(MXPIG),
     &              puffmass(4,MXPIG),fmspig(MXPIG),agepig(MXPIG),
     &              thetapig(MXPIG),npigon(MXGRID),
     &              idpigvar(20), kntpig, pufftop(MXPIG), 
     &              puffbot(MXPIG),
     &              cncradp(MXPIG,MXRECTR,MXRADCL), nreactr,
     &              pgmserr(MXSPEC),puffpol(MXSPEC,MXRECTR,MXPIG),
     &              kopfail
