c-----CAMx v4.02 030709
c  
c     DEPOSIT.COM contains arrays for dry deposition 
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c        none  
c 
c-----------------------------------------------------------------------
c     Arrays defining the Wesely (1989) resistance model:
c
c     z0lu    -- surface roughness length by landuse (meter)
c     istress -- vegetation stress code by landuse/season
c     rj      -- baseline minimum stomatal resistance (s/m) 
c     rlu     -- baseline upper canopy (cuticle) resistance (s/m) 
c     rac     -- baseline canopy height/density resistance (s/m) 
c     rlcs    -- baseline SO2 lower canopy resistance (s/m) 
c     rlco    -- baseline O3 lower canopy resistance (s/m) 
c     rgss    -- baseline SO2 ground surface resistance (s/m) 
c     rgso    -- baseline O3 ground surface resistance (s/m)
c-----------------------------------------------------------------------
c
      common /deposit/ z0lu(NLU),istress(NLU,5),rj(NLU,5),rlu(NLU,5),
     &                 rac(NLU,5),rlcs(NLU,5),rlco(NLU,5),rgss(NLU,5),
     &                 rgso(NLU,5)
c
c     sectional deposition for M4 (11/15/03)
c
      parameter (nsecfin = 5)
      parameter (nseccrs = 3)
      parameter (nsecdep = nsecfin + nseccrs)
c
      real diadep(nsecdep)    ! log-mean diameter for each bin
      real wfdep(nsecdep)     ! weighting-factor for each bin
c
c     wfdep(1) + ... + wfdep(nsecfin) = 1
c     wfdep(nsecfin+1) + ... + wfdep(nsecdep) = 1
c
      data diadep /0.055013, 0.12847,     0.3, 0.70057,  1.6360,
     &               3.8204,  8.9214,  20.833/
      data wfdep  /  0.0294,  0.2094,  0.4061,  0.2216,  0.1335,
     &               0.4078,  0.4205,  0.1717/

