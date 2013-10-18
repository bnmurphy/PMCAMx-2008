c-----CAMx v4.02 030709
c  
c     AHOMAP.COM contains albedo/haze/ozone index maps and value classes
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c        none  
c
c-----------------------------------------------------------------------
c     Variables to contain the albedo/haze/ozone data:
c
c     albcl   -- albedo column lookup table
c     hazcl   -- haze column lookup table
c     ozcl    -- ozone column lookup table
c     icdalb  -- index of the albedo column to use for each cell
c     icdhaz  -- index of the haze column to use for each cell
c     icdozn  -- index of the ozone column to use for each cell
c     lrdalb  -- flag to determine if data for grid was read
c
c-----------------------------------------------------------------------
c
      real      albcl(NALB)
      real      hazcl(NHAZE)
      real      ozcl(NOZN)
      integer   icdalb(MXVEC2D)
      integer   icdhaz(MXVEC2D)
      integer   icdozn(MXVEC2D)
      logical   lrdalb(MXGRID)
c 
      common /ahomap/ albcl, hazcl, ozcl, icdalb, icdhaz, icdozn,
     &                lrdalb
