c-----CAMx v4.02 030709
c 
c     FLAGS.COM contains all model flags except those related to chemistry
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        4/26/99   Removed LSMO, LBOT flags and replaced with IADVCT flag 
c        7/05/02   Removed lpig and replaced it with ipigflg
c        4/2/03    Removed option for UAM-V type cloud adjustment
c
c-----------------------------------------------------------------------
c     Integer flags:
c
c     ipigflg   -- PiG flag - GREASD or IRON
c     iadvct  -- Advection solver flag
c                (1 = Smolarkiewicz, 2 = Bott, 3 = Piecewise Parabolic Method)
c
c     Logical flags:
c
c     lrstrt  -- simulation restart flag
c     lchem   -- chemistry flag
c     ldry    -- dry deposition flag
c     lwet    -- wet deposition flag
c     lutm    -- uTM projection flag
c     llatlon -- lat/lon coordinates flag
c     lpolar  -- polar stereographic projection flag
c     lambrt  -- Lambert conformal projection flag
c     lstagw  -- staggered input wind field flag
c     larsrc  -- area source flag
c     lptsrc  -- point source flag
c     le1day  -- 1-day emission input flag
c     lairqul -- output average file flag
c     l3davg  -- 3-D output average file flag
c-----------------------------------------------------------------------
c
      integer ipigflg
      logical lrstrt,lchem,ldry,lwet,lutm,llatlon,lpolar,lambrt,
     &        lstagw,larsrc,lptsrc,le1day,lairqul,l3davg
c
      common /cflgs/ pigflg
      common /flags/ lrstrt,lchem,ldry,lwet,lutm,llatlon,lpolar,
     &               lambrt,lstagw,larsrc,lptsrc,le1day,lairqul,l3davg,
     &               iadvct,ipigflg
