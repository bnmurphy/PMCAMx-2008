c-----CAMx v4.02 030709
c
c     PROCAN.COM contains general model variables
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c-----------------------------------------------------------------------
c     Parameters for Process Analysis:
c
c     MXPADOM     -- maximum number of domains allowed in Process Analysis
c     MXPACEL     -- maximum number of Process Analysis grid cells 
c     NPAPRC      -- number of processes tracked in Integrated Process Rates 
c     MXCPA       -- maximum number of "species" in Chemical Process Analysis (CPA)
c
c     ID of processes tracked in Integrated Process Rates (IPR)
c     (used in array CIPR):
c
c     IPR_INIT    -- initial concentration of each process analysis interval
c     IPR_CHEM    -- concentration change from chemistry
c     IPR_AEMIS   -- concentration change from area emissions
c     IPR_PTEMIS  -- concentration change from point source emissions
c     IPR_PIGEMIS -- concentration change from plume-in-grid (PiG)
c     IPR_WADV    -- concentration change due to advection at
c                    west boundary of grid cell
c     IPR_EADV    -- concentration change due to advection at
c                    east boundary of grid cell
c     IPR_SADV    -- concentration change due to advection at
c                    south boundary of grid cell
c     IPR_NADV    -- concentration change due to advection at
c                    north boundary of grid cell
c     IPR_BADV    -- concentration change due to entrainment
c                    from lower grid cell
c     IPR_TADV    -- concentration change due to entrainment
c                    from upper grid cell
c     IPR_DADV    -- concentration change due to diffusion of
c                    current grid cell
c     IPR_WDIF    -- concentration change due to diffusion at
c                    west boundary of grid cell
c     IPR_EDIF    -- concentration change due to diffusion at
c                    east boundary of grid cell
c     IPR_SDIF    -- concentration change due to diffusion at
c                    south boundary of grid cell
c     IPR_NDIF    -- concentration change due to diffusion at
c                    north boundary of grid cell
c     IPR_BDIF    -- concentration change due to diffusion at
c                    bottom boundary of grid cell
c     IPR_TDIF    -- concentration change due to diffusion at
c                    top boundary of grid cell
c     IPR_DDEP    -- concentration change from dry deposition
c     IPR_WDEP    -- concentration change from wet deposition
c     IPR_FAERO   -- concentration change from fast aero
c     IPR_FINAL   -- final concentration of each process analysis interval
c     IPR_CONV    -- umole/m3 <=> ppm conversion
c     IPR_VOL     -- average volume of the cell
c
c-----------------------------------------------------------------------
c
      integer   MXPADOM
      integer   MXPACEL  
      integer   NPAPRC
      integer   MXCPA
 
      integer   IPR_INIT
      integer   IPR_CHEM
      integer   IPR_AEMIS
      integer   IPR_PTEMIS
      integer   IPR_PIGEMIS
      integer   IPR_WADV
      integer   IPR_EADV
      integer   IPR_SADV
      integer   IPR_NADV
      integer   IPR_BADV
      integer   IPR_TADV
      integer   IPR_DADV
      integer   IPR_WDIF
      integer   IPR_EDIF
      integer   IPR_SDIF
      integer   IPR_NDIF
      integer   IPR_BDIF
      integer   IPR_TDIF
      integer   IPR_DDEP
      integer   IPR_WDEP
      integer   IPR_FAERO
      integer   IPR_FINAL
      integer   IPR_CONV
      integer   IPR_VOL

      parameter ( MXPADOM     =    4   )
      parameter ( MXPACEL     =    1   )
      parameter ( NPAPRC      =   24   )
      parameter ( MXCPA       =   35   )
c
      parameter ( IPR_INIT    =    1 )
      parameter ( IPR_CHEM    =    2 )
      parameter ( IPR_AEMIS   =    3 )
      parameter ( IPR_PTEMIS  =    4 )
      parameter ( IPR_PIGEMIS =    5 )
      parameter ( IPR_WADV    =    6 )
      parameter ( IPR_EADV    =    7 )
      parameter ( IPR_SADV    =    8 )
      parameter ( IPR_NADV    =    9 )
      parameter ( IPR_BADV    =   10 )
      parameter ( IPR_TADV    =   11 )
      parameter ( IPR_DADV    =   12 )
      parameter ( IPR_WDIF    =   13 )
      parameter ( IPR_EDIF    =   14 )
      parameter ( IPR_SDIF    =   15 )
      parameter ( IPR_NDIF    =   16 )
      parameter ( IPR_BDIF    =   17 )
      parameter ( IPR_TDIF    =   18 )
      parameter ( IPR_DDEP    =   19 )
      parameter ( IPR_WDEP    =   20 )
      parameter ( IPR_FAERO   =   21 )
      parameter ( IPR_FINAL   =   22 )
      parameter ( IPR_CONV    =   23 )
      parameter ( IPR_VOL     =   24 )
c
c-----------------------------------------------------------------------
c     Variables for for determing which method to 
c     use for Process Analysis:
c-----------------------------------------------------------------------
c
c     technological type
c
c     STRPA   --  string to request both
c                 Integrated Process Rates (IPR)
c                 Integrated Reaction Rates (IRR)
c     STRIPR  --  string to request Integrated Process Rates (IPR)
c     STRIRR  --  string to request Integrated Reaction Rates (IRR)
c
c     logical statements
c
c     lproca  --  logical, determine whether to do any Process Analysis
c     lipr    --  logical, if true, do Integrated Process Rates (IPR)
c     lirr    --  logical, if true, do Integrated Reaction Rates (IRR)
c     lcpacum --  logical, if true, accumulate CPA variables over hours
c
c-----------------------------------------------------------------------
c
      character*10 STRPA
      character*10 STRIPR
      character*10 STRIRR
c
      parameter ( STRPA  =  'PA        ' )
      parameter ( STRIPR =  'IPR       ' )
      parameter ( STRIRR =  'IRR       ' )
c
      logical   lproca
      logical   lipr
      logical   lirr
      logical   lcpacum
c
      common /irmb_var/ lproca, lipr, lirr, lcpacum
c
c-----------------------------------------------------------------------
c     Variables for Process Analysis data structures:
c-----------------------------------------------------------------------
c
c     npa_cels  -- total number of Process Analysis grid cells in simulation
c     nirrrxn   -- number of chemical reactions used in simulation
c     ipax      -- array to determine Process Analysis grid cell column id
c     ipay      -- array to determine Process Analysis grid cell row id
c     ipaz      -- array to determine Process Analysis layer id
c     ipagrd    -- array to determine which nesting grid is to 
c                  perform Process Analysis
c     ipanst    -- array to trace current Process Analysis grid cell back 
c                  to its nesting
c     ipadom    -- array to trace current Process Analysis grid cell 
c                  back to sub-domain
c     cipr      -- array for cumulative mass balance of each species
c                  NPAPRC  - number of PA processes
c                  MXPACEL - number of PA cells 
c                  MXSPEC  - number of species         (from Inc/camx.prm)
c
c     npastep   -- array for the number of steps in this output period
c     cirr      -- array for cumulative mass balance of each species
c                  MXPACEL - number of PA cells
c                  MXRXN   - number of reactions       (from Inc/camx.prm)
c
c     ipcacl_3d -- 3-D array containing the subdomain value for each
c                  cell in each grid
c     ipcacl_2d -- 2-D array containing the subdomain value for each
c                  cell in each grid (if any layer in the cell is used)
c     npadom    -- number of Process Analysis domain
c     i_sw      -- southwest corner column index of Process Analysis domain
c     j_sw      -- southwest corner row index of Process Analysis domain
c     i_ne      -- northeast corner column index of Process Analysis domain
c     j_ne      -- northeast corner row index of Process Analysis domain
c     b_lay     -- bottom layer of Process Analysis domain
c     t_lay     -- top layer of Process Analysis domain
c
c-----------------------------------------------------------------------
c
      integer   npa_cels
      integer   nirrrxn

      integer  ipacl_3d ( MXVEC3D )
      integer  ipacl_2d ( MXVEC2D )
      integer  ipagrd   ( MXPADOM )
      integer  ipadom   ( MXPACEL )
      integer  ipax     ( MXPACEL )
      integer  ipay     ( MXPACEL )
      integer  ipaz     ( MXPACEL )
      integer  ipanst   ( MXPACEL )
c
      integer  npadom
      integer  i_sw   ( MXPADOM )
      integer  j_sw   ( MXPADOM )
      integer  i_ne   ( MXPADOM )
      integer  j_ne   ( MXPADOM )
      integer  b_lay  ( MXPADOM )
      integer  t_lay  ( MXPADOM )
c
      real*4   cipr   ( NPAPRC,  MXPACEL, MXSPEC  )
      real*4   cirr   ( MXPACEL, MXRXN  )
      integer  npastep( MXPACEL, MXSPEC )
c
      common /irmb_dom/ npadom, i_sw, j_sw, i_ne, j_ne, b_lay, t_lay
c
      common /irmb_dim/ npa_cels, nirrrxn,
     &                  ipacl_3d, ipacl_2d, 
     &                  ipagrd, ipadom,
     &                  ipax, ipay, ipaz, ipanst, npastep
c
      common /irmb_con/ cipr, cirr
