c**** TRACER
c
c-----CAMx v4.03 031205
c
c
c----------------------------------------------------------------------
c
c    Include file for tracer parameters and data structures used
c    in the source apportionment version of the CAMx
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     01/03/96   -gwilson-  original development
c     02/10/97   -cemery-   added SA output file root name
c     04/28/97   -gwilson-  added OPPAT and GOAT
c     05/22/97   -gwilson-  added APCA
c     11/17/99   -gwilson-  added DDM
c     01/17/02   -gwilson-  added RTRAC
c     07/31/02   -gwilson-  removed OPPAT
c     10/11/03   -gyarwood- changed I/O unit numbers to avoid conflicts
c
c***********************************************************************
c    >>>>>>>>> MXTRSP Parameter for number of tracer species.
c              The value of this parameter greatly affects the
c              amount of memory required by the model.
c
c    OSAT>>>>> For OSAT, use the following formula to determine the
c              appropriate value:
c
c       MXTRSP = 2*ntime*nreg*(ndays+1) + 4*(1+nbound+ngroup*nreg)
c
c             where:
c
c        ntime  = number of timing releases per day (usually 0)
c        nreg   = number of source regions
c        ndays  = number of days in complete simulation
c        nbound = 1 if "stratify boundary" is off, 5 otherwise
c        ngroup = number of source groupings
c
c     For a simulation without source apportionment, set MXTRSP = 1.
c
c   DDM>>>>>>> In the DDM, MXTRSP represents the number of species
c              times the number of parameters for which sensitivities
c              are calculated.  Use the following formula
c              for calculating MXTRSP:
c
c       MXTRSP = nspec * number of parameters
c              = nspec * ( nicddm + (nbcddm * nbound) +
c                           (nemddm*ngroup*nreg) )
c
c             where:
c
c       nspec  = total number of modeled species
c       nicddm = number of IC DDM species groups
c       nbcddm = number of BC DDM species groups
c       nbound = 1 if "stratify boundary" is off, 5 otherwise
c       nemddm = number of emission DDM species groups
c       nreg   = number of source regions
c       ngroup = number of source groupings
c
c   PA>>>>>>>> For Process Analysis, MXTRSP must be greater than or
c              equal to the parameter MXCPA which is set in procan.com
c              The current value for MXCPA in CAMx3 is 35
c
c***********************************************************************
c
c   MXTRSP   I   maximum number of total tracer species in simulation
c
      integer   MXTRSP
c
      parameter( MXTRSP =    1)
c
c-----------------------------------------------------------------------
c    Parameters for Source Apportionment:
c-----------------------------------------------------------------------
c
c   VEROSAT  C  character string for the version of the model
c   VERGOAT  C  character string for the version of the model
c   VERAPCA  C  character string for the version of the model
c   VERDDM   C  character string for the version of the model
c   VERPA    C  character string for the version of the model
c   VERRTRAC C  character string for the version of the model
c
      character*80 VEROSAT
      character*80 VERGOAT
      character*80 VERAPCA
      character*80 VERDDM
      character*80 VERPA
      character*80 VERRTRAC
c
      parameter(VEROSAT='Ozone Source Apportionment Technology, OSAT 030
     &709')
      parameter(VERGOAT='Goegraphic Ozone Apportionment Technology, GOAT
     & 030709')
      parameter(VERAPCA='Anthropogenic Precursor Culpability Analysis, A
     &PCA 030709')
      parameter(VERDDM='De-coupled Direct Method, DDM 030709')
      parameter(VERPA='Process Anaylsis, PA 030709')
      parameter(VERRTRAC='Reactive Tracers, RTRAC 030709')
c
c-----------------------------------------------------------------------
c    Parameters for I/O units of source grouping emissions files:
c-----------------------------------------------------------------------
c
c    MXTEMF   I  maximum number of source grouping emission files
c    IORMAP   I  unit number of the source area mapping file
c    IORRCP   I  unit number of the file with the receptor definition
c    IORIC    I  unit number of the IC file for DDM
c    IORBC    I  unit number of the BC file for DDM
c    IORTC    I  unit number of the BC file for DDM
c    IORCHM   I  unit mumber of the chemistry parameters for RTRAC
c    IORINI   I  unit number of instantaneous file used to initalize
c    IORTEM   I  base unit number of the source grouping files (surface)
c    IORTPT   I  base unit number of the source grouping files (elevated)
c    IOWCN1   I  unit number of the first output instantaneous file
c    IOWCN2   I  unit number of the second output instantaneous file
c    IOWSFC   I  unit number of the surface tracer concentrations file
c    IOWAVG   I  unit number of output average concentrations
c
c    IDXCRS   I  index of course grid in arrays
c    IDXFIN   I  index of fine grid in arrays
c
      integer   MXTEMF
      integer   IORMAP
      integer   IORRCP
      integer   IORINI
      integer   IORIC
      integer   IORBC
      integer   IORTC
      integer   IORCHM
      integer   IORTEM
      integer   IORTPT
      integer   IOWCN1
      integer   IOWCN2
      integer   IOWSFC
      integer   IOWAVG
c
      integer   IDXCRS
      integer   IDXFIN
c
      parameter( MXTEMF = 10 )
      parameter( IORMAP = 150 )
      parameter( IORRCP = 150 )
      parameter( IORIC  = 150 )
      parameter( IORBC  = 150 )
      parameter( IORTC  = 150 )
      parameter( IORCHM = 150 )
      parameter( IORINI = 150 )
      parameter( IOWAVG = 151 )
      parameter( IOWCN1 = 152 )
      parameter( IOWCN2 = 155 )
      parameter( IOWSFC = 158 )
      parameter( IORTPT = 161 )
      parameter( IORTEM = IORTPT + MXTEMF + 1 )
c
      parameter( IDXCRS = 1 )
      parameter( IDXFIN = 2 )
c
c-----------------------------------------------------------------------
c    Variables for filenames:
c-----------------------------------------------------------------------
c
c   flrtsa   C   output file root name
c   mapfil   C   filename for the source region map
c   lmapfl   L   flag to dertimine if source region map has been provided
c   rcpfil   C   filename for the receptor definition file
c   icfil    C   filename for the IC file for DDM or RTRAC
c   bcfil    C   filename for the BC file for DDM or RTRAC
c   tcfil    C   filename for the TOPCONC file for DDM or RTRAC
c   chmfil   C   filename for the chemistry parameters file for RTRAC
c   inifil   C   filename for the instantaneous file used for initializing
c                (one for each grid, zero element is for course grid)
c   temfil   C   array of filenames for source group emissions (surface)
c   ltemfl   L   flag to determine if emissions file was supplied for 
c                the source grouping (position 0 is for the regular
c                emisions file)
c   tptfil   C   array of filenames for source group emissions (elevated)
c   ltptfl   L   flag to determine if emissions file was supplied for 
c                the source grouping (position 0 is for the regular
c                emisions file)
c   cn1fil   C   filename for first output instantaneous concentrations
c   cn2fil   C   filename for second output instantaneous concentrations
c   sfcfil   C   filename for the surface tracer concentration file
c   lsfcfl   L   flag to determine if surface filename was supplied
c   avgfil   C   filename for the average tracer concentrations
c   lfirst   L   flag for firs time interval of average concs
c   verson   C   string containing the version of the model being run
c
      character*200 flrtsa
      character*200 mapfil(MXGRID)
      logical       lmapfl(MXGRID)
      character*200 rcpfil
      character*200 icfil
      character*200 bcfil
      character*200 tcfil
      character*200 chmfil
      character*200 inifil(IDXCRS:IDXFIN)
      character*200 temfil(MXGRID,MXTEMF)
      logical       ltemfl(MXGRID,0:MXTEMF)
      character*200 tptfil(MXTEMF)
      logical       ltptfl(0:MXTEMF)
      character*200 cn1fil(IDXCRS:IDXFIN)
      character*200 cn2fil(IDXCRS:IDXFIN)
      character*200 sfcfil(IDXCRS:IDXFIN)
      logical       lsfcfl(IDXCRS:IDXFIN)
      character*200 avgfil
      logical       lfirst
      character*80 verson
c     
      common /filchr/ flrtsa, mapfil, rcpfil, icfil, bcfil, tcfil, 
     &                chmfil, inifil, temfil, tptfil, cn1fil, cn2fil,
     &                sfcfil, avgfil, verson 
      common /fildat/ lmapfl, ltemfl, ltptfl, lsfcfl, lfirst
c
c-----------------------------------------------------------------------
c    Parameters for array bounds:
c-----------------------------------------------------------------------
c
c   MXSA2D   I   maximum number of entries in the 2-D concentration arrays
c   MXSA3D   I   maximum number of entries in the 3-D concentration arrays
c   MXFDDM   I   maximum number of families for DDM 
c                  (calculated from the max numver of tracer 
c                   species and max number of modeled species)
c                
c
      integer   MXSA2D
      integer   MXSA3D
      integer   MXFDDM
c
      parameter( MXSA2D = MXVEC2D*MXTRSP )
      parameter( MXSA3D = MXVEC3D*MXTRSP )
      parameter( MXFDDM = MXTRSP/MXSPEC + 1)
c
c-----------------------------------------------------------------------
c    Parameters for releasing timing tracers:
c-----------------------------------------------------------------------
c
c   TIMEMS   R   concentration (PPM) of timing tracer release emissions
c
      real      TIMEMS
c
      parameter( TIMEMS = 10.0E-07 )
c
c-----------------------------------------------------------------------
c    Parameters for array indices:
c-----------------------------------------------------------------------
c
c    ITRNOX   I   base index of NOx emissions tracer
c    ITRVOC   I   base index of VOC emissions tracer
c    ITRONX   I   base index of NOx ozone reactivity tracer
c    ITROVC   I   base index of VOC ozone reactivity tracer
c    ITRTNX   I   base index of NOx timing tracers     
c    ITRTVC   I   base index of VOC timing tracers     
c                 (NOTE:  New timing tracer species will be created 
c                         throughout the simulation, so these species
c                         will appear at the end of the list)
c
      integer   ITRNOX 
      integer   ITRVOC 
      integer   ITRONX 
      integer   ITROVC 
      integer   ITRTNX 
      integer   ITRTVC 
c
      parameter( ITRNOX = 1 )
      parameter( ITRVOC = 2 )
      parameter( ITRONX = 3 )
      parameter( ITROVC = 4 )
      parameter( ITRTNX = 1 )
      parameter( ITRTVC = 2 )
c
c-----------------------------------------------------------------------
c    Variables for gridded tracer emissions data:
c-----------------------------------------------------------------------
c
c   saemis    R   gridded array of tracer emissions from the surface
c                 emissions files
c   sapnts    R   array of tracer emissions for each point source
c   xlocpt    R   location of the point source in X direction
c   ylocpt    R   location of the point source in Y direction
c   ntrtim    I   number of times per day that a timing release is done
c   lreles    L   flag for determing if a new timing tracer should
c                 be released
c   nreles    I   number of the current timing release
c   ipigsp    I   index into tracer species array of the region and
c                 group where each PiG puff originated
c   lpigsa    L   flag to determine if source is a PiG source
c                   
      integer   nreles
      integer   ntrtim
      integer   ipigsp(MXPTSRC)
      real      saemis(MXSA2D)
      real      sapnts(MXPTSRC,MXTRSP)
      real      xlocpt(MXPTSRC)
      real      ylocpt(MXPTSRC)
      logical   lreles
      logical   lpigsa(MXPTSRC)
c
      common /empdat/  saemis, sapnts, ipigsp, xlocpt, ylocpt, lreles, 
     &                 nreles, ntrtim, lpigsa
c
c-----------------------------------------------------------------------
c    Parameters for gridded tracer concentration data:
c-----------------------------------------------------------------------
c
c   BNDLPT    R   lower bound for tracer concentrations
c
      real      BNDLPT
c
      parameter( BNDLPT = 1.0E-16 )
c
c-----------------------------------------------------------------------
c    Variables for gridded tracer concentration data:
c-----------------------------------------------------------------------
c
c   ptconc    R   gridded array of tracer concentrations by layer(OSAT)
c                 or gridded array of sensitivities by layer (DDM)
c   ptavrg    R   average tracer concentrations at surface (used for
c                 hourly peak receptor)
c   modavg    R   average regular model species at species 
c   ptloft    R   tracer concentrations at the top of the model
c   ptvdep    R   diffusion velocities for tracer species
c
      real   ptconc(MXSA3D)
      real   ptavrg(MXSA2D)
      real   modavg(MXCOL1*MXROW1*MXSPEC)
      real   ptloft(MXTRSP)
      real   ptvdep(MXCOLA*MXROWA*MXTRSP)
c
      common /ptgdat/ ptconc, ptavrg, ptloft, modavg, ptvdep
c
c-----------------------------------------------------------------------
c    Variables for tracer names:
c-----------------------------------------------------------------------
c
c   ptname   C    name of passive tracer species (imbedded into this name
c                 will be the type of tracer, the source region from 
c                 where it originated and the source grouping emissions
c                 file from which it originated)
c
c   nsaspc   I    number of tracer species currently in the simulation
c   npttim   I    number of timing tracer species (including future timings)
c   ntotsp   I    number of total tracer species (including future timings)(OSAT)
c                 or number of species times number of parameters (DDM)
c   iptnox   I    index into tracer list of the beginning of the NOx
c                 concentration tracers
c   iemnox   I    index into tracer list of the begginning of the NOx
c                 emissions tracers
c   iptvoc   I    index into tracer list of the beginning of the VOC
c                 emissions tracers
c   iemvoc   I    index into tracer list of the begginning of the VOC
c                 emissions tracers
c   ipto3n   I    index into tracer list of the beginning of the NOx 
c                  reaction tracers
c   ipto3v   I    index into tracer list of the beginning of the VOC
c                  reaction tracers
c   ipttim   I    index into tracer list of the beginning of the
c                 timing tracers
c   iemtim   I    index into tracer list of the begginning of the 
c                 timing emissions tracers
c   ipto3n   I    index into tracer list of the beginning of the NOx 
c   nbdic    I    number of boundary/initial conditions tracer species
c                 (depends on wether the boundary is stratified)
c   lsamap   I    dummy mapping of arrays species (lsamap(i) = i)
c
      character*10 ptname(MXTRSP)
      integer      nsaspc
      integer      npttim
      integer      ntotsp
      integer      iptnox
      integer      iemnox
      integer      iptvoc
      integer      iemvoc
      integer      ipto3n
      integer      ipto3v
      integer      ipttim 
      integer      iemtim 
      integer      nbdic
      integer      lsamap(MXTRSP)
c
      common /ptnchr/ ptname
      common /ptndat/ nsaspc, npttim, ntotsp, iptnox, iemnox, iptvoc, 
     &                iemvoc, ipto3n, ipto3v, ipttim, iemtim, nbdic,
     &                lsamap
c
c-----------------------------------------------------------------------
c    Parameters for user options and flags:
c-----------------------------------------------------------------------
c
c   OSAT   C   code for OSAT technology
c   GOAT   C   code for GOAT technology
c   APCA   C   code for APCA technology
c   DDM    C   code for DDM technology
c   RTRAC  C   code for RTRAC technology
c
      character*10 OSAT
      character*10 GOAT
      character*10 APCA
      character*10 DDM
      character*10 RTRAC
c
      parameter( OSAT  = 'OSAT      ')
      parameter( GOAT  = 'GOAT      ')
      parameter( APCA  = 'APCA      ')
      parameter( DDM   = 'DDM       ')
      parameter( RTRAC = 'RTRAC     ')
c
c-----------------------------------------------------------------------
c    Variables for user options and flags:
c-----------------------------------------------------------------------
c
c   ltrace   L   flag for determining if the passive tracer algorithm 
c                should be used
c   lddm     L   flag for determining if DDM is being used 
c   lrestrt  L   flag for determining if the simulation is a first day
c   leftovr  L   flag for determining if the left-over group should be 
c                used
c   lbndry   L   flag to determine if the boundary conditions should
c                be stratified by edge
c   ngroup   I   number of source groupings for emissions
c   tectyp   C   flag for determining which type of technology will
c                be performed.
c   loutsa   L   flag for determining if the species should
c                be output to average file
c
      integer      ngroup
      integer      nchar
      logical      ltrace
      logical      lddm
      logical      lrestrt
      logical      leftovr
      logical      lbndry
      character*10 tectyp
      logical      loutsa(MXTRSP)
c
      common /usrchr/ tectyp
      common /usrdat/ ngroup, nchar, ltrace, lrestrt, leftovr, lbndry,
     &                lddm, loutsa
c
c-----------------------------------------------------------------------
c    Variables for region mapping:
c-----------------------------------------------------------------------
c
c   nregin   I   number of source regions 
c   nxcell   I   number of cells in X-direction
c   nycell   I   number of cells in X-direction
c   igrmap   I   grid that maps cell to source region
c
      integer   nregin
      integer   nxcell(MXGRID)
      integer   nycell(MXGRID)
      integer   igrmap(MXGRID,MXCOLA,MXROWA)
c
      common /mapdat/ nregin, nxcell, nycell, igrmap
c
c-----------------------------------------------------------------------
c    Parmeters for receptor variables:
c-----------------------------------------------------------------------
c
c   MXRECP   I   maximum number of receptor locations
c   MXCELR   I   maximum number of cells in CELL AVERAGE type of receptor
c   CDPNT    C   string for indicating POINT type of recptor
c   IDPNT    I   id code for POINT type of receptor
c   CDCEL    C   string for indicating SINGLE CELL type of recptor
c   IDCEL    I   id code for SINGLE CELL type of receptor
c   CDAVG    C   string for indicating CELL AVERAGE type of recptor
c   IDAVG    I   id code for CELL AVERAGE type of receptor
c   CDWAL    I   string indicating WALL OF CELLS type of receptors
c   IDWAL    I   id code for WALL OF CELLS type of receptors
c
      character*20 CDPNT
      character*20 CDCEL
      character*20 CDAVG
      character*20 CDWAL
      integer      MXRECP
      integer      MXCELR
      integer      IDPNT
      integer      IDCEL
      integer      IDAVG
      integer      IDWAL
c
      parameter( CDPNT  = 'POINT          ' )
      parameter( CDCEL  = 'SINGLE CELL    ' )
      parameter( CDAVG  = 'CELL AVERAGE   ' )
      parameter( CDWAL  = 'WALL OF CELLS  ' )
      parameter( MXRECP = 100 )
      parameter( MXCELR = 70  )
      parameter( IDPNT  = 1   )
      parameter( IDCEL  = 2   ) 
      parameter( IDAVG  = 3   ) 
      parameter( IDWAL  = 4   ) 
c
c-----------------------------------------------------------------------
c    Variables for receptor data:
c-----------------------------------------------------------------------
c
c   rcpnam    C   names of the receptors
c   nrecep    I   number of receptors
c   idrcp     I   array of id codes for each receptor
c   igrdrcp   I   grid index for each receptor
c   irecep    I   array of I-cell locations for CELL type of receptors
c   jrecep    I   array of I-cell locations for CELL type of receptors
c   nclrcp    I   number of cells to average for CELL type of receptors
c   recepx    R   array of X-ccordinates for POINT type of receptors
c   recepy    R   array of Y-ccordinates for POINT type of receptors
c   conrcp    R   array of tracer surface concentrations at each receptor
c   rcpvoc    R   array for regular model VOC concentration at receptor
c   rcpnox    R   array for regular model NOx concentration at receptor
c   rcpo3     R   array for regular model O3 concentration at receptor
c   ipekcl    I   index into gridded array of the hourly peak cell
c   iwalbg    I   beginning column for WALL OF CELLS receptor
c   iwalnd    I   ending column for WALL OF CELLS receptor
c   jwalbg    I   beginning row for WALL OF CELLS receptor
c   jwalnd    I   ending row for WALL OF CELLS receptor
c   kwalbg    I   beginning layer for WALL OF CELLS receptor
c   kwalnd    I   ending layer for WALL OF CELLS receptor
c   lwalls    L   flag to determine if any WALL OF CELLS receptors 
c                 were specified
c
      character*10 rcpnam(MXRECP)
      integer      nrecep
      integer      idrcp(MXRECP)
      integer      igrdrcp(MXRECP)
      integer      irecep(MXRECP,MXCELR)
      integer      jrecep(MXRECP,MXCELR)
      integer      nclrcp(MXRECP)
      integer      ipekcl(2)
      integer      iwalbg(MXRECP)
      integer      iwalnd(MXRECP)
      integer      jwalbg(MXRECP)
      integer      jwalnd(MXRECP)
      integer      kwalbg(MXRECP)
      integer      kwalnd(MXRECP)
      real         recepx(MXRECP)
      real         recepy(MXRECP)
      real         conrcp(MXTRSP,MXRECP)
      real         rcpvoc(MXRECP)
      real         rcpnox(MXRECP)
      real         rcpo3(MXRECP)
      logical      lwalls
c
      common /rcpchr/ rcpnam
      common /rcpdat/ nrecep, idrcp, igrdrcp, irecep, jrecep, nclrcp, 
     &                ipekcl, recepx, recepy, conrcp, rcpvoc, rcpnox, 
     &                rcpo3,  iwalbg, iwalnd, jwalbg, jwalnd, kwalbg, 
     &                kwalnd, lwalls
c
c-----------------------------------------------------------------------
c    Parameters for boundary conditions:
c-----------------------------------------------------------------------
c
c   IDXBWS   I   index of the WEST boundary in arrays
c   IDXBES   I   index of the EAST boundary in arrays
c   IDXBST   I   index of the SOUTH boundary in arrays
c   IDXBNT   I   index of the NORTH boundary in arrays
c   IDXBTP   I   index of the TOP boundary in arrays
c
      integer   IDXBWS
      integer   IDXBES
      integer   IDXBST
      integer   IDXBNT
      integer   IDXBTP
c
      parameter( IDXBWS = 1 )
      parameter( IDXBES = 2 )
      parameter( IDXBST = 3 )
      parameter( IDXBNT = 4 )
      parameter( IDXBTP = 5 )
c
c
c-----------------------------------------------------------------------
c    Variables for species order and species flags:
c-----------------------------------------------------------------------
c
c  NOTE:  In all arrays, position zero is for the regular emissions
c         files.
c
c   lvocsp   L   flag to determine if species is VOC species
c   lnoxsp   L   flag to determine if species is NOx species
c   lo3sp    L   flag to determine if species is O3 species
c   crbnum   R   carbon number of each species 
c   idxems   I   index of each species in each emissions file into
c                species list arrays
c   idxpts   I   index of each species in each point source file
c                into species list arrays
c   nspcem   I   number of species in each emissions file
c   nspcpt   I   number of species in each point source file
c   reacrt   R   the reactivity fraction of each species
c   vocwt    R   the VOC reactivity weighting factor for each tracer 
c                species (really only needed for VOC species, but easier
c                to code this way)
c 
      integer   idxems(MXGRID,0:MXTEMF,MXSPEC)
      integer   nspcem(MXGRID,0:MXTEMF)
      integer   idxpts(0:MXTEMF,MXSPEC)
      integer   nspcpt(0:MXTEMF)
      logical   lvocsp(MXSPEC)
      logical   lnoxsp(MXSPEC)
      logical   lo3sp(MXSPEC)
      real      crbnum(MXSPEC)
      real      reacrt(MXSPEC)
      real      vocwt(MXTRSP)
c
      common /spcdat/ idxems, nspcem, idxpts, nspcpt, lvocsp, lnoxsp, 
     &                lo3sp, crbnum, reacrt, vocwt
c
c-----------------------------------------------------------------------
c  Parameters for DDM flags:
c-----------------------------------------------------------------------
c
c    NAMVOC  C   character string for name of VOC species
c    NAMNOX  C   character string for name of NOX species
c    NAMALL  C   character string for name of ALL species
c
c    IDVOC   I   species ID for VOC species
c    IDNOX   I   species ID for NOX species
c    IDALL   I   species ID for ALL species
c
      character*10 NAMVOC
      character*10 NAMNOX
      character*10 NAMALL
      integer      IDVOC
      integer      IDNOX
      integer      IDALL
c
      parameter( NAMVOC = 'VOC       ')
      parameter( NAMNOX = 'NOX       ')
      parameter( NAMALL = 'ALL       ')
      parameter( IDVOC  = -1 )
      parameter( IDNOX  = -2 )
      parameter( IDALL  = -3 )
c
c-----------------------------------------------------------------------
c  Tracer species map for DDM
c-----------------------------------------------------------------------
c
c   icddmsp  C  species names for initial conditions treated by DDM
c   bcddmsp  C  species names for boundary condition treated by DDM
c   emddmsp  C  species names for emissions treated by DDM
c   nicddm   I  number of initial conditions groups in DDM
c   nbcddm   I  number of boundary conditions groups in DDM
c   nemddm   I  number of emissions groups in DDM
c 
      character*10 icddmsp(MXSPEC)
      character*10 bcddmsp(MXSPEC)
      character*10 emddmsp(MXSPEC)
      integer      nicddm
      integer      nbcddm
      integer      nemddm
c
      common /ddmchr/ icddmsp, bcddmsp, emddmsp
      common /ddmdat/ nicddm, nbcddm, nemddm
c
c-----------------------------------------------------------------------
c  Species list and pointers into arrays for DDM species:
c-----------------------------------------------------------------------
c
c   ptlong   C  long names of the DDM species
c   iptddm   I  index into the gridded arrays of the DDM families
c   nddmsp   I  total number of DDM parameters
c
      character*14 ptlong(MXTRSP) 
      integer      iptddm(MXSPEC)
      integer*4    nddmsp
c
      common /dspchr/ ptlong
      common /dspdat/ iptddm, nddmsp
c
c-----------------------------------------------------------------------
c  Tracer species map for netCDF I/O
c-----------------------------------------------------------------------
c
      common /tracemap/ linsamap(MXTRSP,MXGRID),lavsamap(MXTRSP,MXGRID),
     &                  lptsamap(MXSPEC,0:MXTEMF),
     &                  larsamap(MXSPEC,0:MXTEMF,MXGRID)
      common /tracefil/ ncemsa(0:MXTEMF,MXGRID),ncptsa(0:MXTEMF),
     &                  ncinsa(MXGRID),ncsfsa(MXGRID),ncrstsa(MXGRID)
