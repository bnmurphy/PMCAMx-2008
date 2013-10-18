c**** RTRACCHM
c
c-----CAMx v4.02 030709
c
c
c----------------------------------------------------------------------
c
c    Include file for parameters and data structures used in the source 
c    apportionment version of the CAMx.  This file is specific to the
c    RTRAC technology.
c    The file TRACER.COM must be included prior to including this file.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     01/23/02   -gwilson-  original development
c
c-----------------------------------------------------------------------
c   Arrays for defining species in the RTRAC technology:
c-----------------------------------------------------------------------
c
c   rtracsp   C   array of species names for the RTRAC technology
c   nrtgas    I   number of gas species used
c   lrtgas    L   .TRUE. if any gas species are used
c   nrtaero   I   number of aerosol (PM) species used
c   lrtaero   L   .TRUE. if any aerosol (PM) species are used
c   nrtrac    I   number of species total
c   nrtphot   I   number of photolysis reactions
c   nrtherm   I   number of thermal reactions
c   lsecnd    L   flag for determining if species is primary or secondary
c   ksec      I   index into regular model arrays of species which acts 
c                 primary to the secondary species
c   rtlbnd    R   lower bound for RTRAC species
c   rthlaw    R   Henry's Law constant for RTRAC gas species
c   rttfact   R   temperature dependence of Henry's Law constant
c   rtdrate   R   diffusivity ratio for RTRAC gas species
c   rtreact   R   reactivity parameter for RTRAC gas species
c   rtscale   R   surface resistance scaling factor for RTRAC gas species
c   rtdens    R   density of RTRAC PM species
c   rtlcut    R   lower size cut for RTRAC PM species
c   rtucut    R   upper size cut for RTRAC PM species
c
      character*10 rtracsp(MXTRSP)
      integer      nrtgas
      logical      lrtgas
      integer      nrtaero
      logical      lrtaero
      integer      nrtrac
      integer      nrtphot
      integer      nrtherm
      logical      lsecnd(MXTRSP)
      integer      ksec(MXTRSP)
      real         rtlbnd(MXTRSP)
      real         rthlaw(MXTRSP)
      real         rttfact(MXTRSP)
      real         rtdrate(MXTRSP)
      real         rtreact(MXTRSP)
      real         rtscale(MXTRSP)
      real         rtdens(MXTRSP)
      real         rtlcut(MXTRSP)
      real         rtucut(MXTRSP)
c
      common /rtspchr/ rtracsp
      common /rtspdat/ nrtgas, lrtgas, nrtaero, lrtaero, nrtrac, 
     &                 nrtphot, nrtherm, lsecnd, ksec, rtlbnd, rthlaw, 
     &                 rttfact, rtdrate, rtreact, rtdens, rtlcut, 
     &                 rtucut, rtscale
c
c-----------------------------------------------------------------------
c   Parameters for reactions in RTRAC technology:
c-----------------------------------------------------------------------
c
c  NAMOH    C   string for identifying the OH species
c  NAMNO3   C   string for identifying the NO3 species
c  NAMO3    C   string for identifying the O3 species
c
      character*10 NAMOH
      character*10 NAMNO3
      character*10 NAMO3
c
      parameter( NAMOH  = 'OH        ' )
      parameter( NAMNO3 = 'NO3       ' )
      parameter( NAMO3  = 'O3        ' )
c
c-----------------------------------------------------------------------
c   Variables for gas-phase reactions in RTRAC technology:
c      thermal rate constant = A * (T/Tref)**B * exp(Ea/T)
c-----------------------------------------------------------------------
c
c  jnum     I   photolysis reaction number
c  rtjfact  R   scale factor used in phololysis reaction
c  aoh      R   array of A  coefficients for OH
c  eaoh     R   array of Ea coefficients for OH
c  boh      R   array of B  coefficients for OH
c  troh     R   array of Tref coefficients for OH
c  ano3     R   array of A  coefficients for NO3
c  eano3    R   array of Ea coefficients for NO3
c  bno3     R   array of B  coefficients for NO3
c  trno3    R   array of Tref coefficients for NO3
c  ao3      R   array of A  coefficients for O3
c  eao3     R   array of Ea coefficients for O3
c  bo3      R   array of B  coefficients for O3
c  tro3     R   array of Tref coefficients for O3
c
      integer jnum(MXTRSP)
      real    rtjfact(MXTRSP)
      real    aoh(MXTRSP)
      real    eaoh(MXTRSP)
      real    boh(MXTRSP)
      real    troh(MXTRSP)
      real    ano3(MXTRSP)
      real    eano3(MXTRSP)
      real    bno3(MXTRSP)
      real    trno3(MXTRSP)
      real    ao3(MXTRSP)
      real    eao3(MXTRSP)
      real    bo3(MXTRSP)
      real    tro3(MXTRSP)
c
      common /rtrxndat/ jnum, rtjfact, aoh, eaoh, boh, troh, ano3, 
     &                  eano3, bno3, trno3, ao3, eao3, bo3, tro3
c
c-----------------------------------------------------------------------
c   Arrays for gridded fields for RTRAC species:
c-----------------------------------------------------------------------
c
c   vdeprt   R   gridded array of deposition velocities 
c   iptrvd   R   pointer into the array for each grid
c
      real    vdeprt( MXVEC2D * MXTRSP )
      integer iptrvd(MXGRID)
c
      common /rtgrd/ vdeprt, iptrvd
c
c-----------------------------------------------------------------------
c   Parameters and Arrays for receptor definitions:
c-----------------------------------------------------------------------
c
c    MXRTCEL  I  maximum number of receptor cells 
c
      integer MXRTCEL
c
      parameter( MXRTCEL = 1000 )
c
c    ircprt   I  cell index in X direction for receptor location
c    jrcprt   I  cell index in Y direction for receptor location
c    krcprt   I  cell index in Z direction for receptor location
c    idomrt   I  grid containing the receptor
c    nrcprt   I  number of receptors requested
c    rcpdcy   R  gas tracer decay rates at each receptor location
c
      integer ircprt( MXRTCEL )
      integer jrcprt( MXRTCEL )
      integer krcprt( MXRTCEL )
      integer idomrt( MXRTCEL )
      integer nrcprt
      real*4  rcpdcy( MXRTCEL, MXSPEC )
c
      common /rtrcp/ ircprt, jrcprt, krcprt, idomrt, nrcprt, rcpdcy
