c-----CAMx v4.02 030709
c
c     SOAP.COM contains common parameters and variables for SOAP
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c
c-----------------------------------------------------------------------
c     Parameters
c
c     NSOAP   -- number of CG/SOA species pairs
c-----------------------------------------------------------------------
c
c      include  'camx.prm'
c
      integer      NSOAP
      parameter  ( NSOAP = 34 )
c
      integer      pflag       ! bkoo (03/09/03)
      parameter  ( pflag = 1 ) ! 1 if there is pre-existing organic aerosol; 0 if not
c-----------------------------------------------------------------------
c     Variables that are initialized in soapdat.f
c
c     mwsoap  -- molecular weights of CG/SOA species (g/mol)
c     csat    -- saturation concentrations of CG/SOA species (ug/m3)
c     cstemp  -- temperatures corresponding to saturation concentrations
c                of CG/SOA species (K)
c     deltah  -- enthalpy of vaporization of CG/SOA species (kJ/mol)
c     flagsoap-- set to 1 if CG/SOA species forms solutions; 0 if not
c     lae3    -- true to emulate CMAQ AE3 algorithm (no evporation)
c-----------------------------------------------------------------------
c
      REAL         mwsoap(NSOAP)
      REAL         csat(NSOAP)
      REAL         cstemp(NSOAP)
      REAL         deltah(NSOAP)
      INTEGER      flagsoap(NSOAP)
      LOGICAL      lae3
c
      common /soapx/ mwsoap, csat, cstemp, deltah, flagsoap,
     &               lae3
c
c-----------------------------------------------------------------------
c  BNM
c    Variables for Enthalpy of Vaporization Calculation
c
c     ntemp -- number of temperatures for which Hvap values exist
c     ncstar -- number of volatilities for which HVAP values exist
c     dhtemp(ntemp) -- array of temperatures in the Htemp file
c     dhcstar(ncstar) -- array of volatilities in the Hcstar file
c     poadhvap(ntemp,ncstar) -- array of temp-, vol- dependent Hvap values
c    
      integer ntemp, ncstar
      parameter(ntemp = 211)
      parameter (ncstar = 109)

      real dhtemp(ntemp), dhcstar(ncstar)
      real poadhvap(ntemp,ncstar), soadhvap(ntemp,ncstar)
c
      common /Hvap/ dhtemp, dhcstar, poadhvap, soadhvap
c
c----------------------------------------------------------------------
