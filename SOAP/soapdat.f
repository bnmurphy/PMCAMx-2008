      Block Data soapdat
c
c-----------------------------------------------------------------------
c     Initialize common block variables for SOAP
c
c     mwsoap   - molecular weights of CG/SOA species (g/mol)
c     csat     - saturation concentrations of SOA species (ug/m3)
c     cstemp   - temperatures corresponding to saturation concentrations
c                of CG/SOA species (K)
c     deltah   - enthalpy of vaporization of CG/SOA species (J/mol)
c     flagsoap - 1 if SOA species forms solutions; 0 if not
c     lae3     - flag to set emulation of the CMAQ AE3 algorithm
c                which allows no evaporation of SOA (recommend false)
c-----------------------------------------------------------------------
c
      include 'soap.com'


c--------------------------------------
c BNM
c   Aerosol Map Table
c
c     / APO1,  APO2,  APO3,  APO4,
c       APO5,  APO6,  APO7,  APO8,
c       APO9,  APO10, AOO1,  AOO2,
c       AOO3,  AOO4,  AOO5,  AOO6,
c       AOO7,  AOO8,  AOO9,  AOO10,
c       ABS1,  ABS2,  ABS3,  ABS4,
c       AAS1,  AAS2,  AAS3,  AAS4   /
c
c--------------------------------------

c
      data mwsoap   /250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0,
     $               180.0, 180.0, 180.0, 180.0,
     $               150.0, 150.0, 150.0, 150.0 /
      data csat     /0.01, 0.1, 1., 10., 100., 1000., 10000., 100000.,
     &               1000000., 0.0,
     &               0.01, 0.1, 1., 10., 100., 1000., 10000., 100000., 
     &               1000000., 0.0,
     &               1., 10., 100., 1000.,
     &               1., 10., 100., 1000. /
      data cstemp   /300., 300., 300., 300., 
     &               300., 300., 300., 300., 
     &               300., 300., 300., 300.,
     &               300., 300., 300., 300., 
     &               300., 300., 300., 300., 
     &               300., 300., 300., 300., 
     &               300., 300., 300., 300. /
      data deltah   /112000., 106000., 100000., 94000., 88000., 82000., 
     &               76000., 70000., 64000., 112000, 112000, 106000,
     &               100000, 94000., 88000., 82000., 76000., 70000., 
     &               64000., 112000.,
     &               30000., 30000., 30000., 30000., 
     &               30000., 30000., 30000., 30000. /
      data flagsoap /28*1/
      data lae3     /.false./
c
      end
