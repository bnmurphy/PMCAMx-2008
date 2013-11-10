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
c BNM - Aerosol Map Table
c
c     / APO1,  APO2,  APO3,  APO4,
c       APO5,  APO6,  APO7,  APO8,
c       AOO1,  AOO2,  AOO3,  AOO4,  
c	AOO5,  AOO6,  AOO7,  AOO8,
c       ABS1,  ABS2,  ABS3,  ABS4,  ABS5,
c       AAS1,  AAS2,  AAS3,  AAS4   AAS5,
c       ANS1,  ANS2,  ANS3,  ANS4,  
c	ANS5,  ANS6,  ANS7,  ANS8  /
c
c--------------------------------------

c
      data mwsoap   /250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0,
     $               180.0, 180.0, 180.0, 180.0, 180.0,
     $               150.0, 150.0, 150.0, 150.0, 150.0, 
     $               250.0, 250.0, 250.0, 250.0,
     $               250.0, 250.0, 250.0, 250.0 /
      data csat     /0.1, 1., 10., 100., 1000., 10000., 100000.,1000000.,
     &               0.1, 1., 10., 100., 1000., 10000., 100000.,1000000.,
     &               0.1, 1., 10., 100., 1000.,
     &               0.1, 1., 10., 100., 1000.,
     &               0.1, 1., 10., 100., 1000., 10000., 100000.,1000000. /
      data cstemp   /300., 300., 300., 300., 
     &               300., 300., 300., 300., 
     &               300., 300., 300., 300.,
     &               300., 300., 300., 300., 
     &               300., 300., 300., 300., 300., 
     &               300., 300., 300., 300., 300.,
     &               300., 300., 300., 300.,
     &               300., 300., 300., 300. /
      data deltah   /106000., 100000., 94000., 88000., 82000., 
     &               76000., 70000., 64000., 
     &		     106000,  100000,  94000., 88000., 82000., 
     &		     76000., 70000., 64000.,
     &               30000., 30000., 30000., 30000., 30000.,
     &               30000., 30000., 30000., 30000., 30000.,
     &		     106000,  100000,  94000., 88000., 82000., 
     &		     76000., 70000., 64000. /
   
      data flagsoap /1, 1, 1, 1, 1, 1, 1, 1,
     &               1, 1, 1, 1, 1, 1, 1, 1,
     &               1, 1, 1, 1, 1,
     &               1, 1, 1, 1, 1,
     &               1, 1, 1, 1, 1, 1, 1, 1/

c      data flagsoap / 34*1 /



      data lae3     /.false./

      end
