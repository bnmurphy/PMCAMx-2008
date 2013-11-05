

----CAMx v6.00 130506

**HADVPPM** performs advection using the one-dimensional implementation
of the piecewise parabolic method of Colella and Woodward (1984).
A piecewise continuous parabola is used as the intepolation polynomial.
The slope of the parabola at cell edges is computed from a cumulative
function of the advected quantity.  These slopes are further modified
so that the interpolation function is monotone.

This version based on CMAQ HPPM.F, v 1.1.1.1 9/14/98, written by
M.T. Odman (10/5/93), NCSC.  This version assumes constant grid cell
size.

The following definitions are used:  

      |-----------> Positive direction  

|Boundary|<-----------------Domain----------------->|Boundary| 

| CON(1) | CON(2) |  ...  | CON(I) |  ...  |CON(N-1)| CON(N) |  

VEL(1)-->|     VEL(I-1)-->|        |-->VEL(I)       |-->VEL(N-1)  

FP(1)-->|      FP(I-1)-->|        |-->FP(I)        |-->FP(N-1)  

FM(2)<--|        FM(I)<--|        |<--FM(I+1)      |<--FM(N)  

                    -->|   DX   |<-- 

Copyright 1996 - 2013
ENVIRON International Corporation

Modifications:  
5/17/00   small modification to flux1,2 to improve mass accounting
10/13/03   area weighting applied to winds rather than conc
11/04/03   only single FLXARR vector passed back for use in
          ZRATES and probling tools
 9/1/09   Revised area weighting, removed from wind vector

Input arguments:  
nn                  Number of cells  
dt                  Time step (s)  
dx                  Length of cell (m)
con                 Concentration vector (umol/m3) 
vel                 Wind speed vector (m/s)  
area                Cell area adjustment vector (1/m2) 
areav               Interfacial area adjustment vector (m2) 

Output arguments:  
con                 Concentration vector (umol/m3)  
flxarr              Conc change from interfacial mass flux (umol/m3)

Routines called:  
none  

Called by:  
XYADVEC  
ZRATES

--------------------------------------------------------------------------------

**XYADVEC** drives 2-D advection of concentrations.  This version also
performs the 2-D advection on DDM sensitivies, if DDM is implemented

Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
ENVIRON International Corporation

Modifications:
4/23/99   Coarse grid outflow BC's set equal to concentration in
      the first inner computational cells
4/26/99   Added Piecewise Parabolic Method for horizontal advection
10/30/01   Revised map scale factor application to OSAT fluxes to be
      more consistent with how the fluxes are used.
12/07/01   added instructions for OMP and rearranged some loops
      to facilitate parallel processing
01/22/02   now only calls OSAT routines if not doing RTRAC
01/30/02   Added code for RTRAC probing tool
4/10/03   X/Y advection now uses layer-dependent timestep

Input arguments:
igrd              grid index
xyordr            order of x & y advection
ncol              number of columns
nrow              number of rows
nlay              number of layers
nspc              number of species
nsen              number of species times number of DDM parameters
nadv              number of sub-steps per timestep
deltat            time step (s)
dx                cell size in x-direction (m)
dy                cell size in y-direction (m)
windu             wind speed in x-direction (m/s)
windv             wind speed in y-direction (m/s)
depth             layer depth (m)
mapscl            map-scale factor at cell centroids
conc              species concentrations (umol/m3)
sens              sensitivity coefficient (umol/m3/parameter unit)
tarray2           CPU timing arguments (s)
isaptr            pointer into tracer conc array for this grid
ipa_cel           gridded array to identify if cell is
              in a IPRM sub-domain

Output arguments:
conc              species concentrations (umol/m3)
fluxes            fluxes across the boundaries (umol)
sens              sensitivity coefficient (umol/m3/parameter unit)

Routines Called:
HADVBOT
HADVPPM
BOTTDDM

Called by:
EMISTRNS
