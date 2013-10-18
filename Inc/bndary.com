c-----CAMx v4.02 030709
c 
c     BNDARY.COM contains all coarse grid boundary information
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none 
c
c
c-----------------------------------------------------------------------
c     Variables to define the computational domain:
c
c     ibeg   --  first modeled cell in each row on the West boundary
c     iend   --  last modeled cell in each row on the East boundary
c     jbeg   --  first modeled cell in each column on the South boundary
c     jend   --  last modeled cell in each column on the North boundary
c
c     NOTE:  A value of -999 indicates the entire row or column is in
c            the boundary, and is not modeled.
c-----------------------------------------------------------------------
c
      integer   ibeg(MXROWA)
      integer   iend(MXROWA)
      integer   jbeg(MXCOLA)
      integer   jend(MXCOLA)
c
      common /bndary/ ibeg, iend, jbeg, jend
c
c-----------------------------------------------------------------------
c     Variables to define concentrations on the TOP boundary:
c
c     caloft   -- invariant TOP boundary concentrations (ppm)
c-----------------------------------------------------------------------
c
      real   caloft(MXSPEC)
c
      common /aloft/ caloft
c
c-----------------------------------------------------------------------
c     Variables to define concentrations on the LATERAL boundaries:
c
c     bc -- spatially varying boundary concentrations on each lateral
c           edge (ppm)
c-----------------------------------------------------------------------
c
      real   bc(MX1D,MXLAYA,MXSPEC,4)
c
      common /lateral/ bc
