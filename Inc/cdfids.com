c
c-----CAMx v4.02 030709
c
c     CDFIDS.COM contains NetCDF file ID's for met data
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c-----------------------------------------------------------------------
c     Variables that identify the parameter list in the NetCDF met file:
c
c     idwindu  --  pointer for U-component of the wind
c     idwindv  --  pointer for V-component of the wind
c     idtempk  --  pointer for temperature
c     idrkv    --  pointer for darkness flag
c     idpress  --  pointer for pressure field
c     idwater  --  pointer for water vapor concentration field
c     idrain   --  pointer for rainfall field
c     idcloud  --  pointer for cloud cover field
c     idprss0  --  pointer for initial pressure field
c     idtopo   --  pointer for topography field
c     idhght   --  pointer for layer heights field
c     iwater   --  flag to determine if water vapor is provided as ppm
c                  or mixing ratio
c-----------------------------------------------------------------------
c
      integer   idwindu(MXGRID)
      integer   idwindv(MXGRID)
      integer   idtempk(MXGRID)
      integer   idrkv(MXGRID)
      integer   idpress(MXGRID)
      integer   idwater(MXGRID)
      integer   idrain(MXGRID)
      integer   idcloud(MXGRID)
      integer   idprss0(MXGRID)
      integer   idtopo(MXGRID)
      integer   idhght(MXGRID)
      integer   iwater
c
      common /varids/ idwindu, idwindv, idtempk, idrkv, idpress, 
     &                idwater, idrain, idcloud,idprss0, idtopo, idhght, 
     &                iwater

