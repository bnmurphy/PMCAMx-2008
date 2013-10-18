c*** PIGSA
c
      subroutine pigsa(ncols,nrows,nlays,nspc,conc,
     &                                icell,jcell,kcell,idxpig,delnox)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine puts the NOx mass from the PiG into the tracer
c     concentration array.  The mass is added to the tracer species
c     from the region/group from which the PiG originated.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       ncols   I  number of columns
c       nrows   I  number of rows
c       nlays   I  number of layers
c       nspc    I  number of species
c       conc    R  gridded array of concentrations
c       icell   I  the X grid location of current cell
c       jcell   I  the X grid location of current cell
c       kcell   I  the vertical grid location of current layer
c       idxpig  I  the index of the puff in PiG arrays
c       delnox  R  change in NO concentrations 
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c       12/08/96  --gwilson--   Original development
c       11/10/97  --gwilson--   Removed unused argument: IGRID
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'pigsty.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   ncols
      integer   nrows
      integer   nlays
      integer   nspc
      real      conc(ncols,nrows,nlays,nspc)
      integer   icell
      integer   jcell
      integer   kcell
      integer   idxpig
      real      delnox
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   idx, ispc
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- get the id of the point source ----
c
      idx = idpig(idxpig)
c
c  ---- get the species offset for this source ----
c
      ispc = iemnox - 1 + ipigsp(idx)
c
c  --- put the mass into the gridded concentrations array ---
c
      conc(icell,jcell,kcell,ispc) = 
     &                          conc(icell,jcell,kcell,ispc) + delnox
      conc(icell,jcell,kcell,ispc) = 
     &                     MAX( BNDLPT,conc(icell,jcell,kcell,ispc) )
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
