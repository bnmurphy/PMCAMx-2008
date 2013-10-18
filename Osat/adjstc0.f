c*** ADJSTC0
c
      subroutine adjstc0(igrd,icl,jcl,kcl,delo3,delnox,delvoc,
     &                conco3,concnox,concvoc,c0o3tr,c0noxtr,c0voctr,
     &                                      c0o3,c0nox,c0voc,c0trac)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine makes the adjustments to the tracer 
c     species.  The adjustments are based on the relative difference in 
c     concentrations of the regular model species due to wet dep.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       igrd    I  grid number
c       icl     I  the X grid location of current cell
c       jcl     I  the X grid location of current cell
c       kcl     I  the vertical grid location of current layer
c       delo3   I  change in cell concentration for O3 
c       delnox  I  change in cell concentration for NOx
c       delvoc  I  change in cell concentration for VOC
c       conco3  I  current cell concentration for O3 
c       concnox I  current cell concentration for NOx
c       concvoc I  current cell concentration for VOC
c       c0o3tr  I  change in rain concentration for O3 
c       c0noxtr I  change in rain concentration for NOx
c       c0voctr I  change in rain concentration for VOC
c       c0o3    I  current rain concentration for O3 
c       c0nox   I  current rain concentration for NOx
c       c0voc   I  current rain concentration for VOC
c     Outputs 
c       c0trac  R  updated rain concentrations for tracer species
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c       04/08/03  --gwilson--   Original development
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrd
      integer icl
      integer jcl
      integer kcl
      real    delo3
      real    delnox
      real    delvoc
      real    conco3
      real    concnox
      real    concvoc
      real    c0o3tr
      real    c0noxtr
      real    c0voctr
      real    c0o3
      real    c0nox
      real    c0voc
      real    c0trac(MXTRSP)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer idxcel, idx, i
      real    relc, relc0
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  icl + ncol(igrd)*(jcl-1) + 
     &                      ncol(igrd)*nrow(igrd)*(kcl-1)
c
c   --- if this is a NOx, loop over the NOx tracer species ----
c
      do i=iptnox,iptvoc-1
          relc = 0.
          if( concnox .NE. 0. ) relc = delnox / concnox
          relc0 = 0.
          if( c0nox .NE. 0. ) relc0 = c0noxtr / c0nox
          idx = ipsa3d(igrd)-1+idxcel + 
     &                        ncol(igrd)*nrow(igrd)*nlay(igrd)*(i-1)
          c0trac(i) = (1.0-relc0) * c0trac(i)
          c0trac(i) = c0trac(i) + relc * ptconc(idx)
          ptconc(idx) = (1.0-relc) * ptconc(idx)
          ptconc(idx) = MAX(BNDLPT,ptconc(idx))
      enddo
c
c   --- if this is a VOC, loop over the VOC species ----
c
      do i=iptvoc,ipto3n-1
          relc = 0.
          if( concvoc .NE. 0. ) relc = delvoc / concvoc
          relc0 = 0.
          if( c0voc .NE. 0. ) relc0 = c0voctr / c0voc
          idx = ipsa3d(igrd)-1+idxcel + 
     &                        ncol(igrd)*nrow(igrd)*nlay(igrd)*(i-1)
          c0trac(i) = (1.0-relc0) * c0trac(i)
          c0trac(i) = c0trac(i) + relc * ptconc(idx)
          ptconc(idx) = (1.0-relc) * ptconc(idx)
          ptconc(idx) = MAX(BNDLPT,ptconc(idx))
      enddo
c
c   --- if this is ozone, loop over the O3 tracers ---
c
      do i=ipto3n,ipttim-1
          relc = 0.
          if( conco3 .NE. 0. ) relc = delo3 / conco3
          relc0 = 0.
          if( c0o3 .NE. 0. ) relc0 = c0o3tr / c0o3
          idx = ipsa3d(igrd)-1+idxcel + 
     &                        ncol(igrd)*nrow(igrd)*nlay(igrd)*(i-1)
          c0trac(i) = (1.0-relc0) * c0trac(i)
          c0trac(i) = c0trac(i) + relc * ptconc(idx)
          ptconc(idx) = (1.0-relc) * ptconc(idx)
          ptconc(idx) = MAX(BNDLPT,ptconc(idx))
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
