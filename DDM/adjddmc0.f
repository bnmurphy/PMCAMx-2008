c*** ADJDDMC0
c
      subroutine adjddmc0(ispc,igrd,icl,jcl,kcl,relc,relc0,c0trac)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine makes the adjustments to the DMM sensitivities
c     species.  The adjustments are based on the relative difference in 
c     concentrations of the regular model species due to wet dep.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       ispc    I  index into regular model species order
c       igrd    I  grid number
c       icl     I  the X grid location of current cell
c       jcl     I  the X grid location of current cell
c       kcl     I  the vertical grid location of current layer
c       relc    R  relative change in cell concentrations (multiplier)
c       relc0   R  relative change in rain concentrations (multiplier)
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
      integer   ispc
      integer   igrd
      integer   icl
      integer   jcl
      integer   kcl
      real      relc
      real      relc0
      real      c0trac(MXTRSP)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer idxcel, idx
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
c   --- loop over all sensitivities associated with this species ---
c
      do iddm=1,nddmsp
        idxspc = iptddm(ispc)+iddm-1
        idx = ipsa3d(igrd)-1+idxcel + 
     &                     ncol(igrd)*nrow(igrd)*nlay(igrd)*idxspc
        c0trac(idxspc) = (1.0-relc0) * c0trac(idxspc)
        c0trac(idxspc) = c0trac(idxspc) + relc * ptconc(idx)
        ptconc(idx) = (1.0-relc) * ptconc(idx)
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
