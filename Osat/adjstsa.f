c*** ADJSTSA
c
      subroutine adjstsa(igrid,icell,jcell,kcell,delo3,delnox,delvoc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine makes the adjustments to the tracer 
c     species.  The adjustments are based on the differences in 
c     concentrations of the regular model species before and after
c     PiG slaughter is performed.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       igrid   I  grid number
c       icell   I  the X grid location of current cell
c       jcell   I  the X grid location of current cell
c       kcell   I  the vertical grid location of current layer
c       delo3   R  change in O3 concentrations 
c       delnox  R  change in NO concentrations 
c       delvoc  R  change in VOC concentrations (carbon weighted sum 
c                  of VOC species)
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c       12/08/96  --gwilson--   Original development
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
      integer   igrid
      integer   icell
      integer   jcell
      integer   kcell
      real      delo3
      real      delnox
      real      delvoc
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   idxcel, idx, i
      real      sumnox, sumvoc, sumo3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  icell + ncol(igrid)*(jcell-1) + 
     &                      ncol(igrid)*nrow(igrid)*(kcell-1)
c
c   --- loop over the NOx tracer species, calculate total NOx tracer ---
c
      sumnox = 0.
      do 10 i=iptnox,iptvoc-1
         idx = ipsa3d(igrid)-1+idxcel + 
     &                       ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         sumnox = sumnox + ptconc(idx)
   10 continue
c
c   ---- apply the wet depostition to the tracer species, apportioning 
c        by contribution to total ----
c
      if( sumnox .GT. 0. .AND. delnox .NE. 0. ) then
         do 20 i=iptnox,iptvoc-1
             idx = ipsa3d(igrid)-1+idxcel + 
     &                        ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
             ptconc(idx) = ptconc(idx) + delnox * ptconc(idx) / sumnox
             ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   20    continue
      endif
c
c   --- loop over the VOC tracer species, calculate total VOC tracer ---
c
      sumvoc = 0.
      do 30 i=iptvoc,ipto3n-1
         idx = ipsa3d(igrid)-1+idxcel + 
     &                        ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         sumvoc = sumvoc + ptconc(idx) * vocwt(i)
   30 continue
c
c   ---- apply the wet deposition to each tracer species, apportioning 
c        by contribution to total ----
c
      if( sumvoc .GT. 0. .AND. delvoc .NE. 0. ) then
         do 40 i=iptvoc,ipto3n-1
             idx = ipsa3d(igrid)-1+idxcel + 
     &                        ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
             ptconc(idx) = ptconc(idx) + delvoc * ptconc(idx) * 
     &                                               vocwt(i) / sumvoc
             ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   40    continue
      endif
c
c   --- loop over the O3 tracer species, calculate total O3 tracer ---
c
      sumo3 = 0.
      do 50 i=ipto3n,ipttim-1
         idx = ipsa3d(igrid)-1+idxcel + 
     &                     ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         sumo3 = sumo3 + ptconc(idx)
   50 continue
c
c   ---- apply the wet deposition to each tracer species, apportioning
c        by contribution to total ----
c
      if( sumo3 .GT. 0. .AND. delo3 .NE. 0. ) then
         do 60 i=ipto3n,ipttim-1
             idx = ipsa3d(igrid)-1+idxcel + 
     &                        ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
             ptconc(idx) = ptconc(idx) + delo3 * ptconc(idx) / sumo3
             ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   60    continue
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
