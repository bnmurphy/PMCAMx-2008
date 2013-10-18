c*** GOATSA
c
      subroutine goatsa(ldrk,igrid,icell,jcell,kcell,
     &                              delo3,delno,delno2,delvoc,dtime)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine makes the "chemistry" adjustments to the tracer 
c     species.  The adjustments are based on the differences in 
c     concentrations of the regular model species before and after
c     the regular model chemistry.  This is essentially an adjustment
c     for production or decay. In GOAT (Geographic Ozone Apportionment
c     Technology), allocation of ozone production is based solely on 
c     location.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       ldrk    L  flag to deterine if current hour is in darkness
c       igrid   I  grid number
c       icell   I  the X grid location of current cell
c       jcell   I  the X grid location of current cell
c       kcell   I  the vertical grid location of current layer
c       delo3   R  change in O3 concentrations 
c       delno   R  change in NO concentrations 
c       delno2  R  change in NO2 concentrations 
c       delvoc  R  change in VOC concentrations (carbon weighted sum 
c                  of VOC species)
c       dtime   R  change in time for current time step 
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     05/22/97   --gwilson--  Now re-calculates the sum of VOC after
c                             adjustment using vocwt.
c     07/19/02   --gwilson--  Added seperate source area map for each grids.
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
      logical   ldrk
      integer   igrid
      integer   icell
      integer   jcell
      integer   kcell
      real      delo3
      real      delno
      real      delno2
      real      delvoc
      real      dtime
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer*4 idxcel, idx, imap, iregn, i
      real*4    sumnox, sumvoc, sumo3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the source region for this cell ---
c
      imap = igrmap(igrid,icell,jcell)
      if( imap .LT. 0 .OR. imap .GT. nregin ) goto 9999
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  ipsa3d(igrid)-1+ icell + ncol(igrid)*(jcell-1) + 
     &                               ncol(igrid)*nrow(igrid)*(kcell-1)
c
c   --- loop over the NOx tracer species, calculate total NOx tracer ---
c
      sumnox = 0.
      do 10 i=iptnox,iptvoc-1
         idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         sumnox = sumnox + ptconc(idx)
   10 continue
c
c   ---- apply the decay to each tracer species, apportioning by 
c        contribution to total ----
c
      if( sumnox .GT. 0. ) then
         do 20 i=iptnox,iptvoc-1
             idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
             ptconc(idx) = ptconc(idx) + 
     &                    (delno + delno2) * ptconc(idx) / sumnox
             ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   20    continue
      endif
c
c  --- adjust the sum of NOx species to account for change ---
c
      sumnox = sumnox + (delno + delno2)
c
c   --- loop over the VOC tracer species, calculate total VOC tracer ---
c
      sumvoc = 0.
      do 30 i=iptvoc,ipto3n-1
         idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         sumvoc = sumvoc + ptconc(idx) * vocwt(i)
   30 continue
c
c   ---- apply the decay to each tracer species, apportioning by 
c        contribution to total ----
c
      if( sumvoc .GT. 0. ) then
         do 40 i=iptvoc,ipto3n-1
             idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
             ptconc(idx) = ptconc(idx) + 
     &               delvoc * ptconc(idx) * vocwt(i) / sumvoc
             ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   40    continue
      endif
c
c  --- recalculate the sum of VOC tracers ----
c
      sumvoc = 0.
      do i=iptvoc,ipto3n-1
         idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         sumvoc = sumvoc + ptconc(idx) * vocwt(i)
      enddo
c
c   --- do the daytime case ----
c
      if( .NOT. ldrk ) then
c
c   --- ozone production, add to local ozone tracers,
c
          if( delo3 .GT. 0 ) then
c
c   --- loop over all the O3N and O3V species ---
c
            do 50 i=ipto3n,ipttim-1
              if( ptname(i)(7:8) .EQ. 'IC' .OR.
     &                     ptname(i)(7:8) .EQ. 'BC' ) goto 50
c
c   ---- the tracer names are of the the form O3N001rrr and
c        O3VN001rrr, where rrr is the region number,  This
c        will match two times, one for O3N and 1 for O3V ---
c
              read(ptname(i)(7:9),'(I3)') iregn
              if( iregn .EQ. imap ) then
                  idx = idxcel + ncol(igrid)*nrow(igrid)*
     &                                               nlay(igrid)*(i-1)
                  ptconc(idx) = ptconc(idx) + delo3/2.0
                  ptconc(idx) = MAX(BNDLPT,ptconc(idx))
              endif
   50       continue
          endif
      endif
c
c   --- if ozone in decay or hour is a darkness hour
c       allocate to all Ozone tracers ---
c
      if(ldrk .OR. delo3 .LT. 0. ) then
         sumo3 = 0.
         do 70 i=ipto3n,ipttim-1
            idx = idxcel + ncol(igrid)*nrow(igrid)* nlay(igrid)*(i-1)
            sumo3 = sumo3 + ptconc(idx)
   70    continue
         if( sumo3 .GT. 0. ) then
            do 80 i=ipto3n,ipttim-1
               idx = idxcel + ncol(igrid)*nrow(igrid)* nlay(igrid)*(i-1)
               ptconc(idx) = ptconc(idx) + 
     &                                delo3 * ptconc(idx) / sumo3
               ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   80       continue
         endif
      endif
c
c   --- decay the deacying tracer species ---
c
      if(ntrtim .GT. 0 ) then
         do 41 i=ipttim+1,nsaspc,2
            idx = idxcel + ncol(igrid)*nrow(igrid)* nlay(igrid)*(i-1)
            ptconc(idx) = ptconc(idx) * 
     &                                     EXP( -0.08333333 * dtime )
            ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   41    continue
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
