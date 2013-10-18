c*** APCASA
c
      subroutine apcasa(ldrk,igrid,icell,jcell,kcell,
     &                  delo3,delno,delno2,delvoc,delh22,delhn3,dtime)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine makes the "chemistry" adjustments to the tracer 
c     species.  The adjustments are based on the differences in 
c     concentrations of the regular model species before and after
c     the regular model chemistry.  This routine does the attribution
c     using the APCA.
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
c       delh22  R  change in H2O2 concentrations 
c       delhn3  R  change in HNO3 concentrations 
c       dtime   R  change in time for current time step 
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     05/20/97   --gwilson--  Original development
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
      real      delh22
      real      delhn3
      real      dtime
c
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   idxcel, idx, jdx, ivoc, inox, i
      real      sumnox, sumvoc, sumo3
      real      bionox, biovoc, facvoc, facnox, facbio, o3used
      real       useall
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- calculate the index of the cell in the grid ---
c
      idxcel =  ipsa3d(igrid)-1 + icell + ncol(igrid)*(jcell-1) + 
     &                               ncol(igrid)*nrow(igrid)*(kcell-1)
      o3used = 0.
      useall = 0.
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
      bionox = 0.
      if( sumnox .GT. 0. ) then
         do 20 i=iptnox,iptvoc-1
             idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
             ptconc(idx) = ptconc(idx) + 
     &                         (delno + delno2) * ptconc(idx) / sumnox
             ptconc(idx) = MAX(BNDLPT,ptconc(idx))
             if(ptname(i)(4:6) .EQ. '001') bionox = bionox + ptconc(idx)
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
      biovoc = 0.
      if( sumvoc .GT. 0. ) then
         do 40 i=iptvoc,ipto3n-1
             idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
             ptconc(idx) = ptconc(idx) + 
     &                        delvoc * ptconc(idx) * vocwt(i) / sumvoc
             ptconc(idx) = MAX(BNDLPT,ptconc(idx))
             if( ptname(i)(4:6) .EQ. '001' ) biovoc = biovoc + 
     &                                         ptconc(idx) * vocwt(i)
   40    continue
      endif
c
c  --- adjust the sum of VOC species to account for change ---
c
      sumvoc = 0.
      do i=iptvoc,ipto3n-1
         idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         sumvoc = sumvoc + ptconc(idx) * vocwt(i)
      enddo
c
c  --- set the biogenics contribution factor ---
c
      facnox = 0.
      facvoc = 0.
      if( sumnox .NE. 0. ) facnox = bionox/sumnox
      if( sumvoc .NE. 0. ) facvoc = biovoc/sumvoc
      facbio = MIN( facnox, facvoc )
c
c   --- do the daytime case ----
c
      if( .NOT. ldrk ) then
c
c   --- if ozone in production check for NOx or VOC limited ---
c
          if( delo3 .GT. 0 ) then
c
c   --- VOC limited, increase the O3V tracers first ----
c
              if( delhn3 .NE. 0. .AND. delh22/delhn3 .LE. 0.35 ) then
                 ivoc = iptvoc - 1
                 do 50 i=ipto3v,ipttim-1
                    idx = idxcel + 
     &                       ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
                    ivoc = ivoc + 1
                    jdx = idxcel + 
     &                      ncol(igrid)*nrow(igrid)*nlay(igrid)*(ivoc-1)
c
c   --- add O3 to biogenics group (group 1) based on it's
c       mimimum potential contribution ---
c
                    if( ptname(i)(4:6) .EQ. '001' ) then
                       if( biovoc .GT. 0. ) then
                           ptconc(idx) = ptconc(idx) + 
     &                      (ptconc(jdx) * vocwt(ivoc) / biovoc) * 
     &                                                 facbio * delo3
                           o3used = o3used + 
     &                      (ptconc(jdx) * vocwt(ivoc) / biovoc) * 
     &                                                  facbio * delo3
                           useall = useall + 
     &                      (ptconc(jdx) * vocwt(ivoc) / biovoc) * 
     &                                                    facbio * delo3
                       endif
c
c  --- add O3 to other groups based on contribution of VOC conc ---
c
                    else
                        if( sumvoc .GT. 0. ) then
                           ptconc(idx) = ptconc(idx) + (ptconc(jdx) * 
     &                                    vocwt(ivoc) / sumvoc) * delo3
                           o3used = o3used + (ptconc(jdx) * vocwt(ivoc) 
     &                                                / sumvoc) * delo3
                           useall = useall + (ptconc(jdx) * vocwt(ivoc) 
     &                                                / sumvoc) * delo3
                        endif
                     endif
                     ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   50            continue
c
c  --- add any ozone not accounted for to O3N tracers, based on
c      contribution to anthropognic NOx ---
c
                 if( sumnox .NE. bionox .AND. o3used .LT. delo3 ) then 
                    inox = iptnox - 1
                    do 60 i=ipto3n,ipto3v-1
                       idx = idxcel + 
     &                         ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
                       inox = inox + 1
                       jdx = idxcel + 
     &                      ncol(igrid)*nrow(igrid)*nlay(igrid)*(inox-1)
                       if( ptname(i)(4:6) .NE. '001' ) then
                           ptconc(idx) = ptconc(idx) + 
     &                      (ptconc(jdx) / (sumnox - bionox)) * 
     &                                                    (delo3-o3used)
                           useall = useall + 
     &                      (ptconc(jdx) / (sumnox - bionox)) * 
     &                                                    (delo3-o3used)
                       endif
   60               continue
                 endif
c
c   --- NOx limited, increase the O3N tracers first ----
c
              else
                 inox = iptnox - 1
                 do 70 i=ipto3n,ipto3v-1
                    idx = idxcel + 
     &                        ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
                    inox = inox + 1
                    jdx = idxcel + 
     &                      ncol(igrid)*nrow(igrid)*nlay(igrid)*(inox-1)
c
c   --- add O3 to biogenics group (group 1) based on it's
c       mimimum potential contribution ---
c
                    if( ptname(i)(4:6) .EQ. '001' ) then
                       if( bionox .GT. 0. ) then
                           ptconc(idx) = ptconc(idx) + 
     &                         (ptconc(jdx) / bionox) * facbio * delo3
                           o3used = o3used + 
     &                          (ptconc(jdx) / bionox) * facbio * delo3
                           useall = useall + 
     &                          (ptconc(jdx) / bionox) * facbio * delo3
                       endif
c
c  --- add O3 to other groups based on contribution of NOx conc ---
c
                    else
                        if( sumnox .GT. 0. ) then
                           ptconc(idx) = ptconc(idx) + 
     &                              (ptconc(jdx) / sumnox) * delo3
                           o3used = o3used + 
     &                              (ptconc(jdx) / sumnox) * delo3
                           useall = useall + 
     &                              (ptconc(jdx) / sumnox) * delo3
                        endif
                     endif
                     ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   70            continue
c
c  --- add any ozone not accounted for to O3V tracers, based on
c      contribution to anthropognic VOC ---
c
                 if( sumvoc .NE. biovoc .AND. o3used .LT. delo3 ) then 
                    ivoc = iptvoc - 1
                    do 80 i=ipto3v,ipttim-1
                       idx = idxcel + 
     &                     ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
                       ivoc = ivoc + 1
                       jdx = idxcel + 
     &                     ncol(igrid)*nrow(igrid)*nlay(igrid)*(ivoc-1)
                       if( ptname(i)(4:6) .NE. '001' ) then
                           ptconc(idx) = ptconc(idx) + 
     &                        (ptconc(jdx) * vocwt(ivoc) / 
     &                             (sumvoc - biovoc) ) * (delo3-o3used)
                           useall = useall +
     &                        (ptconc(jdx) * vocwt(ivoc) / 
     &                             (sumvoc - biovoc) ) * (delo3-o3used)
                       endif
   80               continue
                 endif
              endif
          endif
      endif
c
c   --- if ozone in decay or hour is a darkness hour
c       allocate to all Ozone tracers ---
c
      if(ldrk .OR. delo3 .LT. 0. ) then
         sumo3 = 0.
         do 90 i=ipto3n,ipttim-1
            idx = idxcel + ncol(igrid)*nrow(igrid)* nlay(igrid)*(i-1)
            sumo3 = sumo3 + ptconc(idx)
   90    continue
         if( sumo3 .GT. 0. ) then
            do 11 i=ipto3n,ipttim-1
               idx = idxcel + ncol(igrid)*nrow(igrid)* nlay(igrid)*(i-1)
               ptconc(idx) = ptconc(idx) + delo3 * ptconc(idx) / sumo3
               ptconc(idx) = MAX(BNDLPT,ptconc(idx))
   11       continue
         endif
      endif
c
c   --- decay the deacying tracer species ---
c
      if(ntrtim .GT. 0 ) then
         do 21 i=ipttim+1,nsaspc,2
            idx = idxcel + ncol(igrid)*nrow(igrid)* nlay(igrid)*(i-1)
            ptconc(idx) = ptconc(idx) * EXP( -0.08333333 * dtime )
            ptconc(idx) = MAX(BNDLPT,ptconc(idx) )
   21    continue
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
