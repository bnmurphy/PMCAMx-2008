c*** CHEMRT
c
      subroutine chemrt(igrid, icell, jcell, kcell, pcell, tcell,
     &                  cold, cnew, oh, o3, no3, dtime, convfac,
     &                  irt_cel)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine performs chemistry adjustments on the reactive
c     tracer gasses using information from the host chemistry solution.  
c
c     Copyright 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       igrid   I  number of the grid containing the cell
c       icell   I  X index of the cell
c       jcell   I  Y index of the cell
c       kcell   I  Z index of the cell
c       pcell   R  cell pressure (mb)
c       tcell   R  cell temperature (K)
c       cold    R  species concentrations before host chemistry (ppm)
c       cnew    R  species concentrations after host chemistry (ppm)
c       oh      R  OH radical concentration (ppm)
c       o3      R  ozone concentration (ppm)
c       no3     R  no3 radical concentration (ppm)
c       dtime   R  time step (hours)
c       convfac R  conversion factor: umol/m3 = ppm * convfac
c       irt_cel I  index for RTRAC receptor cells
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'chmstry.com'
      include 'grid.com'
      include 'tracer.com'
      include 'rtracchm.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer igrid
      integer icell
      integer jcell
      integer kcell
      real    pcell
      real    tcell
      real    cold(MXTRSP)
      real    cnew(MXTRSP)
      real    oh
      real    o3
      real    no3
      real    dtime
      real    convfac
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   i
      integer   idx
      integer   idxcel
      real      drate(MXTRSP)
      real      factor
      real      prod
      real      rkoh
      real      rko3
      real      rkno3
      real      rkspc
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c --- calculate the index of the cell in the grid ---
c
      idxcel =  ipsa3d(igrid)-1+ icell + ncol(igrid)*(jcell-1) + 
     &                               ncol(igrid)*nrow(igrid)*(kcell-1)
c
c --- factor is a temp/press adjustment for second order
c     rate constants that were converted from cm3 molec-1 sec-1 to
c     ppm-1 min-1 at STP, also convert to per hour ---
c
      factor = (298.0/tcell)*(pcell/1013.) * 60.0
c
c --- calculate the total decay rate and any secondary production
c     rate for each tracer.  get host photolysis rates from rk in the
c     chmstry.com include file ---
c
      do i=1,nrtgas
         if( jnum(i) .GT. 0 ) then
            rkspc = rk(jnum(i))
         else
            rkspc = 0.0
         endif 
         rkoh  = arrhen(aoh(i),eaoh(i),boh(i),troh(i),tcell) * factor
         rko3  = arrhen(ao3(i),eao3(i),bo3(i),tro3(i),tcell) * factor
         rkno3 = arrhen(ano3(i),eano3(i),bno3(i),trno3(i),tcell)
     &           * factor
         drate(i) = oh*rkoh + o3*rko3 + no3*rkno3 + rkspc*rtjfact(i)
         prod = 0.0
         if ( lsecnd(i) ) then
            prod = cnew(ksec(i)) - cold(ksec(i))*EXP(-drate(i)*dtime)
            prod = prod * convfac
            prod = MAX(prod,0.0)
         endif
         idx = idxcel + ncol(igrid)*nrow(igrid)*nlay(igrid)*(i-1)
         ptconc(idx) = ptconc(idx)*EXP(-drate(i)*dtime) + prod
         ptconc(idx) = MAX(rtlbnd(i),ptconc(idx))
      enddo
c
c --- if this is a receptor cell, accumulate average decay rates
c
      if (irt_cel .GT. 0 ) then 
         do i=1,nrtgas
            rcpdcy(irt_cel, i) = rcpdcy(irt_cel,i) + 
     &                              drate(i)*dtime/(dtout/60.)
         enddo
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
c
      function arrhen(a,ea,b,tref,temp)
c
c --- Calculate a rate constant using k = A*(T/Tref)^B*exp(Ea/T)
c
      real arrhen
c
      arrhen = a*((temp/tref)**b)*exp(ea/temp)
c
      end
