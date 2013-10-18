      subroutine average(losat,igrd,dt,ncol,nrow,nlay,nlayav,
     &               nspav,nspc,lmap,tempk,press,conc,avcnc,ipa_cel)
c
c-----CAMx v4.02 030709
c
c     AVERAGE computes time-averaged concentrations
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications: 
c        01/30/02    --gwilson--  Added code for RTRAC probing tool
c
c     Input arguments:
c        losat              .TRUE. if concentrations are tracer species
c        igrd               grid index
c        dt                 time step for present grid concentration (s)
c        ncol               number of columns
c        nrow               number of rows
c        nlay               number of layers in instantaneous array
c        nlayav             number of layers in average array
c        nspav              number of average species
c        nspc               number of species in conc array
c        lmap               mapping array for average species
c        tempk              temperature field (K)
c        press              pressure field (mb)
c        conc               instant species concentration (umol/m3)
c        avcnc              average species concentration (gas=ppm,
c                                                          other=ug/m3)
c        ipa_cel            gridded array to identify if cell is
c                            in a IPRM sub-domain
c
c     Output arguments:
c        avcnc              average species concentration (gas=ppm,
c                                                          other=ug/m3)
c
c     Routines Called:
c        none
c
c     Called by:
c        CAMx
c        FGAVRG
c
      include "camx.prm"
      include "camx.com"
      include "bndary.com"
      include "chmstry.com"
c
c========================= Source Apportion Begin ==============================
c
      include "tracer.com"
      include "rtracchm.com"
c
c========================= Source Apportion End ==============================
c
c
c========================= Process Analysis Begin ==============================
c
      include "procan.com"
c
      integer ipa_cel(ncol,nrow,nlay)
c
c========================= Process Analysis End ==============================
c
      logical lgas, losat
      dimension tempk(ncol,nrow,nlay),press(ncol,nrow,nlay),
     &          avcnc(ncol,nrow,nlayav,nspav),conc(ncol,nrow,nlay,nspc),
     &          lmap(nspc)
c
c-----Entry point
c
c-----Increment running average
c
      dtfact = dt/(dtout*60.)
      do 40 l = 1,nspav
        lsp = lmap(l) 
        lgas = .true. 
        if( lsp .GT. ngas .AND. .NOT. losat) then 
          convfac = 1.
          lgas = .false. 
        endif 
        if( losat .AND. tectyp .EQ. RTRAC .AND. lsp .GT. nrtgas ) then
          convfac = 1.
          lgas = .false. 
        endif 
        do 30 j = 2,nrow-1
          i1 = 2
          i2 = ncol-1
          if (igrd.eq.1) then
            if (ibeg(j).eq.-999) goto 30
            i1 = ibeg(j)
            i2 = iend(j)
          endif
          do i = i1,i2
c
            do k=1,nlayav
                if (lgas) then
                  tmp = 273./tempk(i,j,k)*press(i,j,k)/1013.
                  convfac = 1./(densfac*tmp)
                endif
                avcnc(i,j,k,l) = convfac*conc(i,j,k,lsp)*dtfact + 
     &                           avcnc(i,j,k,l)
c
c========================= Process Analysis Begin ==============================
c
                if( .NOT. losat .AND. lipr .AND. 
     &                                ipa_cel(i,j,k) .GT. 0 ) then
                   ipa_idx = ipa_cel(i,j,k) 
c
c-----Save the units conversion factor for use in IPR post-processing
c
                   cipr(IPR_CONV, ipa_idx, lsp) = 
     &                    cipr(IPR_CONV, ipa_idx, lsp)+ convfac * dtfact
                 endif
c
c========================= Process Analysis End ==============================
c
            enddo
          enddo
  30    continue
  40  continue
c
      return
      end
