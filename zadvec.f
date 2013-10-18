      subroutine zadvec(losat,igrd,ncol,nrow,nlay,nspc,nsen,deltat,dx,
     &                  dy,densfac,idfin,tpconc,depth,entrn,
     &                  dilut,tempk,press,species,conc,fluxes,
     &                  tpsens,sens,tarray2,ipa_xy,ipa_lay)
c
c-----CAMx v4.02 030709
c
c     ZADVEC drives vertical transport of concentrations.  The operation
c     is performed on the regular model species (losat = False) or source
c     apportionment concentrations (losat = True).
c     This version also performs vertical transport of sensitivities 
c     if DDM is enabled.
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        12/3/99   Top boundary condition is converted from ppm to umol/m3
c                  using extrapolated density in an overlying layer
c                  that is assumed to be the same thickness as the top
c                  layer of the model
c       12/07/01   added instructions for OMP
c       01/30/02   Added code for RTRAC probing tool
c
c     Input arguments:
c        losat             flag to determine if this call is processing 
c                          OSAT concentrations.
c        igrd              grid index
c        ncol              number of columns
c        nrow              number of rows
c        nlay              number of layers
c        nspc              number of species
c        nsen              number of species times number of DDM parameters
c                          or 1, whichever is larger
c        deltat            time step (s)
c        dx                cell size in x-direction (m)
c        dy                cell size in y-direction (m)
c        densfac           density conversion factor (mol/m3)
c        idfin             map of nested grids in this grid
c        tpconc            top concentrations (ppm)
c        depth             layer depth (m)
c        entrn             entrainment rate (m/s)
c        dilut             dilution rate (m/s)
c        tempk             temperature (K)
c        press             pressure (mb)
c        species           species names
c        conc              species concentrations (umol/m3)
c        tpsens            top sensitivities (ppm/parameter unit)
c        sens              sensitivity coefficients(umol/m3/parameter unit)
c        tarray2           CPU timing arguments (s)
c        ipa_xy            2-D gridded array to identify if cell is
c                          in a IPRM sub-domain
c        ipa_lay           3-D gridded array to identify which IPRM sub-domain
c                          each layer is in
c
c     Output arguments:
c        conc              species concentrations (umol/m3)
c        fluxes            fluxes across the boundaries (umol)
c        sens              sensitivity coefficients(umol/m3/parameter unit)
c
c     Routines Called:
c        VRTSLV
c
c     Called by:
c        EMISTRNS
c
      include "camx.prm"
      include "bndary.com"
      include "chmstry.com"
      include "filunit.com"
      include "tracer.com"
      include "rtracchm.com"
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
c
      integer ipa_xy(ncol,nrow)
      integer ipa_lay(ncol,nrow,nlay)
      real fc1(MXLAYA+1), fc2(MXLAYA+1), fc3(MXLAYA+1)
      logical ldoipts
c
c========================= Process Analysis End =====================================
c
      dimension conc(ncol,nrow,nlay,nspc)
      dimension tpconc(nspc),tarray2(2)
      dimension sens(ncol,nrow,nlay,nsen),tpsens(nsen)
      character*10 species(nspc)
      real entrn(ncol,nrow,nlay),dilut(ncol,nrow,nlay),
     &     depth(ncol,nrow,nlay),tempk(ncol,nrow,nlay),
     &     press(ncol,nrow,nlay),dx(nrow)
      integer idfin(ncol,nrow)
      real c1d(MXLAYA+1),d1d(MXLAYA),ent1d(MXLAYA),dil1d(MXLAYA)
      real sen1d((MXLAYA+1)*MXTRSP)
      real*8 fluxes(nspc,11),fluxtop
      logical losat
c
c-----Entry point
c
c-----Vertical advection, entrainment, and dilution
c
      write(*,'(a20,$)') 'z advection ......' 
      write(iout,'(a20,$)') 'z advection ......' 
c
c
c$omp parallel default(shared)
c$omp&  private(ispc,i,j,k,i1,i2,ent1d,dil1d,d1d,c1d,rhop,
c$omp&          rho,rhom,convfac,ls,isen,ioff,sen1d,fluxtop,
c$omp&          drodz,ldoipts,ipa_idx,fc1,fc2,fc3)
c
c$omp do schedule(dynamic)
c
      do 61 ispc = 1,nspc
        do 60 j = 2,nrow-1
          i1 = 2
          i2 = ncol - 1
          if (igrd.eq.1) then
            if (ibeg(j).eq.-999) goto 60
            i1 = ibeg(j)
            i2 = iend(j)
          endif
          do 50 i = i1,i2
c
c-----Skip cells occupied by child grids
c
            if (idfin(i,j).gt.igrd) goto 50
            do k = 1,nlay
              ent1d(k) = entrn(i,j,k)
              dil1d(k) = dilut(i,j,k)
              d1d(k) = depth(i,j,k)
              c1d(k) = conc(i,j,k,ispc)
            enddo
            rho = press(i,j,nlay)/tempk(i,j,nlay)
            rhom = press(i,j,nlay-1)/tempk(i,j,nlay-1)
            drodz = 2.*(rho - rhom)/(d1d(nlay) + d1d(nlay-1))
            rhop = rho + drodz*d1d(nlay)
            convfac = densfac*273.*rhop/1013.
            if (ispc.gt.ngas .and. .not.losat) then
              convfac = 1.
            endif
            if (losat .AND. tectyp .EQ. RTRAC .AND. ispc.GT.nrtgas) then
              convfac = 1.
            endif
            c1d(nlay+1) = convfac*tpconc(ispc)
c
c================================DDM Begin==============================
c
            if (lddm) then
              ls = 0
              do isen = 1,nddmsp
                ioff = iptddm(ispc)+isen-1
                do k = 1,nlay
                  ls = ls + 1
                  sen1d(ls) = sens(i,j,k,ioff)
                enddo
                ls = ls + 1
                sen1d(ls) = convfac*tpsens(ioff)
              enddo
            endif
c
c================================DDM End================================
c
c
c======================== Process Analysis Begin ====================================
c
            ldoipts = .FALSE.
            if ( .NOT. ltrace .AND. lipr .AND. ipa_xy(i,j) .GT. 0 )
     &                                                 ldoipts = .TRUE.

c
c========================= Process Analysis End =====================================
c
c-----Solve vertical mass adjustments using Crank-Nicholson solver
c
            call vrtslv(nlay,i,j,igrd,deltat,ent1d,dil1d,d1d,c1d,
     &                  fluxtop,sen1d,species(ispc),fc1,fc2,fc3,ldoipts)
c
c======================== Process Analysis Begin ====================================
c
            if ( ldoipts ) then
               do k=1,nlay
                  if( ipa_lay(i,j,k) .GT. 0 ) then
                    ipa_idx = ipa_lay(i,j,k)
c
c-----Concentration change from vertical advection
c
                    cipr(IPR_BADV, ipa_idx, ispc) =
     &                            cipr(IPR_BADV, ipa_idx, ispc) + fc1(k)
c
                    cipr(IPR_TADV, ipa_idx, ispc) =
     &                           cipr(IPR_TADV, ipa_idx, ispc) + fc2(k)
c
                    cipr(IPR_DADV, ipa_idx, ispc) =
     &                           cipr(IPR_DADV, ipa_idx, ispc) + fc3(k)
c
                  endif
               enddo
            endif
c
c========================= Process Analysis End =====================================
c
            do k = 1,nlay
              conc(i,j,k,ispc) = c1d(k)
            enddo
c
            if (fluxtop.lt.0) then
              fluxes(ispc,9)  = fluxes(ispc,9) - fluxtop*dx(j)*dy*deltat
            else
             fluxes(ispc,10) = fluxes(ispc,10) - fluxtop*dx(j)*dy*deltat
            endif
c
c================================DDM Begin=============================
c
            if (lddm) then
              ls = -1
              do isen = 1, nddmsp
                ls = ls + 1
                ioff = iptddm(ispc)+isen-1
                do k = 1, nlay
                  ls = ls + 1
                  sens(i,j,k,ioff) = sen1d(ls)
                enddo
              enddo
            endif
c
c================================DDM End===============================
c
  50      continue
  60    continue
  61  continue
c
c$omp enddo
c
c$omp end parallel
c
      tcpu = dtime(tarray2) 
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1) 
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
c
      call flush(6)
      call flush(iout)
      return
      end
