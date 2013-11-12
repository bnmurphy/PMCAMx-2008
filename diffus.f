      subroutine diffus(igrd,ncol,nrow,nlay,nspc,nsen,deltat,dx,dy,
     &                  idfin,vdep,rkx,rky,rkv,depth,tempk,press,mapscl,
     &                  conc,fluxes,depfld,sens,tarray2,strz,strxy,
     &                  ipa_xy,ipa_lay)
c
c-----CAMx v4.02 030709
c
c     DIFFUS drives 3-D diffusion of concentrations.  This version also
c     performs diffusion of sensitivities if DDM is enabled.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003 
c     ENVIRON International Corporation
c          
c     Modifications:
c        4/17/00   Revised diffusion equations to weight fluxes by density
c       12/07/01   added instructions for OMP
c        1/13/03   added deposited mass array
c
c     Input arguments:
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
c        idfin             map of nested grids in this grid
c        vdep              deposition velocity (m/s)
c        rkx/y             horizontal diffusion coefficient (m2/s)
c        rkv               vertical diffusion coefficient (m2/s)
c        depth             layer depth (m)
c        tempk             temperature (K)
c        press             pressure (mb)
c        mapscl            map-scale factor at cell centroid
c        conc              species concentrations (umol/m3)
c        sens              DDM sensitivities (umol/m3/parameter unit)
c        tarray2           CPU timing arguments (s)
c        strz              string for labeling the z diffusion process
c        strxy             string for labeling the x/y diffusion process
c        ipa_xy            2-D gridded array to identify if cell is
c                          in a IPRM sub-domain
c        ipa_lay           3-D gridded array to identify which IPRM sub-domain
c                          each layer is in
c
c     Output arguments:
c        conc              species concentrations (umol/m3)
c        sens              DDM sensitivities (umol/m3/parameter unit)
c        fluxes            fluxes across the boundaries (umol)
c        depfld            2-D array of dry deposited mass (mol/ha, g/ha)
c
c     Routines Called:
c        VDIFFIMP
c	 SUBDOMAIN 	!BNM 9-23-09
c
c     Called by:
c        EMISTRNS
c
      include "camx.prm"
      include "filunit.com"
      include "bndary.com"
      include "chmstry.com"
c
c=============================DDM Begin================================
c
      include "tracer.com"
c
c=============================DDM End==================================
c
c
c======================== Process Analysis Begin ====================================
c
      include "procan.com"
c
      integer ipa_xy(ncol,nrow)
      integer ipa_lay(ncol,nrow,nlay)
      real fcup(MXLAYA),fcdn(MXLAYA)
      logical ldoipts
c
c========================= Process Analysis End =====================================
c
      dimension conc(ncol,nrow,nlay,nspc),vdep(ncol,nrow,nspc)
      dimension sens(ncol,nrow,nlay,nsen)
      dimension tarray2(2)
      real depfld(ncol,nrow,3*nspc)
      real rkx(ncol,nrow,nlay),rky(ncol,nrow,nlay),
     &     rkv(ncol,nrow,nlay),depth(ncol,nrow,nlay),
     &     tempk(ncol,nrow,nlay),press(ncol,nrow,nlay),
     &     mapscl(ncol,nrow),dx(nrow)
      integer idfin(ncol,nrow)
      dimension c1d(MXLAYA+MXLAYA*MXTRSP),d1d(MXLAYA),rk1d(MXLAYA),
     &          ro1d(MXLAYA)
      real cnc(MXCOLA,MXROWA),sns(MXCOLA,MXROWA,MXTRSP),
     &     rho(MXCOLA,MXROWA)
      real*8 fluxes(nspc*nlay,13),fluxbot
      character*20 strz, strxy

      integer subd(ncol,nrow)	!BNM 9-23-09
c
c-----Entry point
c
      call subdomain(subd)

c-----Vertical diffusion
c
      write(*,'(a20,$)') strz
      write(iout,'(a20,$)') strz
c
c$omp parallel default(shared)
c$omp&  private(ispc,i1,i2,i,j,k,ro1d,d1d,c1d,rk1d,lk,
c$omp&          ioff,fluxbot,isen,rokp1,ldoipts,fcup,
c$omp&          fcdn,ipa_idx)
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
c-----Skip cells occupied by child grids; load 1-D arrays
c
            if (idfin(i,j).gt.igrd) goto 50
            do k = 1,nlay
              ro1d(k) = press(i,j,k)/tempk(i,j,k)
              d1d(k) = depth(i,j,k)
              c1d(k) = conc(i,j,k,ispc)
            enddo
            do k = 1,nlay-1
              rokp1 = press(i,j,k+1)/tempk(i,j,k+1)
              rk1d(k) = rkv(i,j,k)*(ro1d(k) + rokp1)/2.
            enddo
            rk1d(nlay) = 0.
c
c
c=============================DDM Begin================================
c
c------Load sensitivities
c
            if (lddm) then
              lk = nlay
              do isen = 1,nddmsp
                ioff = iptddm(ispc)+isen-1
                do k = 1,nlay
                  lk = lk + 1
                  c1d(lk) = sens(i,j,k,ioff)
                enddo
              enddo
            endif
c
c=============================DDM End==================================
c
c
c======================== Process Analysis Begin ====================================
c
            ldoipts = .FALSE.
            if( .NOT. ltrace .AND. lipr .AND. ipa_xy(i,j) .GT. 0) 
     &                                                 ldoipts = .TRUE.
c
c========================= Process Analysis End =====================================
c
            call vdiffimp(nlay,deltat,vdep(i,j,ispc),d1d,ro1d,rk1d,c1d,
     &                    nddmsp,fcup,fcdn,ldoipts)
c
c======================== Process Analysis Begin ====================================
c
            if( ldoipts ) then
                do k = 1,nlay
                  if( ipa_lay(i,j,k) .GT. 0) then
                     ipa_idx  = ipa_lay(i,j,k)
                     if( k .GT. 1 ) then
                        cipr(IPR_BDIF, ipa_idx, ispc) =
     &                     cipr(IPR_BDIF, ipa_idx, ispc) + fcdn(k)
                     else
                        cipr(IPR_DDEP, ipa_idx, ispc) =
     &                     cipr(IPR_DDEP, ipa_idx, ispc) + fcdn(k)
                     endif
                     cipr(IPR_TDIF, ipa_idx, ispc) =
     &                     cipr(IPR_TDIF, ipa_idx, ispc) + fcup(k)
                  endif
                enddo
            endif
c
c========================= Process Analysis End =====================================
c
            do k = 1,nlay
              conc(i,j,k,ispc) = c1d(k)
            enddo
            fluxbot = -vdep(i,j,ispc)*conc(i,j,1,ispc)
            fluxes(1+(ispc-1)*nlay,11) = fluxes(1+(ispc-1)*nlay,11) 
     &			+ fluxbot*dx(j)*dy*deltat * subd(i,j)  	!i=col, j=row
		!Keep Dry Deposition at Flux 11
            do ll = 1,navspc
              if (ispc .eq. lavmap(ll)) then
                depfld(i,j,ll) = depfld(i,j,ll) - 1.e-2*fluxbot*deltat
                goto 100
              endif
            enddo
 100        continue
c
c=============================DDM Begin================================
c
            if (lddm) then
              lk = nlay
              do isen = 1,nddmsp
                ioff = iptddm(ispc)+isen-1
                do k = 1,nlay
                  lk = lk + 1
                  sens(i,j,k,ioff) = c1d(lk)
                enddo
              enddo
            endif
c
c=============================DDM End==================================
c
  50      continue
  60    continue
  61  continue
c
c$omp end parallel
c$omp master
c
      tcpu = dtime(tarray2) 
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1) 
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
c
c-----Perform explicit horizontal diffusion
c
      write(*,'(a20,$)') strxy
      write(iout,'(a20,$)') strxy
c
c$omp end master
c
c$omp parallel default(shared)
c$omp&  private(ispc,i1,i2,i,j,k,ioff,rho,cnc,sns,rhoavg,scl,
c$omp&          dfxp,dfxm,dfyp,dfym,fxp,fxm,fyp,fym,isen)
c
c$omp do schedule(dynamic)
c
      do 91 ispc = 1,nspc
        do 90 k = 1,nlay
          do 80 j = 1,nrow
            do 70 i = 1,ncol
              rho(i,j) = press(i,j,k)/tempk(i,j,k)
              cnc(i,j) = conc(i,j,k,ispc)/rho(i,j)
  70        continue
  80      continue
c
c=============================DDM Begin================================
c
          if (lddm) then
            do isen = 1,nddmsp
              ioff = iptddm(ispc)+isen-1
              do j=1,nrow
                do i = 1,ncol
                  sns(i,j,isen) = sens(i,j,k,ioff)/rho(i,j)
                 enddo
               enddo
            enddo
          endif
c
c=============================DDM End==================================
c
          do 85 j = 2,nrow-1
            i1 = 2
            i2 = ncol - 1
            if (igrd.eq.1) then
              if (ibeg(j).eq.-999) goto 85
              i1 = ibeg(j)
              i2 = iend(j)
            endif
            do 75 i = i1,i2
              if (idfin(i,j).gt.igrd) goto 75
c
              rhoavg = (rho(i,j) + rho(i+1,j))/2.
              scl = (mapscl(i,j) + mapscl(i+1,j))/2.
              dfxp = scl*rhoavg*rkx(i,j,k)*deltat/dx(j)/dx(j)
              rhoavg = (rho(i,j) + rho(i-1,j))/2.
              scl = (mapscl(i,j) + mapscl(i-1,j))/2.
              dfxm = scl*rhoavg*rkx(i-1,j,k)*deltat/dx(j)/dx(j)
c
              rhoavg = (rho(i,j) + rho(i,j+1))/2.
              scl = (mapscl(i,j) + mapscl(i,j+1))/2.
              dfyp = scl*rhoavg*rky(i,j,k)*deltat/dy/dy
              rhoavg = (rho(i,j) + rho(i,j-1))/2.
              scl = (mapscl(i,j) + mapscl(i,j-1))/2.
              dfym = scl*rhoavg*rky(i,j-1,k)*deltat/dy/dy
c
              fxp  = (cnc(i+1,j) - cnc(i,j))*dfxp
              fxm  = (cnc(i,j) - cnc(i-1,j))*dfxm
              fyp  = (cnc(i,j+1) - cnc(i,j))*dfyp
              fym  = (cnc(i,j) - cnc(i,j-1))*dfym
              conc(i,j,k,ispc) = conc(i,j,k,ispc) + mapscl(i,j)*
     &                           ((fxp - fxm) + (fyp - fym))
c
c======================== Process Analysis Begin ====================================
c
              if( .NOT. ltrace .AND. lipr .AND. 
     &                                      ipa_lay(i,j,k) .GT. 0 ) then
                     ipa_idx = ipa_lay(i,j,k)
                     cipr(IPR_WDIF, ipa_idx, ispc) =
     &                       cipr( IPR_WDIF, ipa_idx, ispc) -
     &                                       mapscl(i,j)*fxm

                     cipr(IPR_EDIF, ipa_idx, ispc) =
     &                       cipr(IPR_EDIF, ipa_idx, ispc) +
     &                                       mapscl(i,j)*fxp

                     cipr(IPR_SDIF, ipa_idx, ispc) =
     &                       cipr(IPR_SDIF, ipa_idx, ispc) -
     &                                       mapscl(i,j)*fym

                     cipr(IPR_NDIF, ipa_idx, ispc) =
     &                       cipr(IPR_NDIF, ipa_idx, ispc) +
     &                                       mapscl(i,j)*fyp
              endif
c
c========================= Process Analysis End =====================================
c
c
c=============================DDM Begin================================
c
              if (lddm) then
                do isen =1,nddmsp
                  fxp  = (sns(i+1,j,isen) - sns(i,j,isen))*dfxp
                  fxm  = (sns(i,j,isen) - sns(i-1,j,isen))*dfxm
                  fyp  = (sns(i,j+1,isen) - sns(i,j,isen))*dfyp
                  fym  = (sns(i,j,isen) - sns(i,j-1,isen))*dfym
                  ioff = iptddm(ispc)+isen-1
                  sens(i,j,k,ioff) = sens(i,j,k,ioff) + mapscl(i,j)*
     &                              ((fxp - fxm) + (fyp - fym))
                enddo
              endif
c
c=============================DDM End==================================
c
  75        continue
  85      continue
  90    continue
  91  continue
c
c$omp end parallel
c$omp master
c
      tcpu = dtime(tarray2) 
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1) 
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
c
c$omp end master
c
      call flush(6)
      call flush(iout)
      return
      end
