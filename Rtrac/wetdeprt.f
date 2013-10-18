      subroutine wetdeprt(igrid,ncol,nrow,nlay,nrtsp,deltat,deltax,
     &                  deltay,depth,tempk,press,cwc,pwc,densfac,idfin,
     &                  conc)
c
c-----CAMx v4.02 030709
c 
c     This version of is for the RTRAC species.
c     WETDEP modifies vertical concentration profiles for a given grid via 
c     precipitation processes.  This subroutine has been completely rewritten
c     for CAMx v4.
c 
c     Copyright 2003
c     ENVIRON International Corporation
c           
c     Modifications:
c        5/01/03   Taken from WETDEP for regular model
c
c     Input arguments:
c        igrid               grid index 
c        ncol                number of columns 
c        nrow                number of rows 
c        nlay                number of layers 
c        nrtsp               number of species in RTRAC array
c        deltat              time step (s)
c        deltax              cell size in x-direction (m)
c        deltay              cell size in y-direction (m)
c        depth               cell depth (m)
c        tempk               temperature field (K)
c        press               pressure field (mb)
c        cwc                 cloud water content (g/m3)
c        pwc                 precipitation water content (g/m3)
c        densfac             factor to convert to umol/m3
c        idfin               map of nested grids in this grid
c        conc                concentration field (umol/m3, ug/m3)
c             
c     Output arguments: 
c        conc                concentration field (umol/m3, ug/m3)
c             
c     Routines called: 
c        none
c             
c     Called by: 
c        EMISTRNS
c 
      include 'camx.prm'
      include 'bndary.com'
      include 'chmstry.com'
      include 'tracer.com'
      include 'rtracchm.com'
c
      real deltax(nrow),tempk(ncol,nrow,nlay),
     &     press(ncol,nrow,nlay),cwc(ncol,nrow,nlay),
     &     pwc(ncol,nrow,nlay),depth(ncol,nrow,nlay)
      integer idfin(ncol,nrow)
      real conc(ncol,nrow,nlay,nrtsp)
      real c0(MXTRSP),rr(MXLAYA),volrat(MXLAYA),tmass(MXTRSP)
      logical lcloud,ltop,lfreez
c
      data rd /287./         ! Dry air gas constant (J/K/kg)
      data rhoh2o /1.e6/     ! water density (g/m3)
c
c-----Entry point
c
c-----Loop over rows and columns
c
      do 10 j = 2,nrow-1 
        i1 = 2 
        i2 = ncol-1 
        if (igrid.eq.1) then 
          if (ibeg(j).eq.-999) goto 10 
          i1 = ibeg(j) 
          i2 = iend(j) 
        endif 
        do 20 i = i1,i2
          if (idfin(i,j).gt.igrid) goto 20
c
c-----Scan column for layers containing precipitation bottom/top
c
          kbot = 0
          ktop = 0
          do k = 1,nlay
            if (pwc(i,j,k).ge.cwmin) then
              kbot = k
              goto 25
            endif
          enddo
          goto 20

  25      continue
          if (kbot.eq.nlay) goto 20
          ncnt = 1
          do k = kbot+1,nlay
            if (pwc(i,j,k).lt.cwmin) then
              ktop = k-1
              goto 26
            endif
            ncnt = ncnt + 1
          enddo
          ktop = nlay
  26      continue
          if (ncnt.eq.1) goto 20
c
c-----Determine rainfall rate profile
c
          do k = 1,nlay
            volrat(k) = 0.
            rr(k) = 0.
          enddo
          do k = kbot,ktop
            volrat(k) = pwc(i,j,k)/rhoh2o            ! drop volume/air volume
            rr(k) = (volrat(k)/1.08e-7)**1.225       ! rainfall rate (mm/hr)
          enddo
c
c-----Loop over layers and species for raining columns
c
          ltop = .true.
          do 30 k = ktop,kbot,-1
            kcl = k
            lcloud = .false.
            lfreez = .false.
            if (cwc(i,j,k).ge.cwmin) lcloud = .true.
            if (tempk(i,j,k).lt.tamin) lfreez = .true.
            cellvol = deltax(j)*deltay*depth(i,j,k)
            rainvol = volrat(k)*cellvol
            rhoair = 100.*press(i,j,k)/(rd*tempk(i,j,k))
c
c-----Calculate scavenging for gas species
c
            do 40 l = 1,nrtgas
              delc = 0.
              delm = 0.
              if (ltop) then
                c0(l) = 0.
                tmass(l) = 0.
              endif
              convfac = densfac*(273./tempk(i,j,k))*(press(i,j,k)/1013.)
              cmin = rtlbnd(l)*convfac
              conc(i,j,k,l) = amax1(cmin,conc(i,j,k,l))

              call scavrat(.false.,lcloud,lfreez,rr(k),
     &                     tempk(i,j,k),cwc(i,j,k),depth(i,j,k),rhoair,
     &                     conc(i,j,k,l),c0(l),rthlaw(l),rttfact(l),
     &                     rtdrate(l),0.,0.,hlaw,cgas,dscav,gscav,
     &                     gdiff,ascav)
c
c-----Calculate the delc limit if Henry's Law equilibrium is reached
c
              dclim = (cgas*hlaw - c0(l))*volrat(k)
c
c-----Combine scavenging of gas disolved in cloudwater (dscav) and
c     dissolution of gas into rainwater (gscav), limited by diffusion
c     factor.  Rapid diffusion defeats dscav and enables gscav
c
              scav = (dscav*(1. - gdiff) + gscav*gdiff)*deltat
              scav = amin1(scav, 16.)
              scav = amax1(scav,-16.)
c
              delc = conc(i,j,k,l)*(1. - exp(-scav))
              if (scav.gt.0.) then
                delc = amin1(delc,conc(i,j,k,l)-cmin)
                if (dclim.gt.0.) delc = amin1(delc,dclim)
              elseif (scav.lt.0.) then
                delc = amax1(delc,dclim)
              endif
c
              conc(i,j,k,l) = conc(i,j,k,l) - delc
              delm = delc*cellvol
              tmass(l) = tmass(l) + delm
              c0(l) = c0(l) + delm/rainvol
              if (k.gt.1 .and. volrat(k-1).gt.0.)
     &          c0(l) = c0(l)*volrat(k)/volrat(k-1)
 40         continue
c
c-----Calculate scavenging for particulate species
c
            if (naero .gt. 0) then
              do 50 l = nrtgas+1,nrtrac
                delc = 0.
                delm = 0.
                if (k.eq.ktop) tmass(l) = 0.
                cmin = rtlbnd(l)
                conc(i,j,k,l) = amax1(cmin,conc(i,j,k,l))
                psize = sqrt(rtlcut(l)*rtucut(l))*1.e-6

                call scavrat(.TRUE.,lcloud,lfreez,rr(k),
     &                      tempk(i,j,k),cwc(i,j,k),depth(i,j,k),
     &                      rhoair,0.,0.,0.,0.,0.,psize,rtdens(l),hlaw,
     &                      cgas,dscav,gscav,gdiff,ascav)

                delc = conc(i,j,k,l)*(1. - exp(-ascav*deltat))
                delc = amin1(delc,conc(i,j,k,l)-cmin)
c
                conc(i,j,k,l) = conc(i,j,k,l) - delc
                delm = delc*cellvol
                tmass(l) = tmass(l) + delm
 50           continue
            endif
            ltop = .false.
 30       continue
c
c-----If rain evaporates before reaching the ground, return all mass back
c     to layer KBOT
c
          if (kbot.gt.1) then
            cellvol = deltax(j)*deltay*depth(i,j,kbot)
            do l = 1,nrtrac
              conc(i,j,kbot,l) = conc(i,j,kbot,l) + tmass(l)/cellvol
            enddo
          endif
c
 20     continue
 10   continue
c
      return
      end
