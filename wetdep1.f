      subroutine wetdep(igrid,ncol,nrow,nlay,nspcs,deltat,deltax,deltay,
     &                  depth,tempk,press,cwc,pwc,densfac,idfin,
     &                  conc,fluxes,depfld,dtout,ipa_cel)
c
c-----CAMx v4.02 030709
c 
c     WETDEP modifies vertical concentration profiles for a given grid via 
c     precipitation processes.  This subroutine has been completely rewritten
c     for CAMx v4.
c 
c     Copyright 2003
c     ENVIRON International Corporation
c           
c     Modifications:
c      06/15/03   Fixed bug in scaling of tmass
c
c     Input arguments:
c        igrid               grid index 
c        ncol                number of columns 
c        nrow                number of rows 
c        nlay                number of layers 
c        nspcs               number of total species
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
c        dtout               output frequency (minutes)
c        ipa_cel             gridded array to identify if cell is
c                            in a IPRM sub-domain
c             
c     Output arguments: 
c        conc                concentration field (umol/m3, ug/m3)
c        fluxes              boundary mass fluxes (umol, ug)
c        depfld              2-D array of wet deposited mass (mol/ha, g/ha)
c                            and surface liquid concentrations (mol/l, g/l)
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

      include 'deposit.com'
      include 'section.inc'
c
c======================== Source Apportion Begin =======================
c======================== DDM Begin ====================================
c
      include 'tracer.com'
      real fc2r, fr2c, delnul, c0trac(MXTRSP), tmtrac(MXTRSP)
      real delnox, delvoc, delo3, c0nox, c0voc, c0o3
      real concnox, concvoc, conco3, dnlnox, dnlvoc, dnlo3
      logical lzerc, levap, lddmopt2
      data lddmopt2 /true/
c
c======================== DDM End ======================================
c========================= Source Apportion End ========================
c
c
c======================== Process Analysis Begin ====================================
c
      include 'procan.com'
c
      integer ipa_cel(ncol,nrow,nlay)
c
c========================= Process Analysis End =====================================
c
      real deltax(nrow),tempk(ncol,nrow,nlay),
     &     press(ncol,nrow,nlay),cwc(ncol,nrow,nlay),
     &     pwc(ncol,nrow,nlay),depth(ncol,nrow,nlay)
      integer idfin(ncol,nrow)
      real conc(ncol,nrow,nlay,nspcs),depfld(ncol,nrow,3*nspcs)
      real*8 fluxes(nspcs,11)
      real c0(MXSPEC),rr(MXLAYA),volrat(MXLAYA),tmass(MXSPEC)
      real delr(MXSPEC)
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
        jcl = j
        i1 = 2 
        i2 = ncol-1 
        if (igrid.eq.1) then 
          if (ibeg(j).eq.-999) goto 10 
          i1 = ibeg(j) 
          i2 = iend(j) 
        endif 
        do 20 i = i1,i2
          icl = i
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
            rr(k) = (volrat(k)/1.0e-7)**1.14         ! rainfall rate (mm/hr)
          enddo
c
c======================== Source Apportion Begin =======================
c
          c0o3 = 0.
          c0nox = 0.
          c0voc = 0.
          dnlo3 = 0.
          dnlnox = 0.
          dnlvoc = 0.
c
c======================== Source Apportion End =======================
c
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
c======================== Source Apportion Begin =======================
c======================== DDM Begin ====================================
c
            if( lddm .OR. ltrace .AND. tectyp .NE. RTRAC ) then
               if (ltop) then
                  do l=1,MXTRSP
                     c0trac(l) = 0.
                     tmtrac(l) = 0.
                  enddo
               else
                  do l=1,MXTRSP
                     c0trac(l) = tmtrac(l) / rainvol
                  enddo
               endif
               delo3 = 0.
               delnox = 0.
               delvoc = 0.
               conco3 = 0.
               concnox = 0.
               concvoc = 0.
            endif
c
c======================== DDM End ======================================
c======================== Source Apportion End =======================
c
c
c-----Calculate scavenging for soluble gas species
c
            do 40 l = 1,ngas
              if( henry0(l).LT.1.e-6 ) goto 40
              delc = 0.
              delm = 0.
              if (ltop) tmass(l) = 0.
              c0(l) = tmass(l) / rainvol
              convfac = densfac*(273./tempk(i,j,k))*(press(i,j,k)/1013.)
              cmin = bdnl(l)*convfac
              conc(i,j,k,l) = amax1(cmin,conc(i,j,k,l))
c
              call scavrat(.false.,lcloud,lfreez,rr(k),
     &                     tempk(i,j,k),cwc(i,j,k),depth(i,j,k),rhoair,
     &                     conc(i,j,k,l),c0(l),henry0(l),tfact(l),
     &                     diffrat(l),0.,0.,hlaw,cgas,dscav,gscav,
     &                     gdiff,ascav)
c
c-----Calculate the delc limit if Henry's Law equilibrium is reached
c
              dclim = (cgas*hlaw - c0(l))*volrat(k)
              dclim2 = (2.0*cgas*hlaw - c0(l))*volrat(k)
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
c
c-----Range checks on scavenging if delc > 0 (i.e., cell to rain)
c       delc < conc(i,j,k,l)-cmin keeps cell concentration > cmin
c       delc < dclim because scavenging does not super-saturate
c     Range checks on scavenging if delc < 0 (i.e., rain to cell)
c       decreasing rain volume may super-saturate rain
c       delc < dclim2 prevents super-saturation exceeding factor of 2
c       delc > -c0(l)*volrat(k) keeps rain concentration > 0.
c       delc > MIN(0.,dclim) prevents rain going below saturation
c
              if (delc.gt.0.) then
                delc = MIN(delc,conc(i,j,k,l)-cmin,dclim)
              else
                delc = MIN(delc,dclim2)
                delc = MAX(delc,-c0(l)*volrat(k),MIN(0.,dclim))
              endif
c
c======================== Source Apportion Begin =======================
c
              if( ltrace .AND. tectyp .NE. RTRAC ) then
c
c   --- call scavrate again with nothing in the rain to get the
c       concentration moving into the rain ---
c
                 call scavrat(.false.,lcloud,lfreez,rr(k),
     &                     tempk(i,j,k),cwc(i,j,k),depth(i,j,k),rhoair,
     &                     conc(i,j,k,l),0.,henry0(l),tfact(l),
     &                     diffrat(l),0.,0.,hlaw,cgas,dscav,gscav,
     &                     gdiff,ascav)
                 dnlim = (cgas*hlaw - cmin)*volrat(k)
                 scav = (dscav*(1. - gdiff) + gscav*gdiff)*deltat
                 scav = amin1(scav, 16.)
                 scav = amax1(scav,-16.)
                 delnul = conc(i,j,k,l)*(1. - exp(-scav))
                 delnul = MIN(delnul,conc(i,j,k,l)-cmin,dnlim)
                 delnul = MAX(delnul,0.)
c
c  --- calculate the amount moving into the rain ---
c
                 if( lo3sp(l) ) then
                     delo3 = delo3  + delc
                     conco3 = conco3 + conc(i,j,k,l)
                     dnlo3 = dnlo3  + delnul
                     c0o3 = c0o3  + c0(l)
                 else if( lnoxsp(l) ) then
                     delnox = delnox  + delc
                     concnox = concnox + conc(i,j,k,l)
                     dnlnox = dnlnox  + delnul
                     c0nox = c0nox  + c0(l)
                 else if( lvocsp(l) ) then
                     delvoc = delvoc  + delc * crbnum(l)
                     concvoc = concvoc + conc(i,j,k,l) * crbnum(l)
                     dnlvoc = dnlvoc  + delnul * crbnum(l)
                     c0voc = c0voc  + c0(l) * crbnum(l)
                 endif
              endif
c
c======================== Source Apportion End =======================
c
c
c======================== DDM Begin =======================
c
              if( lddm ) then
c
c   --- call scavrate again with nothing in the rain to get the
c       flux into the rain ---
c
c                 if( c0(l) .EQ. 0.0 ) then
c                    delnul = delc
c                 else
c                    call scavrat(.false.,lcloud,lfreez,rr(k),
c     &                     tempk(i,j,k),cwc(i,j,k),depth(i,j,k),rhoair,
c     &                     conc(i,j,k,l),0.,henry0(l),tfact(l),
c     &                     diffrat(l),0.,0.,hlaw,cgas,dscav,gscav,
c     &                     gdiff,ascav)
c                    dnlim = (cgas*hlaw - cmin)*volrat(k)
c                    scav = (dscav*(1. - gdiff) + gscav*gdiff)*deltat
c                    scav = amin1(scav, 16.)
c                    scav = amax1(scav,-16.)
c                    delnul = conc(i,j,k,l)*(1. - exp(-scav))
c                    if (delnul.gt.0.) then
c                      delnul = MIN(delnul,conc(i,j,k,l)-cmin,dnlim)
c                    else
c                      delnul = 0.
c                    endif
c                 endif
c
c  --- calculate the relative fluxes between cell and rain and
c      determine whether the cell concentration is at lower bound ---
c
c                 fc2r = MAX( 0., delnul/conc(i,j,k,l) )
c                 if( c0(l) .LT. 1.0e-20 ) then
c                    fr2c = 0.
c                 else
c                   fr2c = MAX( 0. , (delnul-delc)/c0(l) )
c                 endif
c
c  --- Convert delc to flux ---
c
                 if( delc.GT.0. ) then
                    fc2r = delc/conc(i,j,k,l)
                    fr2c = 0.
                 elseif( c0(l) .GT. 1.0e-20 ) then
                    fc2r = 0.
                    fr2c = -delc/c0(l)
                 else
                    fc2r = 0.
                    fr2c = 0.
                 endif
c
c  --- check for saturation ---
c
                 if ( ((tmass(l)+delc*cellvol)/rainvol)/
     &                 ((conc(i,j,k,l)-delc)*hlaw) .GT. 0.999 ) then
                    fc2r = 0.
                 endif
c
c  --- check for special situations ---
c
                 levap = .FALSE.
                 if( kbot.GT.1 .and. k.EQ.kbot ) levap = .TRUE.
                 lzerc = .FALSE.
                 if( conc(i,j,k,l) .LE. 1.1*cmin ) lzerc = .TRUE.
c
c  --- call routine to apply the changes to sensitivities ---
c 
                 call adjddmc0(l,igrid,icl,jcl,kcl,fc2r,fr2c,
     &                         c0trac,tmtrac,cellvol,rainvol,lzerc,
     &                         levap)
              endif
c
c======================== DDM End =======================
c
c
c-----Update the cell concentration and rain mass
c
              conc(i,j,k,l) = conc(i,j,k,l) - delc
              delm = delc*cellvol
              tmass(l) = tmass(l) + delm
c
c======================== Process Analysis Begin ====================================
c
              if( lipr .AND. ipa_cel(i,j,k) .GT. 0 ) then
                ipa_idx = ipa_cel(i,j,k)
c
c-----Change from wet deposition
c
                cipr(IPR_WDEP, ipa_idx, l) =
     &                              cipr(IPR_WDEP, ipa_idx, l) - delc
              endif
c
c========================= Process Analysis End =====================================
c
 40         continue
c
c======================== Source Apportion Begin =======================
c
            if( ltrace .AND. tectyp .NE. RTRAC ) then
               levap = .FALSE.
               if( kbot.GT.1 .and. k.EQ.kbot ) levap = .TRUE.
               call adjstc0(igrid,icl,jcl,kcl,delo3,delnox,delvoc,
     &                 conco3,concnox,concvoc,dnlo3,dnlnox,dnlvoc,
     &                 c0o3,c0nox,c0voc,c0trac,tmtrac,cellvol,rainvol,
     &                 levap)
            endif
c
c======================== Source Apportion End =======================
c
c
c-----Calculate scavenging for particulate species
c
            if (naero .gt. 0) then
c->   calculate wet diameter - bkoo (11/05/03)
              if (kph2o.ne.nspec+1) then ! mechanism 4
                isempty = 1
                qt = 0.0
                roprta = 0.0
                do l = ngas+1,nspec
                  if (dcut(l,2).lt.2.51 .and. l.ne.kph2o) then ! dry PM2.5
                    if (conc(i,j,k,l).gt.bdnl(l)) isempty = 0
                    qt = qt + conc(i,j,k,l)
                    roprta = roprta + conc(i,j,k,l) / roprt(l)
                  endif
                enddo
                rfin = 0.
                if (isempty.eq.0) then ! avoid empty particles
                  roprta = ( qt + conc(i,j,k,kph2o) ) / ( roprta + 
     &                       conc(i,j,k,kph2o) / roprt(kph2o) )
                  do isec = 1, nsecfin
                    psize = ( 1. + conc(i,j,k,kph2o) / qt )**0.33333
     &                     * diadep(isec) * 1.e-6
                    call scavrat(.true.,lcloud,lfreez,rr(k),
     &                       tempk(i,j,k),cwc(i,j,k),depth(i,j,k),
     &                       rhoair,0.,0.,0.,0.,0.,psize,roprta,hlaw,
     &                       cgas,dscav,gscav,gdiff,ascav)
                    rfin = rfin + (1.-exp(-ascav*deltat)) * wfdep(isec)
                  enddo
                endif
                do l = ngas+1,nspec
                  if (dcut(l,2).lt.2.51) then
                    delr(l) = rfin
                  else
                    rcrs = 0.
                    do isec = nsecfin + 1, nsecdep
                      psize = diadep(isec) * 1.e-6
                      call scavrat(.true.,lcloud,lfreez,rr(k),
     &                       tempk(i,j,k),cwc(i,j,k),depth(i,j,k),
     &                       rhoair,0.,0.,0.,0.,0.,psize,roprt(l),hlaw,
     &                       cgas,dscav,gscav,gdiff,ascav)
                      rcrs = rcrs + (1.-exp(-ascav*deltat))*wfdep(isec)
                    enddo
                    delr(l) = rcrs
                  endif
                enddo
              elseif (kph2o_1.ne.nspec+1) then ! mechanism 6
                kwtr = (kph2o_1 - ngas) / nsec + 1
                if (nsec.eq.1) kwtr = kph2o_1 - ngas
                do isec = 1, nsec
                  isempty = 1
                  qt = 0.0
                  roprta = 0.0
                  do iaero = 1, naero
                    if (iaero.ne.kwtr) then
                      l = ngas + (iaero-1)*nsec + isec
                      if (conc(i,j,k,l).gt.bdnl(l)) isempty = 0
                      qt = qt + conc(i,j,k,l)
                      roprta = roprta + conc(i,j,k,l) / roprt(l)
                    endif
                  enddo
                  ascav = 0.
                  if (isempty.eq.0) then ! avoid empty particles
                    psize = sqrt(dcut(ngas+isec,1)*dcut(ngas+isec,2))
                    psize = ( 1. + conc(i,j,k,kph2o_1-1+isec) / qt
     &                       )**0.33333 * psize * 1.e-6
                    roprta = ( qt + conc(i,j,k,kph2o_1-1+isec) ) /
     &                       ( roprta + conc(i,j,k,kph2o_1-1+isec) /
     &                                       roprt(kph2o_1) )
                    call scavrat(.true.,lcloud,lfreez,rr(k),
     &                       tempk(i,j,k),cwc(i,j,k),depth(i,j,k),
     &                       rhoair,0.,0.,0.,0.,0.,psize,roprta,hlaw,
     &                       cgas,dscav,gscav,gdiff,ascav)
                  endif
                  do iaero = 1, naero
                    l = ngas + (iaero-1)*nsec + isec
                    delr(l) = 1. - exp(-ascav*deltat)
                  enddo
                enddo
              else
                do l = ngas+1,nspec
                  psize = sqrt(dcut(l,1)*dcut(l,2))*1.e-6
                  call scavrat(.true.,lcloud,lfreez,rr(k),
     &                       tempk(i,j,k),cwc(i,j,k),depth(i,j,k),
     &                       rhoair,0.,0.,0.,0.,0.,psize,roprt(l),hlaw,
     &                       cgas,dscav,gscav,gdiff,ascav)
                  delr(l) = 1. - exp(-ascav*deltat)
                enddo
              endif
c<-
              do 50 l = ngas+1,nspec
                delc = 0.
                delm = 0.
                if (ltop) tmass(l) = 0.
                cmin = bdnl(l)
                conc(i,j,k,l) = amax1(cmin,conc(i,j,k,l))
cbk                psize = sqrt(dcut(l,1)*dcut(l,2))*1.e-6

cbk                call scavrat(.true.,lcloud,lfreez,rr(k),
cbk     &                       tempk(i,j,k),cwc(i,j,k),depth(i,j,k),
cbk     &                       rhoair,0.,0.,0.,0.,0.,psize,roprt(l),hlaw,
cbk     &                       cgas,dscav,gscav,gdiff,ascav)

cbk                delc = conc(i,j,k,l)*(1. - exp(-ascav*deltat))
c->
                delc = conc(i,j,k,l) * delr(l)
c<-
                delc = amin1(delc,conc(i,j,k,l)-cmin)
c
                conc(i,j,k,l) = conc(i,j,k,l) - delc
                delm = delc*cellvol
                tmass(l) = tmass(l) + delm
c
c======================== Process Analysis Begin ====================================
c
                if( lipr .AND. ipa_cel(i,j,k) .GT. 0 ) then
                  ipa_idx = ipa_cel(i,j,k)
c
c-----Change from wet deposition
c
                  cipr(IPR_WDEP, ipa_idx, l) =
     &                              cipr(IPR_WDEP, ipa_idx, l) - delc
                endif
c
c========================= Process Analysis End =====================================
c
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
            do l = 1,nspec
              conc(i,j,kbot,l) = conc(i,j,kbot,l) + tmass(l)/cellvol
              tmass(l) = 0.
            enddo
c
c-----Otherwise rain reaches the ground, increment deposition flux arrays
c
          else
            do l = 1,nspec
              fluxes(l,11) = fluxes(l,11) - tmass(l)
              do ll = 1,navspc
                if (l .eq. lavmap(ll)) then
                  depfld(i,j,navspc+ll) = depfld(i,j,navspc+ll) +
     &                          1.e-2*tmass(l)/(deltax(j)*deltay)
                  depfld(i,j,2*navspc+ll) = depfld(i,j,2*navspc+ll) +
     &                    1.e-9*(tmass(l)/rainvol)*deltat/(60.*dtout)
                  goto 100
                endif
              enddo
 100          continue
            enddo
          endif
c
 20     continue
 10   continue
c
      return
      end
