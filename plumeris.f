      subroutine plumeris(nlay,hght,temp,dtdz,wind,hstk,dstk,tstk,
     &                    vstk,prise)
c
c-----CAMx v4.02 030709
c
c     PLUMERIS calculates plume rise above a given elevated point source.
c     Methodology based on the EPA TUPOS Gaussian Plume Model
c     (Turner and Catalano, 1986).
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications: 
c        none
c
c     Input arguments:
c        nlay                number of layers
c        hght                layer interface heights (m)
c        temp                temperature (K)
c        dtdz                potential temperature lapse rate (K/m)
c        wind                total horizontal wind (m/s)
c        hstk                stack height (m)
c        dstk                stack diameter (m)
c        tstk                stack temperature (K)
c        vstk                stack velocity (m/s)
c
c     Output arguments:
c        prise               plume rise (m)
c
c     Routines called:
c        none
c
c     Called by:
c        EMISS
c        PIGINIT
c
      dimension hght(nlay),temp(nlay),dtdz(nlay),wind(nlay)
      logical lfirst
c
      data grav /9.8/
c
c-----Entry point
c
      lfirst = .true.
      rise = 0.
c
c-----Find beginning layer; determine vertical coordinates relative
c     to stack-top 
c
      do k = 1,nlay
        if (hstk.lt.hght(k)) then
          kstk = k
          goto 5
        endif
      enddo
 5    ztop = hght(kstk) - hstk
      zbot = 0.
      stkt = amax1(tstk,temp(kstk)+1.)
c
c-----Determine downwash factor as a function of stack Froude number
c
      fr = temp(kstk)*vstk*vstk/(grav*dstk*(stkt - temp(kstk)))
      dwfact = 1.
      if (fr.ge.3.) then
        if (wind(kstk).ge.vstk) then
          goto 999
        elseif (wind(kstk).gt.vstk/1.5 .and. wind(kstk).lt.vstk) then
          dwfact = 3.*(vstk - wind(kstk))/vstk
        endif
      endif
c
c-----Neutral-unstable momentum rise and stack buoyancy flux
c
      umrise = 3.*dstk*vstk/wind(kstk)
      if (umrise.gt.ztop) then
        wsum = wind(kstk)*ztop
        do k = kstk+1,nlay 
          wsum = wsum + wind(k)*(hght(k) - hght(k-1))
          wavg = wsum/(hght(k) - hstk)
          umrise = 3.*dstk*vstk/wavg
          if (umrise.lt.hght(k)-hstk) goto 10
        enddo
      endif
  10  bflux0 = grav*vstk*dstk*dstk*(stkt - temp(kstk))/(4.*stkt)
      bflux = bflux0
c
c-----Top of layer loop; determine stability
c
  20  continue
      kstab = kstk
      if (lfirst) then
        if( (kstk.gt.1 .and. hstk.lt.((hght(kstk)+hght(kstk-1))/2.))
     &      .or. (kstk.eq.nlay) ) kstab = kstk-1
      else
        kstab = kstk-1
        bflux = rflux
        ztop = hght(kstk) - hstk
        zbot = hght(kstk-1) - hstk
      endif
      if (dtdz(kstab).gt.1.5e-3) goto 30 
c
c-----Neutral-unstable buoyancy rise
c
      ubrise1 = 30.*(bflux/wind(kstk))**0.6 + zbot
      ubrise2 = 24.*(bflux/wind(kstk)**3)**0.6 *
     &          (hstk + 200.*(bflux/wind(kstk)**3))**0.4 + zbot
      ubrise = amin1(ubrise1,ubrise2)
      iuse = 1
      if (ubrise.eq.ubrise2) iuse = 2
c
c-----Find maximum of neutral-unstable momentum and buoyancy rise
c
      if (lfirst .and. umrise.gt.ubrise) then
        rise = umrise
        goto 999
      else
        rise = ubrise
        if (rise.le.ztop) goto 999
      endif
c
c-----Neutral-unstable residual buoyancy flux 
c
      if (iuse.eq.1) then
        rflux = wind(kstk)*((rise - ztop)/30.)**(5./3.)
      else
        rflux = 0.0055*(rise - ztop)*wind(kstk)**3 /
     &          (1. + hstk/(rise - ztop))**(2./3.)
      endif
      kstk = kstk + 1
      if (kstk.gt.nlay .or. rflux.le.0.) goto 999
      lfirst = .false.
      goto 20
c
c-----Stable momentum rise
c
  30  continue
c
      smrise = 0.646*(vstk*vstk*dstk*dstk/(stkt*wind(kstk)))**(1./3.) *
     &         sqrt(temp(kstk))/(dtdz(kstab)**(1./6.))
      smrise = amin1(smrise,umrise)
c
c-----Stable buoyancy rise
c
      sbrise1 = (1.8*bflux*temp(kstk)/(wind(kstk)*dtdz(kstab)) + 
     &           zbot*zbot*zbot)**(1./3.)
      sbrise2 = (4.1*bflux*temp(kstk)/(bflux0**(1./3.)*dtdz(kstab)) +
     &           zbot**(8./3.))**(3./8.)
      sbrise = amin1(sbrise1,sbrise2)
      iuse = 1
      if (sbrise.eq.sbrise2) iuse = 2
c
c-----Find maximum between stable momentum and (2/3)*bouyancy rise
c
      if (lfirst .and. smrise.gt.(2.*sbrise/3.)) then
        rise = smrise
        goto 999
      elseif (sbrise.le.ztop) then
        rise = 2.*sbrise/3.
        goto 999
      else
        rise = sbrise
      endif
c
c-----Stable residual buoyancy flux
c
      if (iuse.eq.1) then
        rflux = bflux - (0.56*dtdz(kstab)*wind(kstk)/temp(kstk))*
     &          (ztop*ztop*ztop - zbot*zbot*zbot) 
      else
        rflux = bflux - 0.24*dtdz(kstab)*bflux0**(1./3.)/temp(kstk)*
     &          (ztop**(8./3.) - zbot**(8./3.))
      endif
      kstk = kstk + 1
      if (kstk.gt.nlay .or. rflux.le.0.) then
        rise = 2.*rise/3.
        goto 999
      endif
      lfirst = .false.
      goto 20
c
c-----Add stack height to downwash-corrected plume rise and return
c
 999  prise = hstk + dwfact*rise 
c
      return
      end
