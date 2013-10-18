      subroutine scavrat(laero,lcloud,lfreez,pmm,temp,cwc,dz,
     &                   rhoair,conc,c0,rk0,tfact,difrat,prtdia,rhoprt,
     &                   hlaw,cgas,dscav,gscav,gdiff,ascav)
c
c-----CAMx v4_hg 030825
c 
c     SCAVRAT calculates wet (liquid) scavenging rates for gases and aerosols.
c     Rates are determined for both in-cloud and below-cloud removal.  Gas 
c     removal rates are calculated for dissolved gas in cloud water and for 
c     uptake of ambient gas.  All in-cloud aerosols are assumed to be in cloud
c     water (liquid and frozen); below cloud scavenging of aerosols is 
c     dependent on particle size.
c 
c     Copyright 2003
c     ENVIRON International Corporation
c           
c     Modifications:
c        02/11/04   Updated empirical rainfall relationships
c
c     Input arguments:
c        laero               aerosol flag
c        lcloud              in-cloud flag (F=no cloud water)
c        lfreez              frozen precip flag
c        pmm                 precip rate (mm/hr)
c        temp                temperature (K)
c        cwc                 cloud water content (g/m3)
c        dz                  cell depth (m)
c        rhoair              atmospheric density (kg/m3)
c        conc                cell gas or aerosol concentration (umol/m3)
c        c0                  initial gas concentration in rain (umol/m3)
c        rk0                 Henry's Law constant (M/atm)
c        tfact               H-Law temperature dependence (K)
c        difrat              Ratio of H2O to gas diffusivity
c        prtdia              mean aerosol size (m)
c        rhoprt              aerosol density (g/m3)
c             
c     Output arguments: 
c        hlaw                Scale Henry's Law constant for temperature (M/M)
c        cgas                Equilibrium gas-phase concentration (umol/m3)
c        dscav               Gas scavenging rate by droplet removal (1/s)
c        gscav               Gas scavenging rate by dissolution in rain (1/s)
c        gdiff               Diffusion factor for gases (dimensionless, 0 to 1)
c        ascav               Aerosol scavenging rate (1/s)
c             
c     Routines called: 
c        none
c             
c     Called by: 
c        WETDEP
c 
      real cwc
      real nuair,muair,muh2o,kc
      logical laero,lcloud,lfreez
c
c-----Constants
c
      data pi /3.1415927/
      data rconst /8.206e-2/ !l.atm/mol.K
      data rhoh2o /1.e6/     !g/m3
      data difh2o /2.3e-5/   !m2/s
      data muair /1.8e-5/    !kg/ms
      data muh2o /1.e-3/     !kg/ms
      data boltz /1.38e-23/  !J/K
      data xmfp /6.5e-8/     !m
      data cldeff /0.9/
c
c-----Entry point
c
      dscav = 0.
      gscav = 0.
      gdiff = 1.
      ascav = 0.
c
c-----Calculate environmental parameters
c
      drpdia = 9.0e-4*(pmm**0.21)      !rain drop diameter (m)
      drpvel = 3100.*drpdia            !rain drop fall speed (m/s)
      nuair = muair/rhoair             !air molecular diffusivity (m2/s)
      cscav = 4.2e-7*pmm*cldeff/drpdia !cloud droplet scavenging rate (1/s)

      if (laero) goto 1000 
c
c-----Gas scavenging
c
      if (lfreez) return
      cgas = conc
      caq  = 0.
      hlaw = rk0*rconst*temp*exp(tfact*(1./298. - 1./temp))
c
c-----If in cloud, partition total gas into equilibrium aqueous and gas phase,
c     and calculate scavenging rate for dissolved gasses in cloud droplets
c
      if (lcloud) then
        cgas = conc/(1. + hlaw*cwc/rhoh2o)
        caq  = conc - cgas
        dscav = cscav*caq/conc
      endif
c
c-----Calculate scavenging rate for ambient gas dissolving into rain
c
      diff  = difh2o/difrat
      term1 = (drpvel*drpdia/nuair)**0.5
      term2 = (nuair/diff)**0.333
      kc    = diff/drpdia*(2. + 0.6*term1*term2)
      expo  = -6.*kc*dz/(drpdia*drpvel*hlaw)
      gdiff = 1. - exp(expo)
      gscav = 2.8e-7*pmm/(dz*conc)*(cgas*hlaw - c0)
c
      return
c
c-----Aerosol scavenging
c
 1000 continue
c
c-----Scavenging rate for aerosol in cloud droplets is set equal to
c     cloud water scavenging rate
c
      if (lcloud) then
        ascav = cscav
      else
c
c-----Calculate scavenging rate of dry aerosols below cloud as f(size)
c
        if (lfreez) return
        reynold = drpdia*drpvel/(2.*nuair)
        power = amin1(7.6,0.55*prtdia/xmfp)
        scf = 1. + (2.514 + 0.8*exp(-power))*xmfp/prtdia
        difbrwn = boltz*temp*scf/(3.*pi*muair*prtdia)
        schmidt = nuair/difbrwn
        stoke = drpvel*prtdia*prtdia*rhoprt*scf/(9000.*muair*drpdia)

        top = 1.2 + alog(1. + reynold)/12.
        bot = 1.0 + alog(1. + reynold)
        star = amin1(top/bot,stoke)

        terma = reynold**0.5 * schmidt**0.333
        termb = reynold**0.5 * schmidt**0.5
        term1 = 4./(reynold*schmidt)*(1. + 0.4*terma + 0.16*termb)
        phi = prtdia/drpdia
        term2 = 4.*phi*(muair/muh2o + (1. + 2.*reynold**0.5)*phi)
        term3 = (stoke - star)/(stoke - star + 2./3.)
        term3 = term3**1.5
        eff = term1 + term2 + term3

        ascav = 4.2e-7*pmm*eff/drpdia
      endif

      return
      end
