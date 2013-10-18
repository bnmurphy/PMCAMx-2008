      subroutine aerochem(h2o,tempk,press,cwc,con,convfac,dt,
     &                    ldoipr,ipa_cel)
c
c-----CAMx v4.02 030709
c
c     AEROCHEM calculates the chemical transformation of:
c        1. Condensible organic gasses to organic aerosol (SOAP)
c        2. Gaseous sulfate to aerosol sulfate (from gas-phase chem)
c        3. SO2 to sulfate via aqueous reactions (RADM-AQ approach)
c        4. Inorganic gas-aerosol equilibrium partitioning for the
c           ammonium/nitrate/sulfate/sodium/chloride system (ISORROPIA)
c     The gas species are in ppm
c     The aerosol species are in (ug/m3)
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications: 
c        1/9/02         Minor bug fixes (units conversions)
c        11/19/02       Incorporated RADM aqueous chemistry and ISORROPIA
c        12/9/02        Incorporated SOAP routine
c
c     Input arguments:
c        h2o                 cell water vapor (ppm)
c        tempk               cell temperature (K)
c        press               cell pressure (mb)
c        cwc                 cloud water content (g/m3)
c        con                 species concentrations (ppm, ug/m3)
c        convfac             conversion factor: umol/m3 = ppm * convfac
c        dt                  time step (hours)
c
c     Output arguments:
c        con                 species concentrations (ppm, ug/m3)
c
c     Routines called: 
c        RAQCHEM
c        ISOROPIA
c        SOAP
c 
c     Called by: 
c        CHEMDRIV 
c
      include 'camx.prm' 
      include 'camx.com' 
      include 'chmstry.com' 
      include 'filunit.com'
c
c========================= Process Analysis Begin ======================
c
      include 'procan.com'
c
      logical ldoipr
      integer ipa_cel
c
c========================== Process Analysis End =======================
c
      real nacl,co2,foa,mhp,paa,caco3,mgco3,a3fe,b2mn,potcl
      real con(MXSPEC+1), cold(MXSPEC)
      real lv,cwc,cw_kgm3
c
c-----Arrays for SOAP
c
      real soa(44),cg(44)
      real csatT(44),cpre,mwpre,mwpre0 ! bkoo (08/29/03)
c
c-----Arrays for RADM aqueous chemistry
c
      real r_gas(11),r_aer(9)
c
c-----Variables for ISOROPIA need double precision
c
      real*8 wi(5),wt(5)
      real*8 rhi,tempi,cntrl(2)
      real*8 gasis(3),aerliq(12),aersld(9),other(6)
      character*15 scasi
c
c-----MW of Primary Organic Aerosol - bkoo (11/13/03)
c     should be consistent with that in /AER/block.f
      data mwpre0 /220.0/
c
c-----Background values if not explicitly modeled
c
      data co2 /330./      ! carbon dioxide, ppm
      data foa /1.e-6/     ! formic acid, ppm
      data mhp /1.e-6/     ! MHP, ppm
      data paa /1.e-6/     ! PAA, ppm
      data nacl /0.05/     ! sea salt, ug/m3
      data caco3 /0./      ! calcium carbonate, ug/m3
      data mgco3 /0./      ! magnesium carbonate, ug/m3
      data a3fe /0.010/    ! Fe+++, ug/m3
      data b2mn /0.005/    ! Mn++, ug/m3
      data potcl /0./      ! potassium chloride, ug/m3
c
c-----Parameters for relative humidity calculation
c
      data eps/0.622/, e0/6.11/, lv/2.5e6/, rv/461./
c 
c-----Entry point 
c
      do ispc = 1, nspec
         cold(ispc) = con(ispc)
      enddo
      con(nspec+1)=0.
c
c-----Calculate relative humidity
c
      qwatr = 1.e-6*h2o*18./28.8 
      ev = qwatr*press/(qwatr + eps) 
      es = e0*exp((lv/rv)*(1./273. - 1./tempk)) 
      rh = 100.0*amin1(1.,ev/es)
c
c-----Partitioning of condensable secondary organics between the gas (CG)
c     and aerosol phases (SOA) using the SOAP semi-volatile scheme.
c     CGs are in ppm, SOAs are in ug/m3, SOAP does the units conversion.
c
      cg(1)  = con(kcpo1)
      cg(2)  = con(kcpo2)
      cg(3)  = con(kcpo3)
      cg(4)  = con(kcpo4)
      cg(5)  = con(kcpo5)
      cg(6)  = con(kcpo6)
      cg(7)  = con(kcpo7)
      cg(8)  = con(kcpo8)
      cg(9)  = con(kcpo9)
      cg(10)  = con(kcpo10)
      cg(11)  = con(kcoo1)
      cg(12)  = con(kcoo2)
      cg(13)  = con(kcoo3)
      cg(14)  = con(kcoo4)
      cg(15)  = con(kcoo5)
      cg(16)  = con(kcoo6)
      cg(17)  = con(kcoo7)
      cg(18)  = con(kcoo8)
      cg(19)  = con(kcoo9)
      cg(20)  = con(kcoo10)
      cg(21)  = con(kcbs1)
      cg(22)  = con(kcbs2)
      cg(23)  = con(kcbs3)
      cg(24)  = con(kcbs4)
      cg(25)  = con(kcas1)
      cg(26)  = con(kcas2)
      cg(27)  = con(kcas3)
      cg(28)  = con(kcas4)
c
      soa(1)  = con(kapo1)
      soa(2)  = con(kapo2)
      soa(3)  = con(kapo3)
      soa(4)  = con(kapo4)
      soa(5)  = con(kapo5)
      soa(6)  = con(kapo6)
      soa(7)  = con(kapo7)
      soa(8)  = con(kapo8)
      soa(9)  = con(kapo9)
      soa(10)  = con(kapo10)
      soa(11)  = con(kaoo1)
      soa(12)  = con(kaoo2)
      soa(13)  = con(kaoo3)
      soa(14)  = con(kaoo4)
      soa(15)  = con(kaoo5)
      soa(16)  = con(kaoo6)
      soa(17)  = con(kaoo7)
      soa(18)  = con(kaoo8)
      soa(19)  = con(kaoo9)
      soa(20)  = con(kaoo10)
      soa(21)  = con(kabs1)
      soa(22)  = con(kabs2)
      soa(23)  = con(kabs3)
      soa(24)  = con(kabs4)
      soa(25)  = con(kaas1)
      soa(26)  = con(kaas2)
      soa(27)  = con(kaas3)
      soa(28)  = con(kaas4)
c
      cpre   = con(kpoa)
      mwpre  = mwpre0

      call soap(28,soa,cg,tempk,convfac,iout,igrdchm,
     &          ichm,jchm,kchm,.true.,cpre,mwpre,csatT) ! bkoo (08/29/03)
c
      con(kcpo1)  = amax1(cg(1),bdnl(kcpo1))
      con(kcpo2)  = amax1(cg(2),bdnl(kcpo2))
      con(kcpo3)  = amax1(cg(3),bdnl(kcpo3))
      con(kcpo4)  = amax1(cg(4),bdnl(kcpo4))
      con(kcpo5)  = amax1(cg(5),bdnl(kcpo5))
      con(kcpo6)  = amax1(cg(6),bdnl(kcpo6))
      con(kcpo7)  = amax1(cg(7),bdnl(kcpo7))
      con(kcpo8)  = amax1(cg(8),bdnl(kcpo8))
      con(kcpo9)  = amax1(cg(9),bdnl(kcpo9))
      con(kcpo10)  = amax1(cg(10),bdnl(kcpo10))
      con(kcoo1)  = amax1(cg(11),bdnl(kcoo1))
      con(kcoo2)  = amax1(cg(12),bdnl(kcoo2))
      con(kcoo3)  = amax1(cg(13),bdnl(kcoo3))
      con(kcoo4)  = amax1(cg(14),bdnl(kcoo4))
      con(kcoo5)  = amax1(cg(15),bdnl(kcoo5))
      con(kcoo6)  = amax1(cg(16),bdnl(kcoo6))
      con(kcoo7)  = amax1(cg(17),bdnl(kcoo7))
      con(kcoo8)  = amax1(cg(18),bdnl(kcoo8))
      con(kcoo9)  = amax1(cg(19),bdnl(kcoo9))
      con(kcoo10)  = amax1(cg(20),bdnl(kcoo10))
      con(kcbs1)  = amax1(cg(21),bdnl(kcbs1))
      con(kcbs2)  = amax1(cg(22),bdnl(kcbs2))
      con(kcbs3)  = amax1(cg(23),bdnl(kcbs3))
      con(kcbs4)  = amax1(cg(24),bdnl(kcbs4))
      con(kcas1)  = amax1(cg(25),bdnl(kcas1))
      con(kcas2)  = amax1(cg(26),bdnl(kcas2))
      con(kcas3)  = amax1(cg(27),bdnl(kcas3))
      con(kcas4)  = amax1(cg(28),bdnl(kcas4))
c
      con(kapo1)  = amax1(soa(1),bdnl(kapo1))
      con(kapo2)  = amax1(soa(2),bdnl(kapo2))
      con(kapo3)  = amax1(soa(3),bdnl(kapo3))
      con(kapo4)  = amax1(soa(4),bdnl(kapo4))
      con(kapo5)  = amax1(soa(5),bdnl(kapo5))
      con(kapo6)  = amax1(soa(6),bdnl(kapo6))
      con(kapo7)  = amax1(soa(7),bdnl(kapo7))
      con(kapo8)  = amax1(soa(8),bdnl(kapo8))
      con(kapo9)  = amax1(soa(9),bdnl(kapo9))
      con(kapo10)  = amax1(soa(10),bdnl(kapo10))
      con(kaoo1)  = amax1(soa(11),bdnl(kaoo1))
      con(kaoo2)  = amax1(soa(12),bdnl(kaoo2))
      con(kaoo3)  = amax1(soa(13),bdnl(kaoo3))
      con(kaoo4)  = amax1(soa(14),bdnl(kaoo4))
      con(kaoo5)  = amax1(soa(15),bdnl(kaoo5))
      con(kaoo6)  = amax1(soa(16),bdnl(kaoo6))
      con(kaoo7)  = amax1(soa(17),bdnl(kaoo7))
      con(kaoo8)  = amax1(soa(18),bdnl(kaoo8))
      con(kaoo9)  = amax1(soa(19),bdnl(kaoo9))
      con(kaoo10)  = amax1(soa(20),bdnl(kaoo10))
      con(kabs1)  = amax1(soa(21),bdnl(kabs1))
      con(kabs2)  = amax1(soa(22),bdnl(kabs2))
      con(kabs3)  = amax1(soa(23),bdnl(kabs3))
      con(kabs4)  = amax1(soa(24),bdnl(kabs4))
      con(kaas1)  = amax1(soa(25),bdnl(kaas1))
      con(kaas2)  = amax1(soa(26),bdnl(kaas2))
      con(kaas3)  = amax1(soa(27),bdnl(kaas3))
      con(kaas4)  = amax1(soa(28),bdnl(kaas4))
c
c-----Do RADM aqueous chemistry if CWC is above threshold
c     all conc units must be mol/mol (mixing ratio)
c
      pres_pa = 100.*press
      dt_sec = dt*3600.
      cw_kgm3 = cwc/1000.

      r_gas(1)  = con(kso2)*1.e-6
      r_gas(2)  = con(khno3)*1.e-6
      r_gas(3)  = con(knxoy)*0.5*1.e-6
      r_gas(4)  = co2*1.e-6
      r_gas(5)  = con(knh3)*1.e-6
      r_gas(6)  = con(kh2o2)*1.e-6
      r_gas(7)  = con(ko3)*1.e-6
      r_gas(8)  = foa*1.e-6
      r_gas(9)  = mhp*1.e-6
      r_gas(10) = paa*1.e-6
      r_gas(11) = con(ksulf)*1.e-6

      r_aer(1)  = (con(kpso4)/96./convfac)*1.e-6
      r_aer(2)  = (con(kpnh4)/18./convfac)*1.e-6
      r_aer(3)  = (con(kpno3)/62./convfac)*1.e-6
      r_aer(4)  = (caco3/100./convfac)*1.e-6
      r_aer(5)  = (mgco3/84./convfac)*1.e-6
      if (kna.eq.nspec+1) then
        r_aer(6) = (nacl/58./convfac)*1.e-6
      else
        if (con(kna).gt.con(kpcl)) then
          r_aer(6) = (con(kpcl)/35./convfac)*1.e-6
        else
          r_aer(6) = (con(kna)/23./convfac)*1.e-6
        endif
      endif
      r_aer(7)  = (a3fe/56./convfac)*1.e-6
      r_aer(8)  = (b2mn/55./convfac)*1.e-6
      r_aer(9)  = (potcl/74./convfac)*1.e-6

      if (cwc.ge.cwmin .and. tempk.ge.tamin) then
        call raqchem(tempk,pres_pa,dt_sec,cw_kgm3,r_gas,r_aer,
     &               idiag,iout,igrdchm,ichm,jchm,kchm)

        con(kso2)  = amax1(r_gas(1)*1.e6,bdnl(kso2))     ! SO2 (ppm)
        con(knxoy) = amax1(r_gas(3)*2.*1.e6,bdnl(knxoy)) ! N2O5 gas (ppm)
        con(kh2o2) = amax1(r_gas(6)*1.e6,bdnl(kh2o2))    ! H2O2 (ppm)
        con(ko3)   = amax1(r_gas(7)*1.e6,bdnl(ko3))      ! O3 (ppm)
        do ispc = 1,9
          r_aer(ispc) = amax1(r_aer(ispc),0.0)
        enddo
      endif
c
c-----Inorganic aerosol equilibrium chemistry with ISOROPIA
c     convert conc units to mol/m3 (double precision)
c
cbk      rhi = rh/100.
      rhi = amin1( rh/100., 0.994 )
      tempi = tempk
      cntrl(1) = 0.d0                  ! 0 = forward problem
      if (cwc.ge.cwmin .and. tempk.ge.tamin) then
        cntrl(2) = 1.d0                ! 1 = metastable (liquid only)
      else
        cntrl(2) = 0.d0                ! 0 = solids and liquid allowed
      endif
      if (kna.eq.nspec+1) then
        wi(1) = r_aer(6)*convfac             ! total sodium
        wi(5) = wi(1)                        ! total chloride
      else
        wi(1) = con(kna)/23./1.e6            ! total sodium
        wi(5) = con(kpcl)/35./1.e6 + 
     &          con(khcl)*convfac/1.e6 +
     &          r_aer(9)*convfac             ! total chloride
      endif
      wi(2) = (r_gas(11) + r_aer(1))*convfac ! total sulfate
      wi(3) = (r_gas(5) + r_aer(2))*convfac  ! total ammonium
      wi(4) = (r_gas(2) + r_aer(3))*convfac  ! total nitrate

      call isoropia(wi,rhi,tempi,cntrl,wt,gasis,aerliq,aersld,
     &              scasi,other)
c 
c-----Load results back to local CON array (ppm for gas, ug/m3 for aerosol)
c 
      con(ksulf) = bdnl(ksulf)                     ! sulfuric acid gas (ppm)
c
      con(kpso4) = wt(2)*96.*1.e6                  ! sulfate aerosol (ug/m3)
cbk      con(kpso4) = amax1(con(kpso4)-bdnl(ksulf),bdnl(kpso4))
      con(kpso4) = amax1(con(kpso4)-con(ksulf)*convfac*96.,bdnl(kpso4))
c
      con(khno3) = gasis(2)*1.e6/convfac           ! nitric acid gas (ppm)
      con(khno3) = amax1(con(khno3),bdnl(khno3))    
c
      con(kpno3) = (aerliq(7) + aerliq(11) +
     &              aersld(1) + aersld(2))*62.*1.e6! nitrate aerosol (ug/m3)
      con(kpno3) = amax1(con(kpno3),bdnl(kpno3))    
c
      con(knh3 ) = gasis(1)*1.e6/convfac           ! ammonia gas (ppm)
      con(knh3 ) = amax1(con(knh3 ),bdnl(knh3 ))    
c
      con(kpnh4) = (aerliq(3) + aerliq(9) +
     &              aersld(2) + aersld(4) + 
     &              2.*aersld(6) + aersld(8) + 
     &              3.*aersld(9))*18.e6            ! ammonium aerosol (ug/m3)
      con(kpnh4) = amax1(con(kpnh4),bdnl(kpnh4))    
c
      if (kna.ne.nspec+1) then
c        con(kna) = wt(1)*23.*1.e6                  ! sodium aerosol (ug/m3)
        con(kna) = amax1(con(kna),bdnl(kna))
        con(khcl) = gasis(3)*1.e6/convfac          ! chlorine gas (ppm)
        con(khcl) = amax1(con(khcl),bdnl(khcl))
        con(kpcl) = (aerliq(4) + aerliq(10) +
     &              aersld(3) + aersld(4))*35.*1.e6! chloride aerosol (ug/m3)
        con(kpcl) = amax1(con(kpcl),bdnl(kpcl))
      endif
c
      if (kph2o.ne.nspec+1) then
        con(kph2o) = aerliq(8)*18.*1.e6
        con(kph2o) = amax1( con(kph2o), bdnl(kph2o) ) ! aerosol water (ug/m3)
      endif
c
c-----Check the lower bounds
c
      do is = 1, nspec
        if (con(is).lt.0.) then
          write(iout,'(//,A)') ' ERROR in AEROCHEM:'
          write(iout,'(/,A,i3)') 
     &              ' Negative concentration for species', is
          write(iout,'(/,a,4i4)') ' igrd, i, j, k = ',
     &                igrdchm,ichm,jchm,kchm
          write(iout,'(/,a,/,4f10.2)')
     &           '  H2O(ppm)  Temp(K)  Press(mBar) CWC(g/m3)',
     &                h2o, tempk, press, cwc
          write(iout,'(//,A)') ' Concentrations are:'
          write(iout,'(A5,A7,A12,A12)') 
     &           ' No.','name','C final','C input'
          do l=1,nspec
            write(iout,'(i3,2x,a7,1p2e12.3)')
     &           l, spname(l), con(l), cold(l)
          enddo
          call camxerr()
        endif
      enddo
c
c========================= Process Analysis Begin ======================
c
      if( ldoipr ) then
         do ispc = 1, nspec
            cipr (IPR_FAERO, ipa_cel, ispc) =
     &         cipr (IPR_FAERO, ipa_cel, ispc) + con(ispc)-cold(ispc)
         enddo
      endif
c
c========================== Process Analysis End =======================
c
      return
c
      end
