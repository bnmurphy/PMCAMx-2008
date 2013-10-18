      subroutine fullaero(water,tempk,press,lwc_c,
     &                    mxspec_c,mxrad_c,nspec_c,ngas_c,
     &                    con,cncrad,convfac,t00,dtaer,ichm,
     &                    jchm,kchm,height)
c
c-----PMCAMx v3.01 020531
c
c     FULLAERO  handles the aerosol modules for CAMx
c
c     Modifications:
c        01/30/2002 (tmg); do not set prevsp to true because each grid cell
c                          will be different
c        03/22/2002 (tmg); do not update inorganic gases after soap module
c                          because they are not changed inside it
c        04/15/2002 (tmg); changed molecular weights in conversion of 
c                          organic gases
c        11/25/2002 (tmg); added include file other.inc to track whether
c                          aerosol or aqueous module is called
c        12/05/2002 (tmg); added countcell to track which grid cell is
c                          current (for debugging purposes)
c        03/09/03 (bkoo)   removed code related with lsoap
c                          (SOAP has been merged with inorganic aerosol module)
c
c     Input arguments:
c                          species concentrations (ugr/m3)
c        t0,t1		   begining/ending time for this step
c     Output arguments:
c        q                 species concentrations (ugr/m3)
c
c     Routines called: 
c        aerochem (hybrid,madm,equil)
c        aqchem
c        soap
c 
c     Called by: 
c        CHEMDRIV 
c 
      include 'aerpar.inc'
      include 'dynamic.inc'
      include 'droppar.inc'
      include 'camx_aero.inc'
      include 'dbg.inc'
cbk      include 'other.inc'    ! tmg (11/25/02)
ckf
      include 'camx.prm'
      include 'camx.com'
      include 'grid.com'
      include 'flags.com'
      include 'filunit.com'
      include "ptemiss.com"
      include "bndary.com"
ckf
c
      real*4 con(mxspec_c+1)
      real*4 cncrad(mxrad_c)
      real*4 lv,lwc_c,tempk,ev,es,water,press,convfac,t00
      real*4 dtaer,hour,tmin,t1_min,dt_min
      real*4 t0_min
      real*4 prs,qwatr,rhumid
      real*4 gas(ngas_aq), aerosol(nsect,naers)   
      real*8 t0,t1
      real*8 q(ntotal)
      logical prevsp
c
      integer modeaero
c
c     Variables for RADM
c
      real  r_gas(11),r_aer(9)
      real  pres_pa,dt_sec,cw_kgm3
      real  co2,foa,mhp,paa,nacl,caco3,mgco3,a3fe,b2mn,potcl
      real  r_sum(5),cut(nsect+1)
      real  dnit,dchlo,dammo,dsulf
      real  s_neg,s_pos
      integer     r_idx(4)
c
ckf
      real*4 hgt, klay
      real*4 n2o5nit, nitbef, nitaf, nitbal     
      integer chtype
      common/hgt/klay,hgt
CKF
c
c     Background values for RADM (see aerochem.f)
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
      data eps/0.622/, e0/6.11/, lv/2.5e6/, rv/461./
      data prevsp /.false./
c 
c-----Entry point 
c
c     pass some variables to aerosol module common blocks
c
cbk        countcell=countcell+1    ! tmg (12/05/02)
      naero_c=nspec_c-ngas_c
      if (naero_c .eq. 0 ) then
	 write(*,*) 'No Aerosol Species !!'
	 return
      endif
c
      modeaero = 0
c
      ng = nsec*nsp
      temp=tempk
      pres=dble(press/1013.25)
      prs=press/1013.25

      hour = aint(t00/100.)
      tmin = t00 - 100.*hour
      if (tmin.ge.59.99 .and. tmin.le.60.01) tmin = 60.
      if (tmin.ge.60.0) then
        tmin = tmin - 60.
        hour = hour + 1.
      endif
c
      t1_min = 60. * hour + tmin ! end time (min) for AQCHEM
      dt_min = dtaer * 60.       ! delta t (min) for AQCHEM
      t1 = DBLE(t1_min * 60.)    ! end time (sec) for AEROCHEM
      dt = DBLE(dtaer * 3600.)   ! delta t (sec) for AEROCHEM
      t0 = t1 - dt               ! start time (sec) for AEROCHEM
      tcom = t0                  ! common t (sec) for AEROCHEM
      t0_min = t1_min - dt_min   ! beginning time (min) for AQCHEM
c
c  Calculate RH
c
      qwatr = 1.e-6*water*18./28.8
      ev = qwatr*press/(qwatr+eps)
      es  = e0*exp((lv/rv)*(1./273.-1./tempk))
      rhumid = amin1(0.99,ev/es)
      rh=rhumid
c
      if (lfrst) then
        do i=1,nsecp1
  	  dsec(i)=dsec_c(i)
	  dsecf(i)=dsecf_c(i)
        enddo
c       Need to initialize fdist for first time n2o5 is converted to nitrate
c       in aqueous module.  The fdist distribution is calculated the first
c       time the aqueous phase module is called.  An alternative would be 
c       to move this conversion inside aqueous module.  tmg (04/12/04)
        data fdist /0.0, 0.0, 0.0, 0.0, 0.0, 0.41, 0.38, 0.17, 0.04, 0.0/
      endif
c
      if ( lwc_c.ge.aqcwmin .and. tempk.ge.aqtamin .and. laq ) then
       if ( chaq.eq.'RADM' ) then
        pres_pa = 100. * press
        dt_sec  = dtaer * 3600.
        cw_kgm3 = lwc_c / 1000.
c     gas species mapping
        r_gas(1)  = con(kso2_c)*1.e-6
        r_gas(2)  = con(khno3_c)*1.e-6
        r_gas(3)  = con(knxoy_c)*0.5*1.e-6
        r_gas(4)  = co2*1.e-6
        r_gas(5)  = con(knh3_c)*1.e-6
        r_gas(6)  = con(kh2o2_c)*1.e-6
        r_gas(7)  = con(ko3_c)*1.e-6
        r_gas(8)  = foa*1.e-6
        r_gas(9)  = mhp*1.e-6
        r_gas(10) = paa*1.e-6
        r_gas(11) = con(kh2so4_c)*1.e-6
c     aerosol species mapping
        do n = 1, 5
          r_sum(n) = 0.0
        enddo
        do ksec = 1, nsect
          r_sum(1) = r_sum(1) + con(kpso4_c+(ksec-1))
          r_sum(2) = r_sum(2) + con(kpnh4_c+(ksec-1))
          r_sum(3) = r_sum(3) + con(kpno3_c+(ksec-1))
          r_sum(4) = r_sum(4) + con(kna_c  +(ksec-1))
          r_sum(5) = r_sum(5) + con(kpcl_c +(ksec-1))
        enddo
        r_aer(1)  = (r_sum(1)/96./convfac)*1.e-6
        r_aer(2)  = (r_sum(2)/18./convfac)*1.e-6
        r_aer(3)  = (r_sum(3)/62./convfac)*1.e-6
        r_aer(4)  = (caco3/100./convfac)*1.e-6
        r_aer(5)  = (mgco3/84./convfac)*1.e-6
        if (kna_c.eq.nspec+1) then
          r_aer(6) = (nacl/58./convfac)*1.e-6
        else
          if (r_sum(4).gt.r_sum(5)) then
            r_aer(6) = (r_sum(5)/35./convfac)*1.e-6
          else
            r_aer(6) = (r_sum(4)/23./convfac)*1.e-6
          endif
        endif
        r_aer(7)  = (a3fe/56./convfac)*1.e-6
        r_aer(8)  = (b2mn/55./convfac)*1.e-6
        r_aer(9)  = (potcl/74./convfac)*1.e-6
c     call RADM
        call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)
        call raqchem(tempk,pres_pa,dt_sec,cw_kgm3,r_gas,r_aer,
     &               idiag,iout,igrdchm,ichm,jchm,kchm)
c     map gas back to con
        con(kso2_c)  = amax1(r_gas(1)*1.e6,0.0)    ! SO2 (ppm)
        con(khno3_c) = amax1(r_gas(2)*1.e6,0.0)
        con(knxoy_c) = amax1(r_gas(3)*2.*1.e6,0.0) ! N2O5 gas (ppm)
        con(knh3_c)  = amax1(r_gas(5)*1.e6,0.0)
        con(kh2o2_c) = amax1(r_gas(6)*1.e6,0.0)    ! H2O2 (ppm)
        con(ko3_c)   = amax1(r_gas(7)*1.e6,0.0)    ! O3 (ppm)
        con(kh2so4_c)= amax1(r_gas(11)*1.e6,0.0)
        do n = 1, 9
          r_aer(n) = amax1(r_aer(n),0.0)
        enddo
c     calculate the differences
        dsulf = amax1(r_aer(1)*1.e6*convfac*96.,0.0) - r_sum(1)
        dammo = amax1(r_aer(2)*1.e6*convfac*18.,0.0) - r_sum(2)
        dnit  = amax1(r_aer(3)*1.e6*convfac*62.,0.0) - r_sum(3)
cbk RADM doesn't change NaCl concentrations/CAMx doesn't have KCl
cbk so set dchlo to zero
cbk currently "no NaCl" case not allowed in the full-science aerosol module 
cbk        if (kna_c.eq.nspec+1) then
cbk          dchlo = amax1(r_aer(6)*1.e6*convfac*35.,0.0) - r_sum(5) ! ???
cbk        else
cbk          dchlo = 0.0
cbk        endif
        dchlo = 0.0
c        write(*,*)'RADM-change: ',dsulf,dnit,dammo
c     calculate fdist
        do ksec = 1, nsect+1
          cut(ksec) = SNGL(dsecf(ksec))
        enddo
        call aqdist(nsect,cut,fdist,fdist2)
c     add the differences to the corresponding section
        do ksec = 1, nsect
          con(kpso4_c+(ksec-1))=con(kpso4_c+(ksec-1))+dsulf*fdist(ksec)
          con(kpnh4_c+(ksec-1))=con(kpnh4_c+(ksec-1))+dammo*fdist(ksec)
          con(kpno3_c+(ksec-1))=con(kpno3_c+(ksec-1))+dnit *fdist(ksec)
          con(kpcl_c +(ksec-1))=con(kpcl_c +(ksec-1))+dchlo*fdist(ksec)
        enddo
c     adjust mass to avoid negative concentration sections
        r_idx(1) = kpso4_c
        r_idx(2) = kpnh4_c
        r_idx(3) = kpno3_c
        r_idx(4) = kpcl_c
        do n = 1, 4
          s_neg = 0.0
          s_pos = 0.0
          do ksec = 1, nsect
            if (con(r_idx(n)+(ksec-1)) .lt. 0.0) then
              s_neg = s_neg - con(r_idx(n)+(ksec-1))
              con(r_idx(n)+(ksec-1)) = 0.0
            elseif (con(r_idx(n)+(ksec-1)) .gt. 1.e-12) then
              s_pos = s_pos + con(r_idx(n)+(ksec-1))
            endif
          enddo
          do ksec = 1, nsect
            if (con(r_idx(n)+(ksec-1)) .gt. 1.e-12)
     &          con(r_idx(n)+(ksec-1))
     &        = amax1( (1.-s_neg/s_pos)*con(r_idx(n)+(ksec-1)), 0.0 )
          enddo
        enddo
c     END of RADM
       else
        do k=1,ngas_aq
          gas(k)=0.0
        enddo
        do k=1,naers
        do ks=1,nsect
          aerosol(ks,k)=0.0
        enddo
        enddo
        if (con(kh2so4_c) .gt. 0.0) then
          tot_sulf=0.0
          do k=1,nsect
            tot_sulf=tot_sulf+con(kpso4_c+(k-1))
          enddo
          if (tot_sulf .gt. 0.0) then
            do k=1,nsect
              sfac=con(kpso4_c+(k-1))/tot_sulf
              con(kpso4_c+(k-1))=con(kpso4_c+(k-1))
     &                          +sfac*con(kh2so4_c)*convfac*96.
            enddo
            con(kh2so4_c)=0.0
          endif
        endif
c
c     map con to gas and aerosol gas in ppm; aerosol in ugr/m3
c   
        gas(nga)       = con(knh3_c)         ! NH3(g) in ppm
        gas(ngn)       = con(khno3_c)        ! HNO3(g) in ppm 
        gas(ngc)       = con(khcl_c)         ! HCl (g) in ppm
        gas(ngso2)     = con(kso2_c)         ! SO2(g) in ppm
        gas(ngh2o2)    = con(kh2o2_c)        ! H2O2(g) in ppm
        gas(nghcho)    = con(kform_c)        ! HCHO(g) in ppm    
        gas(nghcooh)   = 0.1*con(kh2o2_c)    ! HCOOH(g) in ppm  = 0.1*H2O2
        gas(nghno2)    = con(khono_c)        ! HNO2(g) in ppm
        gas(ngo3)      = con(ko3_c)          ! O3 in ppm
        gas(ngoh)      = cncrad(koh_c)       ! OH in ppm
        gas(ngho2)     = cncrad(kho2_c)      ! HO2 in ppm
        gas(ngno3)     = cncrad(kno3_c)      ! NO3 in ppm
        gas(ngno)      = con(kno_c)          ! NO(g) in ppm
        gas(ngno2)     = con(kno2_c)         ! NO2(g) in ppm
        gas(ngpan)     = con(kpan_c)         ! PAN(g) in ppm
        gas(ngch3o2h)  = 0.2*con(kh2o2_c)    ! CH3OOH(g) in ppm  = 0.2*H2O2
        gas(ngch3o2)   = 1.0e-6              ! CH3O2(g) in ppm
        gas(ngch3oh)   = 1.0e-3              ! CH3OH(g) in ppm = 1 ppb
        gas(ngch3co3h) = 0.05*con(kh2o2_c)   ! CH3C(O)OOH(g) in ppm  = 0.05*H2O2
c
        do knsec=1,nsect
          aerosol(knsec,naw)    = con(kph2o_c+(knsec-1))   ! water
          aerosol(knsec,naa)    = con(kpnh4_c+(knsec-1))   ! ammonium
          aerosol(knsec,na4)    = con(kpso4_c+(knsec-1))   ! sulfate
          aerosol(knsec,nan)    = con(kpno3_c+(knsec-1))   ! nitrate
          aerosol(knsec,nas)    = con(kna_c+(knsec-1))     ! sodium
          aerosol(knsec,nac)    = con(kpcl_c+(knsec-1))    ! chloride
          aerosol(knsec,nae)    = con(kpec_c+(knsec-1))    ! elemental carbon
          aerosol(knsec,nao)    = con(kpoc_c+(knsec-1))    ! primary organics
          aerosol(knsec,nar)    = con(kcrst_c+(knsec-1))   ! crustal
          aerosol(knsec,nahso5) = 0.0
          aerosol(knsec,nahmsa) = 0.0
        enddo
c
ckf
c        Nitrogen mass balance (convert n2o5 to nitrate)
c
ckf
         n2o5nit = cncrad(3)*prs*101325*62./8.314/tempk
	 cncrad(3) = 0.0
ckf

ckf        
 	 do knsec=1,nsect
	   aerosol(knsec,nan)    = aerosol(knsec,nan)+fdist(knsec)*n2o5nit
	 enddo
	 
ckf
         nitbef = 0.0
	 do i =1, nsect
	 nitbef = nitbef+aerosol(i,nan)*14./62.
	 enddo
	 nitbef = nitbef+(14*prs/(8.314e-5*tempk))*(gas(ngn)+gas(nghno2)+gas(ngno3)+gas(ngno)
     & +gas(ngno2)+gas(ngpan))    
	 nitbef = nitbef+(2*(cncrad(3)*prs*101325*108.011/8.314/tempk)*14./108.01) 
c
c
         call aqchem(gas,aerosol,rhumid,prs,tempk,lwc_c,t0_min,t1_min,
     &            dt_min,ierr,kchm,height,chtype)
CKF     &            dt_min,ierr,kchm,height)
ckf     &   
ckf
         nitaf = 0.0
	 do i =1, nsect
	 nitaf = nitaf+aerosol(i,nan)*14./62.
	 enddo
	 nitaf = nitaf+(14*prs/(8.314e-5*tempk))*(gas(ngn)+gas(nghno2)+gas(ngno3)+gas(ngno)
     & +gas(ngno2)+gas(ngpan))    
	 nitbal = nitaf/nitbef
	 
c	 write(*,*) nitbal
	 
ctmg	 if (nitbal.lt.0.99 .or. nitbal.gt.1.01) then
ctmg         write(iout,*) 'ERROR IN NITROGEN BALANCE'
ctmg	 write(6,*) 'ERROR IN NITROGEN BALANCE'
ctmg         write(iout,*) nitbal, nitbef, nitaf
ctmg	 write(6,*) n2o5nit, nitbal, nitbef, nitaf
ctmg	 endif
ckf
c
c     CAMx doesn't carry HMSA or HSO5 so we put their mass into sulfate (NA4)
c
        do knsec=1,nsect
          aerosol(knsec,na4)=aerosol(knsec,na4)+
     &                       96./111.*aerosol(knsec,nahmsa)+
     &                       96./113.*aerosol(knsec,nahso5)
        enddo
c
c     map gas and aerosol back con,cncrad
c
        con(knh3_c)    = gas(nga)   
        con(khno3_c)   = gas(ngn)   
        con(khcl_c)    = gas(ngc)   
        con(kso2_c)    = gas(ngso2) 
        con(kh2o2_c)   = gas(ngh2o2)
        con(kform_c)   = gas(nghcho)
        con(khono_c)   = gas(nghno2)
        con(ko3_c)     = gas(ngo3)  
        con(kno_c)     = gas(ngno)  
        con(kno2_c)    = gas(ngno2) 
        con(kpan_c)    = gas(ngpan) 
c
        do knsec=1,nsect
          con(kph2o_c+(knsec-1)) = aerosol(knsec,naw)
          con(kpnh4_c+(knsec-1)) = aerosol(knsec,naa)
          con(kpso4_c+(knsec-1)) = aerosol(knsec,na4)
          con(kpno3_c+(knsec-1)) = aerosol(knsec,nan)
          con(kna_c+(knsec-1))   = aerosol(knsec,nas)
          con(kpcl_c+(knsec-1))  = aerosol(knsec,nac)
          con(kpoc_c+(knsec-1))  = aerosol(knsec,nao)
          con(kpec_c+(knsec-1))  = aerosol(knsec,nae)
          con(kcrst_c+(knsec-1)) = aerosol(knsec,nar)
        enddo
       endif
       modeaero = 1
      endif
c     RADM or VSRM -> call AER even if the aqueous module is called
c     OVSR         -> call SOAP alone if the aqueous module is called
      if (laero) then
        if ( chaq.eq.'OVSR' .and. modeaero.eq.1 ) then
          modeaero = 1
        else
          modeaero = 2
        endif
      endif
c
c     if neither aqchem nor aerchem is called we don't call soap
c     - should we??? 
c
      if (modeaero.ne.0) then
        do kk=1,ntotal
          q(kk)=0.d0
        enddo
c
c     map con to q [ugr]  gas in ppm
c
        do knsec=1,nsec
          q((knsec-1)*nsp+kcl+1)=con(kapo1_c+(knsec-1))
          q((knsec-1)*nsp+kcl+2)=con(kapo2_c+(knsec-1))
          q((knsec-1)*nsp+kcl+3)=con(kapo3_c+(knsec-1))
          q((knsec-1)*nsp+kcl+4)=con(kapo4_c+(knsec-1))
          q((knsec-1)*nsp+kcl+5)=con(kapo5_c+(knsec-1))
          q((knsec-1)*nsp+kcl+6)=con(kapo6_c+(knsec-1))
          q((knsec-1)*nsp+kcl+7)=con(kapo7_c+(knsec-1))
          q((knsec-1)*nsp+kcl+8)=con(kapo8_c+(knsec-1))
          q((knsec-1)*nsp+kcl+9)=con(kapo9_c+(knsec-1))
          q((knsec-1)*nsp+kcl+10)=con(kapo10_c+(knsec-1))
          q((knsec-1)*nsp+kcl+11)=con(kaoo1_c+(knsec-1))
          q((knsec-1)*nsp+kcl+12)=con(kaoo2_c+(knsec-1))
          q((knsec-1)*nsp+kcl+13)=con(kaoo3_c+(knsec-1))
          q((knsec-1)*nsp+kcl+14)=con(kaoo4_c+(knsec-1))
          q((knsec-1)*nsp+kcl+15)=con(kaoo5_c+(knsec-1))
          q((knsec-1)*nsp+kcl+16)=con(kaoo6_c+(knsec-1))
          q((knsec-1)*nsp+kcl+17)=con(kaoo7_c+(knsec-1))
          q((knsec-1)*nsp+kcl+18)=con(kaoo8_c+(knsec-1))
          q((knsec-1)*nsp+kcl+19)=con(kaoo9_c+(knsec-1))
          q((knsec-1)*nsp+kcl+20)=con(kaoo10_c+(knsec-1))
          q((knsec-1)*nsp+kcl+21)=con(kabs1_c+(knsec-1))
          q((knsec-1)*nsp+kcl+22)=con(kabs2_c+(knsec-1))
          q((knsec-1)*nsp+kcl+23)=con(kabs3_c+(knsec-1))
          q((knsec-1)*nsp+kcl+24)=con(kabs4_c+(knsec-1))
          q((knsec-1)*nsp+kcl+25)=con(kaas1_c+(knsec-1))
          q((knsec-1)*nsp+kcl+26)=con(kaas2_c+(knsec-1))
          q((knsec-1)*nsp+kcl+27)=con(kaas3_c+(knsec-1))
          q((knsec-1)*nsp+kcl+28)=con(kaas4_c+(knsec-1))
c
          q((knsec-1)*nsp+kso4)=con(kpso4_c+(knsec-1)) / 96.  * 98.
          q((knsec-1)*nsp+kcl)=con(kpcl_c+(knsec-1))   / 35.5 * 36.5
          q((knsec-1)*nsp+kno3)=con(kpno3_c+(knsec-1)) / 62.  * 63.
          q((knsec-1)*nsp+kna)=con(kna_c+(knsec-1))
          q((knsec-1)*nsp+knh4)=con(kpnh4_c+(knsec-1)) / 18.  * 17.
          q((knsec-1)*nsp+kh2o)=con(kph2o_c+(knsec-1))
          q((knsec-1)*nsp+kec)=con(kpec_c+(knsec-1))
          q((knsec-1)*nsp+kpom)=con(kpoc_c+(knsec-1))
          q((knsec-1)*nsp+kcrus)=con(kcrst_c+(knsec-1))
        enddo
c
cbk   SOAP has been merged with inorganic aerosol module - bkoo (03/09/03)
cbk        if (lsoap) then
cbk   	  q(naer+icg1)    =0.d0
cbk   	  q(naer+icg2)    =0.d0
cbk   	  q(naer+icg3)    =0.d0
cbk   	  q(naer+icg4)    =0.d0
cbk        else
c     use MW instead of 100 and actual pressure and temperature for
c     conversion to ppm for organic gases (tmg, 04/15/02)
c     now organic gases are given in ppm - bkoo (08/25/03)
        q(naer+icpo1)   = con(kcpo1_c)
        q(naer+icpo2)   = con(kcpo2_c)
        q(naer+icpo3)   = con(kcpo3_c)
        q(naer+icpo4)   = con(kcpo4_c)
        q(naer+icpo5)   = con(kcpo5_c)
        q(naer+icpo6)   = con(kcpo6_c)
        q(naer+icpo7)   = con(kcpo7_c)
        q(naer+icpo8)   = con(kcpo8_c)
        q(naer+icpo9)   = con(kcpo9_c)
        q(naer+icpo10)   = con(kcpo10_c)
        q(naer+icoo1)   = con(kcoo1_c)
        q(naer+icoo2)   = con(kcoo2_c)
        q(naer+icoo3)   = con(kcoo3_c)
        q(naer+icoo4)   = con(kcoo4_c)
        q(naer+icoo5)   = con(kcoo5_c)
        q(naer+icoo6)   = con(kcoo6_c)
        q(naer+icoo7)   = con(kcoo7_c)
        q(naer+icoo8)   = con(kcoo8_c)
        q(naer+icoo9)   = con(kcoo9_c)
        q(naer+icoo10)   = con(kcoo10_c)
        q(naer+icbs1)   = con(kcbs1_c)
        q(naer+icbs2)   = con(kcbs2_c)
        q(naer+icbs3)   = con(kcbs3_c)
        q(naer+icbs4)   = con(kcbs4_c)
        q(naer+icas1)   = con(kcas1_c)
        q(naer+icas2)   = con(kcas2_c)
        q(naer+icas3)   = con(kcas3_c)
        q(naer+icas4)   = con(kcas4_c)
cbk        endif
c
        q(naer+ih2so4) = con(kh2so4_c)
        q(naer+inh3)   = con(knh3_c)
        q(naer+ihno3)  = con(khno3_c)
        q(naer+ihcl)   = con(khcl_c)
c
        if (modeaero.eq.1) then     ! call SOAP only
          do i=1,nsec
            qt(i)=0.d0
            do ki=2,nsp ! use dry basis
              qt(i)=qt(i)+q((i-1)*nsp+ki) ! total mass per section(ug/m3)
            enddo  
            qn(i)=qt(i)/dsec(i)**3 ! calculate 0th moment(number of particles)
          enddo ! if density = 1g/cm3 units are particles/cm3
          call wdiameter(q) ! wet diameter
          ntotalx  = 0 ! not used
          nsecx    = 0 ! not used
          ntotalx2 = ntotal
          nsecx2   = nsec
          call eqparto(t1,q) ! equilibrium organic aerosol partitioning
          call ddiameter(q) ! dry diameter
          call newdist(t1,q) ! size distribution mapping
          call step(nsec,q) ! calculate water in each section
        else                        ! call SOAP + AER
          call aerchem(chaero,q,t0,t1,lfrst,ierr)
        endif
c
c     map q back to con 
c
        do knsec=1,nsec
          con(kapo1_c+(knsec-1))=q((knsec-1)*nsp+kcl+1)
          con(kapo2_c+(knsec-1))=q((knsec-1)*nsp+kcl+2)
          con(kapo3_c+(knsec-1))=q((knsec-1)*nsp+kcl+3)
          con(kapo4_c+(knsec-1))=q((knsec-1)*nsp+kcl+4)
          con(kapo5_c+(knsec-1))=q((knsec-1)*nsp+kcl+5)
          con(kapo6_c+(knsec-1))=q((knsec-1)*nsp+kcl+6)
          con(kapo7_c+(knsec-1))=q((knsec-1)*nsp+kcl+7)
          con(kapo8_c+(knsec-1))=q((knsec-1)*nsp+kcl+8)
          con(kapo9_c+(knsec-1))=q((knsec-1)*nsp+kcl+9)
          con(kapo10_c+(knsec-1))=q((knsec-1)*nsp+kcl+10)
          con(kaoo1_c+(knsec-1))=q((knsec-1)*nsp+kcl+11)
          con(kaoo2_c+(knsec-1))=q((knsec-1)*nsp+kcl+12)
          con(kaoo3_c+(knsec-1))=q((knsec-1)*nsp+kcl+13)
          con(kaoo4_c+(knsec-1))=q((knsec-1)*nsp+kcl+14)
          con(kaoo5_c+(knsec-1))=q((knsec-1)*nsp+kcl+15)
          con(kaoo6_c+(knsec-1))=q((knsec-1)*nsp+kcl+16)
          con(kaoo7_c+(knsec-1))=q((knsec-1)*nsp+kcl+17)
          con(kaoo8_c+(knsec-1))=q((knsec-1)*nsp+kcl+18)
          con(kaoo9_c+(knsec-1))=q((knsec-1)*nsp+kcl+19)
          con(kaoo10_c+(knsec-1))=q((knsec-1)*nsp+kcl+20)
          con(kabs1_c+(knsec-1))=q((knsec-1)*nsp+kcl+21)
          con(kabs2_c+(knsec-1))=q((knsec-1)*nsp+kcl+22)
          con(kabs3_c+(knsec-1))=q((knsec-1)*nsp+kcl+23)
          con(kabs4_c+(knsec-1))=q((knsec-1)*nsp+kcl+24)
          con(kaas1_c+(knsec-1))=q((knsec-1)*nsp+kcl+25)
          con(kaas2_c+(knsec-1))=q((knsec-1)*nsp+kcl+26)
          con(kaas3_c+(knsec-1))=q((knsec-1)*nsp+kcl+27)
          con(kaas4_c+(knsec-1))=q((knsec-1)*nsp+kcl+28)
c
          con(kpso4_c+(knsec-1))=q((knsec-1)*nsp+kso4) * 96.  / 98.
          con(kpcl_c +(knsec-1))=q((knsec-1)*nsp+kcl)  * 35.5 / 36.5
          con(kpno3_c+(knsec-1))=q((knsec-1)*nsp+kno3) * 62.  / 63.
          con(kna_c  +(knsec-1))=q((knsec-1)*nsp+kna)
          con(kpnh4_c+(knsec-1))=q((knsec-1)*nsp+knh4) * 18.  / 17.
          con(kph2o_c+(knsec-1))=q((knsec-1)*nsp+kh2o)
          con(kcrst_c+(knsec-1))=q((knsec-1)*nsp+kcrus)
          con(kpec_c +(knsec-1))=q((knsec-1)*nsp+kec)
          con(kpoc_c +(knsec-1))=q((knsec-1)*nsp+kpom)
        enddo
c
cbk   SOAP has been merged with inorganic aerosol module - bkoo (03/09/03)
cbk      if (.not.lsoap) then
        con(kcpo1_c) = q(naer+icpo1)
        con(kcpo2_c) = q(naer+icpo2)
        con(kcpo3_c) = q(naer+icpo3)
        con(kcpo4_c) = q(naer+icpo4)
        con(kcpo5_c) = q(naer+icpo5)
        con(kcpo6_c) = q(naer+icpo6)
        con(kcpo7_c) = q(naer+icpo7)
        con(kcpo8_c) = q(naer+icpo8)
        con(kcpo9_c) = q(naer+icpo9)
        con(kcpo10_c) = q(naer+icpo10)
        con(kcoo1_c) = q(naer+icoo1)
        con(kcoo2_c) = q(naer+icoo2)
        con(kcoo3_c) = q(naer+icoo3)
        con(kcoo4_c) = q(naer+icoo4)
        con(kcoo5_c) = q(naer+icoo5)
        con(kcoo6_c) = q(naer+icoo6)
        con(kcoo7_c) = q(naer+icoo7)
        con(kcoo8_c) = q(naer+icoo8)
        con(kcoo9_c) = q(naer+icoo9)
        con(kcoo10_c) = q(naer+icoo10)
        con(kcbs1_c) = q(naer+icbs1)
        con(kcbs2_c) = q(naer+icbs2)
        con(kcbs3_c) = q(naer+icbs3)
        con(kcbs4_c) = q(naer+icbs4)
        con(kcas1_c) = q(naer+icas1)
        con(kcas2_c) = q(naer+icas2)
        con(kcas3_c) = q(naer+icas3)
        con(kcas4_c) = q(naer+icas4)

cbk      endif
        con(kh2so4_c)=q(naer+ih2so4)
        con(knh3_c)  =q(naer+inh3)
        con(khno3_c) =q(naer+ihno3)
        con(khcl_c)  =q(naer+ihcl)
c
      endif
c
c      if ( .not. lsoap .or. .not.lcond ) then
cbk   SOAP has been merged with inorganic aerosol module - bkoo (03/09/03)
cbk      if ( .not. lsoap ) then
      lfrst = .false.
      return
cbk      endif
      end
