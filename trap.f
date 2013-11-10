      subroutine trap(rxnrate,radslvr,ratejac,rateslow,dtin,ldark,
     &                H2O,atm,O2,CH4,H2,conc,cncrad,avgrad,tcell,
     &                sddm,nfam,nsen,ddmjac,lddm,nirrrxn,titrt,rrxn_irr,
     &                ldoirr)
c
c-----CAMx v4.02 030709
c
c     TRAP solves the chemical reaction ODEs for a given time step
c     using the Crank-Nicolson method and Euler's method
c     ctmp(nspec+1) is used to store unused species
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        01/22/02    Now has array for to save average of radicals 
c                    over time step
c
c     Input arguments:
c        rxnrate           name of the routine which computes reaction
c                          rates
c        radslvr           name of the routine which computes radical 
c                          concentrations
c        ratejac           name of the routine which computes the
c                          reaction rate and the Jacobian matrix of fast
c                          state species
c        rateslow          name of the routine which computes the
c                          reaction rate of slow state species
c        dtin              time duration to be integrated (hr)
c        ldark             logical flag for darkness
c        H2O               water vapor concentration (ppm)
c        atm               total gas concentration (M) in ppm
c        O2                oxygen concentration (ppm)
c        CH4               methane concentration (ppm)
c        H2                hydrogen concentration (ppm)
c        conc              state species concentrations (ppm)
c        avgrad            average radical species concentrations (ppm)
c        tcell             temperature (K)
c        sddm              sensitivity matrix sddm(nfam,nsen)
c        nfam              number of sensitivity families
c        nsen              number of coefficients in each family
c        ddmjac            name of the routine to calc full Jacobian
c                          (unused elements must be zeroed)
c        lddm              logical flag for DDM sensitivities
c        nirrrxn           number of reactions for IPRM
c        ldoirr            logical flag to calculate IRR data
c
c     Output arguments:
c        conc              state species concentrations (ppm)
c        cncrad            radical species concentrations (ppm)
c        titrt             titration value for IPRM
c        rrxn_irr          array for reaction data for IPRM
c
c     Routines called:
c        RXNRATE
c        RADSLVR
c        RATEJAC
c        CPIVOT
c        RATESLOW
c        DDMCHEM
c
c     Called by:
c        CHEMDRIV
c
      include "camx.prm"
      include "chmstry.com"
      include "filunit.com"
      include "flags.com"
c
      parameter (MXDDM = MXSPEC + MXRADCL)
      parameter (WT = 0.5)
      parameter (FUZZ = 10.0)
c
      logical   ldark, lddm, lzero(MXSPEC), ltest
c
c========================== Process Analysis Begin =============================
c
      integer nirrrxn
      real rrxn_irr(MXRXN)
      real titrt
      logical ldoirr
c
c=========================== Process Analysis End ==============================
c
      dimension conc(MXSPEC+1),rate(MXSPEC+1),res(MXSPEC),rrxn(MXRXN),
     &          rloss(MXSPEC+1),xjac(MXSPEC,MXSPEC),cncrad(MXRADCL),
     &          crold(MXRADCL),cncold(MXSPEC),ctmp(MXSPEC+1),avgrad(MXRADCL)
      dimension y(MXDDM+1), prod(MXDDM), stmp(MXDDM), ipvt(MXDDM),
     &          djac(MXDDM+1,MXDDM+1), amx(MXDDM,MXDDM), 
     &          bmx(MXDDM,MXDDM), cmx(MXDDM,MXDDM),
     &          sddm(nfam,nsen)
      external  rxnrate,radslvr,ratejac,rateslow,ddmjac
c
c-----Entry point
c
c
c======================== DDM Begin =======================
c
c   --- Make sure the sensitivities correctly reflect lower bound ---
c
      if (lddm) then
         do l=1,nsen
           if(conc(l) .LE. FUZZ*bdnl(l) ) then
             do is=1,nfam
               sddm(is,l) = 0.
             enddo
           endif
         enddo
      endif
c
c======================== DDM End =======================
c
c
c========================== Process Analysis Begin =============================
c
      if (ldoirr) then
        do ixn = 1, nirrrxn
          rrxn_irr(ixn) = 0.0
        enddo
      endif
c
c=========================== Process Analysis End ==============================
c
c-----Mass of NO3 and N2O5 goes in NO2 because steady state used 
c 
      conc(kno2) = conc(kno2) + conc(knxoy)
      conc(knxoy) = bdnl(knxoy)
c
c-----Titrate NO against O3 at night
c
      nstrt = 1
      tmp = 0.0
      if(ldark) then
        if (conc(kno).le.bdnl(kno)) then
          nstrt=2
        elseif (conc(ko3).le.bdnl(ko3)) then
          nstrt=1
        elseif (conc(ko3).le.conc(kno)) then
          tmp = conc(ko3) - bdnl(ko3)
          conc(ko3) = bdnl(ko3)
          conc(kno2) = conc(kno2) + tmp
          conc(kno) = amax1((conc(kno)-tmp),bdnl(kno))
          if (conc(kno).le.bdnl(kno)) nstrt=2
        else
          tmp  = conc(kno) - bdnl(kno)
          conc(kno) = bdnl(kno)
          conc(kno2) = conc(kno2) + tmp
          conc(ko3) = amax1((conc(ko3)-tmp),bdnl(ko3))
          nstrt = 2
        endif
      endif
c
c======================== DDM Begin =======================
c
c
c========================== Process Analysis Begin =============================
c
      titrt = tmp
c
c=========================== Process Analysis End ==============================
c
c-----Apply NXOY conservation and NO/O3 titration to DDM coefficients
c
      if (lddm) then
        do k=1,nfam
          sddm(k,kno2) = sddm(k,kno2) + sddm(k,knxoy)
          sddm(k,knxoy) = 0.
          if (ldark) then
            if (nstrt.eq.1) then
              sddm(k,kno) = sddm(k,kno) - sddm(k,ko3)
              sddm(k,kno2) = sddm(k,kno2) + sddm(k,ko3)
              sddm(k,ko3) = 0.
            else
              sddm(k,kno2) = sddm(k,kno2) + sddm(k,kno)
              sddm(k,ko3) = sddm(k,ko3) - sddm(k,kno)
              sddm(k,kno) = 0.
            endif
          endif
        enddo
        do i=1,MXDDM+1
          do j=1,MXDDM+1
            djac(i,j) = 0.
          enddo
        enddo
c
      endif
c
c======================== DDM End   =======================
c
c-----Set concentration of unused species to zero
c
      ctmp(nspec+1) = 0.
c
c-----Set initial dt
c
      if(ldark) then
        dt = amin1(dtin,0.13)
      else
        tmp = conc(kno2) + 5.0*conc(kno)
        tmp = amax1(tmp,0.001)
        dt = (0.07*0.002/(tmp+0.001))+0.03
        dt = amin1(dt,dtin)
      endif
      dtmx = dt
c
c-----Initialize clock
c
      time = 0.
      kountmax=20
c
c-----Save old concentration
c
      do l=1,ngas
        cncold(l) = conc(l)
      enddo
c
c-----Initialize time-averaging of radical concentrations
c
      do l=1,nrad
        avgrad(l) = 0.0
        crold(l) = cncrad(l) !VAK
      enddo
c
c-----Start iteration loop for the solution of the fast state species
c
      ktotal = 0
      kount = 0
  10  kount = kount + 1
c
c-----Test for non-convergence
c
      if (dt.lt.1e-5) then 
        write(*,*) ' No Convergence in TRAP'
        goto 900
      endif
      if (dt.lt.1e-3) kountmax=40
c
c-----Set number of fast species
c
      neq1 = 3
      if (conc(kPAN).gt.0.01*conc(kNO2)) neq1 = nspfst
      if (idmech.eq.5) then
        if (conc(kNPHE).lt.1.0e-5) then
            neq1 = nspfst-1
        endif
      endif
c
c-----Update temporary species concentrations
c
      do l=1,ngas
        ctmp(l) = WT*conc(l) + (1.-WT)*cncold(l)
      enddo
c
c-----Compute reaction rate
c
      call rxnrate(H2O,atm,O2,CH4,H2,cncrad,ctmp,rrxn)
c
c-----Solve for radical species concentrations
c
      call radslvr(ldark,H2O,atm,O2,CH4,H2,cncrad,ctmp,rrxn,crold,dt) !VAK


c
c-----Get rate and Jacobian for fast state species
c
      call ratejac(nstrt,neq1,cncrad,ctmp,rrxn,rate,rloss,
     &             xjac)
      do i=nstrt,neq1
        do j=nstrt,neq1
          xjac(i,j) = WT*dt*xjac(i,j)
        enddo
      enddo
      do l=nstrt,neq1
        res(l) = conc(l) - cncold(l) - dt*rate(l)
        xjac(l,l) = 1. + xjac(l,l)
      enddo
c
c-----Solve for fast state species concentrations
c
      call cpivot(nstrt,neq1,MXSPEC,xjac,res,ierr)
      if (ierr.ne.0) then
         write(*,*) 'Determinant zero in CPIVOT at ', ierr
         goto 900
      endif
c
c-----Check for convergence
c
      rerr = 0.
      errbig = 0.
      istbl = 1
      do l=nstrt,neq1
        conc(l) = conc(l) - res(l)
        if (conc(l).lt.0.0) then
           istbl = -1
        elseif(conc(l).gt.bdnl(l)) then
           rerr = amax1(rerr,abs(res(l)/conc(l)))
           errbig = amax1(errbig,abs(res(l)))
        endif
        conc(l) = amax1(conc(l), bdnl(l))
      enddo
c
c-----If negative concentration or many iterations, reduce time step
c
      if (istbl.lt.0.or.kount.gt.kountmax) then
        dt = dt/2.
        do l=1,neq1
          conc(l) = cncold(l)
        enddo
        kount = 0
        goto 10
      endif
c
c-----Next iteration?
c
      if (kount.le.30) then
        if (rerr.gt.0.001) goto 10
      else
        if (rerr.gt.0.005.and.errbig.gt.1e-8) goto 10
      endif
c
c-----Convergence achieved, update species conentrations.
c
      do l=1,ngas
        ctmp(l) = WT*conc(l) + (1.-WT)*cncold(l)
      enddo
      do l=1,nrad
        avgrad(l) = avgrad(l) + (dt/dtin)*cncrad(l)
        crold(l) = cncrad(l) !VAK
      enddo
c
c-----Compute reaction rate
c
      call rxnrate(H2O,atm,O2,CH4,H2,cncrad,ctmp,rrxn)


c-----------------write OH Tsimpidi------------------------
      call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)
c      do irads = 1,nrads
c	if (l3davg.or.kchm.eq.1) then
c	  print *, 'l3davg = ',l3davg
c	  print *,'Time: ',dtin
c	  print *,'Grid cell: i=',ichm,' j=',jchm,' k=',kchm
c	  print *,'iavg = ',iavg
c	  write(iavg+(kchm-1)*(1+nrads)+1) dtin,ichm,jchm,cncrad(kOH)
c	  print *, 'Writing to file unit: ',iavg+(kchm-1)*(1+nrads)+1-1
c	  write(iavg+(kchm-1)*(1+nrads)+2) dtin,ichm,jchm,cncrad(kNO3)
c	  print *, 'Writing to file unit: ',iavg+(kchm-1)*(1+nrads)+2-1
c	endif  
c      enddo

	if (ichm.eq.55.and.jchm.eq.61) then
	print *,'time = ',time, 'cncrad(OH) = ',cncrad(kOH)
	print *,'time = ',time, 'cncrad(NO3) = ',cncrad(kNO3)
	print *,'time = ',time, 'cncrad(HO2) = ',cncrad(kHO2)
	endif

	bnmradcnc(1,ichm,jchm,kchm) = cncrad(kOH)
	bnmradcnc(2,ichm,jchm,kchm) = cncrad(kNO3)
c	bnmradcnc(3,ichm,jchm,kchm) = cncrad(kN2O5)
	bnmradcnc(3,ichm,jchm,kchm) = cncrad(kHO2)

c----------------------------------------------------------


c
c-----Get rates for slow species
c
      call rateslow(neq1,rrxn,rate)
c
c----------------------------------------------------------------------
c-----Added VOC/NOx dependence for biogenic VOC with O3 
c-----T.Lane (3/22/05)
c     I think nflag is zero at night
c
c      rate(kCG4  )= ( 0.136)*rrxn( 97)+( 0.136)*rrxn( 98)+( 0.136)*
c     &              rrxn(99)*(ctmp(kOLE2)/(ctmp(kNO)+ctmp(kNO2)))
c     &              +( 0.136)*rrxn(100)
c
c
c----------------------------------------------------------------------

c
c======================== DDM Begin =======================
c
c-----Organize information and call DDMCHEM
c     Set species to be solved
c     Radicals formed only by photolysis (e.g., O-atoms) are zero at
c     night so their sensitivities are also zero.  It is assumed that
c     the first photolytic radical is O(1D) and the last is O(3P).
c
      if (lddm) then
        nr=nrad
        if (ldark) nr=nr-ko
        nf=neq1
        ns=ngas-neq1
        lzero(1)=.false.
        if (nstrt.eq.2) lzero(1)=.true.
        do l=2,ngas
          lzero(l)=.false.
        enddo
c
        call ddmchem (sddm,nfam,nsen,nrad,ngas,nr,nf,ns,nr+nf,nf+ns,
     &                nr+nf+ns,nrad+ngas,
     &                nreact,y,prod,stmp,ipvt,djac,amx,bmx,cmx,lzero,
     &                dt,WT,cncrad,ctmp,conc,bdnl,rrxn,ddmjac,
     &                knxoy,kno3,kn2o5,ko1d,ko,ierr,ldark)
c
        if (ierr.ne.0) goto 900
c
      endif
c
c======================== DDM end   =======================
c
c
c-----Update slow species concentrations
c
      do l=neq1+1,ngas
        tmp = cncold(l) + dt*rate(l)
        conc(l)=amax1(tmp,bdnl(l))
      enddo
c
c========================== Process Analysis Begin =============================
c
      if (ldoirr) then
        do ixn = 1, nirrrxn
          rrxn_irr(ixn) = rrxn_irr(ixn) + dt*rrxn(ixn)
        enddo
      endif
c
c=========================== Process Analysis End ==============================
c
c
c-----Report mass of NO3 and N2O5 as NXOY and set flag to zero out NXOY
c     sensitivity if NXOY concentration is reset to NO2 concentration.
c
      conc(knxoy) = cncrad(kno3) + 2.*cncrad(kn2o5)
      ltest = .false.
      if (conc(knxoy).gt.conc(kno2)) then
        conc(knxoy) = conc(kno2)
        ltest = .true.
      endif
c
c======================== DDM Begin =======================
c
      if (lddm) then
        do l=neq1+1,ngas
          if (conc(l).le.FUZZ*bdnl(l)) then
            do k=1,nfam
              sddm(k,l) = 0.
            enddo
          endif
        enddo
      endif
c
c======================== DDM end   =======================
c
c
c-----Update clock, start new time step if there is time left
c
      time = time + dt
      if (time .lt. 0.99*dtin) then
        dt = amin1(dt*1.5, dtmx, dtin-time)
        ktotal = ktotal + kount  
        do l=nstrt,ngas
          cncold(l) = conc(l)
        enddo
        kount = 0
        goto 10
      endif
c
c-----Done time step
c 
c-----Remove mass of NXOY from NO2
c
      conc(kno2) = amax1( 
     &             bdnl(kno2), (conc(kno2) - conc(knxoy)))
c
c======================== DDM Begin =======================
c
c-----Complete NXOY conservation for DDM coefficients
c     NXOY was set equal to NO3 + 2*N2O5 in DDMCHEM
c     If NO2 at lower bound, zero final NO2 sensitivities
c
      if (lddm) then
        do k=1,nfam
          if (ltest) sddm(k,knxoy) = 0.
          if (conc(kno2).le.FUZZ*bdnl(kno2)) then
            sddm(k,kno2) = 0. 
          else
            sddm(k,kno2) = sddm(k,kno2) - sddm(k,knxoy)
          endif
        enddo
      endif
c
c======================== DDM End   =======================
c
      return
c
c-----Error reporting on failure to converge
c
 900  write(iout,'(//,a)') 'ERROR in TRAP:'
      write(iout,'(/,a,i4,1p2e11.3)')
     &    ' No Convergence in TRAP: kount, dt, time = ',kount,dt,time
      write(iout,'(a,1p2e10.3)') ' errbig, rerr = ', errbig,rerr
      write(iout,*) 'ldark, temp(K), water = ',
     &     ldark, tcell, H2O
      write(iout,*) 'M, O2, CH4, H2 = ',
     &     atm, O2, CH4, H2
      write(iout,'(a,4i4)') ' igrd, i, j, k = ', igrdchm,ichm,jchm,kchm 
      write(iout,*) 'No  Name     New Conc  Old Conc   Rerr'
      if (nstrt.ne.1) 
     &  write(iout,800) 1,spname(1),conc(1),cncold(1)
      do l=nstrt,neq1
        write(iout,800) l,spname(l),conc(l),cncold(l),
     &               abs(res(l)/conc(l))
      enddo
      do l=neq1+1,ngas
        write(iout,800) l,spname(l),conc(l),cncold(l)
      enddo
      do l=1,nrad
        write(iout,800) l, nmrad(l), cncrad(l)
      enddo
      write(iout,*) 'No   Rate Con  Rxn Rate'
      do l=1,nreact
        write(iout,810) l, rk(l), rrxn(l)
      enddo
      call camxerr()
 800  format(i3,2x,a7,1p3e10.3)
 810  format(i3,2x,1p3e10.3)
      end


