      subroutine aqrates1(nt, yaq, aqprod, aqdest, yaqprime)
        
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
      include 'math.inc'
      common / radical / irradical      
c      
      real*4 yaq(meqn1),yaqprime(meqn1)
      real*4 aqprod(meqn1), aqdest(meqn1)
      real*4 sodium(nsp_aq)
      real*4 con(28),cmet(4), c(46), gascon(ngas_aq)
      real*4 spres(22), gcon(22)
      real*4 fgl(21), flg(21), gfgl(21), gflg(21), rp(28), rl(28)
      real*4 dp(49), dl(49), rr(120)
      real*4 x(numfunc),par(302)
      real fvec(numfunc),diag(numfunc),fjac(ldfjac,numfunc),r(lr)
      real qtf(numfunc),wa1(numfunc),wa2(numfunc),wa3(numfunc)
      real wa4(numfunc)
      real one      
      real*4 ph
c 
      real hc1, hc2, wl, wlm, tlwc, vlwc, hyd, arytm, hno2, af
      real ah, rsrate, sulfrate, chen, rad, wvol
      real form, tnow, sstemp, crustal, salt, cph
      real pres, qsat, watvap, watcont
      real therm, tmin, temp
       
      integer nfev, info, i, j, n, icount, iprint
           
      common/integration/wat1,relh1,therm,dt,dtemp
      common/aqrate/sodium,pres
      common/gascon/gascon                     ! in ppm
      common /istep/istep
      common/time1/tnow
      common/ph/ph
      common/count/icount,iprint
      common/aqchem1/watcont,watvap,crustal,salt  
                ! lwc,wat. vapor in g/m3
                ! crustal species concentration INSIDE DROPLETS (ug/m3)
                ! Na concentration INSIDE DROPLETS (ug/m3)
      common/tprof/sstemp, sspres
      common/sstate/gcon,con,cmet,rad, wvol, cph
      data one/1.0e0/
      
      external state
c
c     CALCULATION OF THE TEMPERATURE FOR THIS TIME (TEMP)
c      (ASSUMING LINEAR TEMPERATURE CHANGE FOR EACH OPERATOR STEPP)
c
      tmin = tnow                                !in seconds
      temp=therm
      n = numfunc
      
      call qsaturation (temp, qsat)              ! qsat is in ug/m3
c
c     ZERO THE RATE OF CHANGE VECTORS  
c
      do i=1, meqn1 
      yaqprime(i)=0.0
      aqprod(i)=0.0
      aqdest(i)=0.0
      enddo
c
c     SET DUMMY PH TO ZERO FOR PRINTING ONLY 
c
      ph=0.0
c
c     FIND TOTAL LWC  
c
      tlwc = watcont*1.e6                            ! in ug/m3
      vlwc=tlwc*1.e-12                               ! vol/vol
c      
c     CHECK FOR NEGATIVE VALUES
c
csp      do i=1, meqn1
csp      if (yaq(i) .le. 0.) yaq(i)=1.e-20
csp      enddo
c
c     RECONSTRUCT THE MATRICES USING THE AVAILABLE YAQ VALUES
c
c     ** GAS PHASE CONCENTRATIONS (IN PPM) **
c      (SOME ARE ASSUMED TO REMAIN CONSTANT DURING THE AQUEOUS-PHASE
c        CHEMICAL PROCESSES, THE REST ARE UPDATED)
c     NOTE : ALL HCHO AND HCOOH ARE STILL IN THE GAS-PHASE
c      
      spres(1)  = yaq(3)*8.314e-2*temp/(1000.*wso2*pres)   ! SO2 (g)
      spres(2)  = 1.e-20                                   ! H2SO4 (g)
      spres(3)  = gascon(nghno2)                           ! HNO2 (g)
      spres(4)  = 1.e-20                                   ! HNO3 (g) [IT HAS ALREADY DISSOLVED]
      spres(5)  = 350.                                     ! CO2 (g)
      spres(6)  = yaq(4)*8.314e-2*temp/(1000.*wh2o2*pres)  ! H2O2 (g)
      hc1=8.314e-2*temp*vlwc*akhen(7)
      hc2=8.314e-2*temp/(1000.*whcho*pres)
      spres(7)  = yaq(1)*hc2/(hc1+1.)                      ! HCHO (g)
      spres(8)  = yaq(2)*8.314e-2*temp/(1000.*whcooh*pres) ! HCOOH (g)
      spres(9)  = gascon(ngno)                             ! NO (g)
      spres(10) = gascon(ngno2)                            ! NO2 (g)
      spres(11) = gascon(ngo3)                             ! O3 (g)
      spres(12) = gascon(ngpan)                            ! PAN (g)
      spres(13) = 1.e-20                                   ! CH3COOOH (g)
      spres(14) = 1.e-20                                   ! CH3OOH (g)
      spres(15) = 1.e-20                                   ! HCl (g)  [IT HAS ALREADY DISSOLVED]
      spres(16) = gascon(ngoh)                             ! OH (g)
      spres(17) = gascon(ngho2)                            ! HO2 (g)
      spres(18) = gascon(ngno3)                            ! NO3 (g)
      spres(19) = yaq(5)*8.314e-2*temp/(1000.*wnh3*pres)   ! NH3 (g)
      spres(20) = gascon(ngch3o2)                          ! CH3O2 (g)
      spres(21) = 1.e-20                                   ! CH3OH (g)
      spres(22) = watvap*8.314e-2*temp/(1000.*18.*pres)    ! H2O (g)
c
c CALCULATION OF GCON VECTOR (GAS-PHASE CONCENTRATIONS IN MOLE/L)
c
      do 5 i=1,22
      gcon(i) = spres(i)*1.e-6/(0.08206*temp)
5     continue
c
c     ** RADIUS AND LWC FOR THE SECTION
c
      rad = 0.5*avdiam                                ! in m
      wl =  watcont                                   ! in g/m3
      wvol= wl*1.e-6                                  ! in vol/vol
      wlm = wl*1.e6                                   ! in ug/m3
c
c     LOADING OF ALL METAL SPECIES
c     NOTE : Na+ IS AN INPUT. THE REST OF THE METAL SPECIES ARE
c            CALCULATED AS A FRACTION OF THE CRUSTAL AEROSOL MASS.
c
      cmet(1) = firon*crustal*1000./(55.85*wlm)  ! Fe(3+) mol/l
      cmet(2) = fman*crustal*1000./(54.9*wlm)    ! Mn(2+) mol/l
      cmet(3) = salt*1000./(23.*wlm)          ! Na(+)  mol/l
      cmet(4) = caratio*crustal*1000./(40.08*wlm)   ! Ca(2+) mol/l
c
c     DO NOT LET THE Fe(3+) AND Mn(2+) CONCENTRATIONS EXCEED A
c     CERTAIN LIMIT BECAUSE THEY CAUSE EXTREME STIFFNESS DUE TO
c     THE REACTION S(IV)->S(VI) 
c
c      if (cmet(1) .gt. 1.0e-5) cmet(1)=1.0e-5
c      if (cmet(2) .gt. 1.0e-5) cmet(2)=1.0e-5
c
c     LOADING OF THE MAIN AQUEOUS CONCENTRATIONS  (in M)
c
      con(1) = yaq(6)*1000./(wmol(1)*wlm)   ! S(IV)
      if (con(1) .lt. 1.e-20) con(1)=1.e-20
      con(2) = yaq(11)*1000./(wmol(2)*wlm)   ! S(VI)
      con(3) = 0.                                ! N(III) (DETERMINED LATER)
      con(4) = yaq(8)*1000./(wmol(4)*wlm)   ! N(V)
      con(5) = 0.                                ! CO2 (DETERMINED LATER)
      con(6) = yaq(7)*1000./(wmol(6)*wlm)      ! H2O2
      if (con(6) .lt. 1.e-20) con(6)=1.e-20
      con(7) = akhen(7)*spres(7)*1.e-6           ! HCHO
      con(8) = 0.                                ! HCOOH (DETERMINED LATER)
      con(9) = 1.0e-6*akhen(9)*spres(9)          ! NO
      con(10) = 1.0e-6*akhen(10)*spres(10)       ! NO2
      con(11) = 1.0e-6*akhen(11)*spres(11)       ! O3
      con(12) = 1.0e-6*akhen(12)*spres(12)       ! PAN
      con(13) = 1.0e-6*akhen(13)*spres(13)       ! CH3COOOH
      con(14) = 1.0e-6*akhen(14)*spres(14)       ! CH3OOH
      con(15) = yaq(9)*1000./(wmol(15)*wlm) ! HCl
      con(16) = 0.                               ! OH (DETERMINED LATER)
      con(17) = 0.                               ! HO2 (DETERMINED LATER)
      con(18) = 0.                               ! NO3 (DETERMINED LATER)
      con(19) = yaq(10)*1000./(wmol(19)*wlm)     ! NH3
      con(20) = 1.0e-6*akhen(20)*spres(20)       ! CH3O2
      con(21) = 1.0e-6*akhen(21)*spres(21)       ! CH3OH
      con(22) = 0.                               ! Cl (DETERMINED LATER)
      con(23) = 0.                               ! ClOH- (DETERMINED LATER)
      con(24) = 0.                               ! SO4- (DETERMINED LATER)
      con(25) = 0.                               ! SO5- (DETERMINED LATER)
      con(26) = yaq(12)*1000./(wmol(26)*wlm)      ! HSO5-
      con(27) = yaq(13)*1000./(wmol(27)*wlm)       ! HOCH2SO3-
      con(28) = 0.                               ! CO3- (DETERMINED LATER)
c
c     SET A MINIMUM CONCENTRATION (TO AVOID DIVISIONS BY ZERO)
c
      do i=1, 28
cgy      if (con(i) .lt. 1.0e-20) con(i)=1.0e-20
        con(i) = amax1(con(i), 1.0e-20)
      enddo
c
c     CALCULATION OF pH and VOLATILE CONCENTRATIONS (CO2, N(III), HCOOH)
c     (SOLVE THE SYSTEM OF EQUATIONS)
c
      call fullequil(con,spres,cmet,akeq,akhen,vlwc,temp,hyd) 
      ph=-ALOG10(hyd)
c      
      ah = 8.314e-2*temp*vlwc*akhen(3)*(1.+akeq(7)/hyd)/pres
      hno2=spres(3)/(1.+ah)
      con(3)=akhen(3)*(1.+akeq(7)/hyd)*1.e-6*hno2
c
      chen=akhen(5)*(1.+akeq(8)/hyd+akeq(8)*akeq(9)/hyd**2)
      con(5)=chen*spres(5)*1.e-6          ! [CO2 T]aq M
c
      af=8.314e-2*temp*vlwc*akhen(8)*(1.+akeq(13)/hyd)/pres
      form=spres(8)/(1.+af)               ! NEW HCOOH(g) ppm
      con(8)=akhen(8)*(1.+akeq(13)/hyd)*1.e-6*form
c
c     WE CALCULATE THE IONIC SPECIES CONCENTRATIONS
c
      call values(hyd, con, cmet, akeq, c)
c
c     BYPASS THE RADICAL CALCULATION IF NECESSARY
c
      if (iradical .eq. 0) go to 270
c
c     CALCULATE THE CONCENTRATIONS OF THE STEADY-STATE SPECIES
c     PSEUDO-STEADY STATE APPROXIMATION FOR Cl
      call steady(rad,temp,c,con,gcon,akeq,akhen,akre)
c
c     FIRST USE AS INITIAL APPROXIMATION THE PREVIOUS VALUES 
      x(1) = con(16)
      x(2) = con(17)
ckf      x(3) = con(18)
      x(3) = con(23)
      x(4) = con(24)
      x(5) = con(25)
      x(6) = con(28)
c
c     SECOND PREPARE VARIABLES USED (THROUGH COMMON BLOCKS) BY SUBROUTINE STATE
      sstemp = temp
      sspres = pres
      cph = hyd
c     
      DO 51 J = 1, N            
         diag(j) = one
 51   CONTINUE      
c
c     CALL HYBRID
      call hybrd(state,numfunc,x,fvec,xtol,par,maxfev,ml,mu,epsfcn,diag,
     *                 mode,factor,nprint,info,nfev,fjac,ldfjac,r,lr,
     *                 qtf,wa1,wa2,wa3,wa4,302)
c      if (info .ne. 1) write(6,*)'INFO = ', info
     
c     hybrd error control
c        if (info .eq. 0) write(80,*) 'hybrd: improper input parameters'
c        if (info .eq. 1) write(80,*) 'hybrd: rel error problem'
c         if (info .eq. 2) write(80,*) 'hybrd: too many calls to fcn'
c         if (info .eq. 3) write(80,*) 'hybrd: xtol is too small'
c         if (info .eq. 4) write(80,*) 'hybrd: err 4 - poor progress'
c         if (info .eq. 5) write(80,*) 'hybrd: err 5 - poor progress'

c     REPLACE THE STEADY STATE VALUES BACK TO THE CON MATRIX
      con(16)=x(1)
      con(17)=x(2)
ckf      con(18)=x(3)
      con(23)=x(3)
      con(24)=x(4)
      con(25)=x(5)
      con(28)=x(6)
c
c     IF HYBRID HAS NOT CONVERGED SET THE CONCENTRATIONS TO ZERO
c
270   continue
      if (info .eq. 2 .OR. info .eq. 3
     &   .OR. info .eq. 4 .OR. info .eq. 5 .OR. iradical .eq. 0) then
      con(16)=1.e-25
      con(17)=1.e-25
      con(18)=1.e-25
      con(23)=1.e-25
      con(24)=1.e-25
      con(25)=1.e-25
      con(28)=1.e-25
      endif
c
c     PSEUDO-STEADY-STATE APPROXIMATION FOR Cl RADICAL
c     NOT USED BECAUSE IT IS OF SECONDARY IMPORTANCE
csp      call values(hyd, con, cmet, akeq,c)
csp      call react(c,cmet,con,akre,rr,arytm)
csp      pro=rr(23)+rr(49)+rr(96)
csp      if (con(22) .le. 0.0)then
csp      con(22) = 1.e-20
csp      go to 20
csp      endif
csp      destr=(rr(24)+rr(25)+rr(26)+rr(27)+rr(28)+rr(29)+rr(30)+rr(42)+
csp     &  rr(56)+rr(61)+rr(69)+rr(109))/con(22)
csp      if (destr .eq. 0.0)then
csp      con(22)= 1.e-10
csp      go to 20
csp      endif
csp      con(22)=pro/destr
csp 20   continue
c
      call values(hyd, con, cmet, akeq,c)
c
c     FINAL CALCULATION OF REACTION RATES
c
      call react(c,cmet,con,akre,rr,arytm)
c
c     CALCULATION OF NET PRODUCTION AND CONSUMPTION RATES
c
      call addit(rr, arytm, rp, rl)
c
c     CALCULATION OF MASS TRANSFER RATES
c
      call mass(wvol,rad,temp,gcon,con,c,akeq,akhen,fgl,flg,gfgl,gflg)
c
c     CALCULATION OF NET RATES OF CHANGE
c
      call differ(rp,rl,fgl,flg,gfgl,gflg,dp,dl)
c
c     CALCULATION OF RIGHT HAND SIDES OF THE DERIVATIVES  
c
c     ** GAS-PHASE SPECIES **
c          (RATES IN UG/M3 AIR S)
      yaqprime(1) = 1.e9*wvol*wmol(7)*(rp(7)-rl(7))   ! HCHO T 
      aqprod(1) = 1.e9*wvol*wmol(7)*rp(7)
      aqdest(1) = 1.e9*wvol*wmol(7)*rl(7)
c      
      yaqprime(2) = 1.e9*wvol*wmol(8)*(rp(8)-rl(8)) ! HCOOH T
      aqprod(2) = 1.e9*wvol*wmol(8)*rp(8)
      aqdest(2) = 1.e9*wvol*wmol(8)*rl(8)
c
      yaqprime(3) = 1.e9*gmol(1)*(dp(29)-dl(29))    ! SO2(g)
      aqprod(3)   = 1.e9*gmol(1)*dp(29)
      aqdest(3)   = 1.e9*gmol(1)*dl(29)
c
      yaqprime(4) = 1.e9*gmol(6)*(dp(34)-dl(34)) ! H2O2(g)
      aqprod(4)   = 1.e9*gmol(6)*dp(34)
      aqdest(4)   = 1.e9*gmol(6)*dl(34)
c
      yaqprime(5) = 1.e9*gmol(19)*(dp(47)-dl(47)) ! NH3(g)
      aqprod(5)   = 1.e9*gmol(19)*dp(47)
      aqdest(5)   = 1.e9*gmol(19)*dl(47)
c
c     ** AQUEOUS-PHASE SPECIES **
c
      yaqprime(6)= 1.e9*wvol*wmol(1)*(dp(1)-dl(1))   ! S(IV)
      aqprod(6)  = 1.e9*wvol*wmol(1)*dp(1)
      aqdest(6)  = 1.e9*wvol*wmol(1)*dl(1)
c
      yaqprime(7)= 1.e9*wvol*wmol(6)*(dp(6)-dl(6))   ! H2O2
      aqprod(7)  = 1.e9*wvol*wmol(6)*dp(6)
      aqdest(7)  = 1.e9*wvol*wmol(6)*dl(6)
c     
ckf      yaqprime(8)= 1.e9*wvol*wmol(4)*(dp(4)-dl(4))   ! N(V)
ckf      aqprod(8)  = 1.e9*wvol*wmol(4)*dp(4)
ckf      aqdest(8)  = 1.e9*wvol*wmol(4)*dl(4)
      
      yaqprime(8)= 1.e9*wvol*wmol(4)*(rp(4)-rl(4))   ! N(V)
      aqprod(8)  = 1.e9*wvol*wmol(4)*rp(4)
      aqdest(8)  = 1.e9*wvol*wmol(4)*rl(4)
c     
      yaqprime(9)= 1.e9*wvol*wmol(15)*(dp(15)-dl(15)) ! Cl-
      aqprod(9)  = 1.e9*wvol*wmol(15)*dp(15)
      aqdest(9)  = 1.e9*wvol*wmol(15)*dl(15)
c     
      yaqprime(10)= 1.e9*wvol*wmol(19)*(dp(19)-dl(19))  ! NH4+
      aqprod(10)  = 1.e9*wvol*wmol(19)*dp(19)
      aqdest(10)  = 1.e9*wvol*wmol(19)*dl(19)
c     
      yaqprime(11)= 1.e9*wvol*wmol(2)*(dp(2)-dl(2))   ! S(VI)
      aqprod(11)  = 1.e9*wvol*wmol(2)*dp(2)
      aqdest(11)  = 1.e9*wvol*wmol(2)*dl(2)
c
      yaqprime(12)= 1.e9*wvol*wmol(26)*(dp(26)-dl(26)) ! HSO5-
      aqprod(12)  = 1.e9*wvol*wmol(26)*dp(26)
      aqdest(12)  = 1.e9*wvol*wmol(26)*dl(26)
c
      yaqprime(13)= 1.e9*wvol*wmol(27)*(dp(27)-dl(27)) ! HMSA
      aqprod(13)  = 1.e9*wvol*wmol(27)*dp(27)
      aqdest(13)  = 1.e9*wvol*wmol(27)*dl(27)
c      
CKF 
      if (gascon(ngso2) .le. minso2) then
      yaqprime(3) = 0.0
      yaqprime(6) = 0.0
      yaqprime(11) = 0.0
      yaqprime(12) = 0.0
      yaqprime(13) = 0.0
      endif
CKF
c
c     CALCULATION OF APPROPRIATE DESTRUCTION RATE
c
      do 50 i=1, meqn1
         if (yaq(i) .le. 1.e-20) then
         aqdest(i) = 0.0
         go to 50
         endif
         aqdest(i) = aqdest(i)/yaq(i)
 50   continue
c
c     CHANGE TO AVOID DIVISIONS BY ZERO IN INTEGRATION
c
      do i=1,meqn1
cgy      if (aqdest(i) .le. 1.e-18) aqdest(i)=1.e-18
        aqdest(i) = amax1(aqdest(i), 1.0e-18)
      enddo
c
c     CALCULATION OF CHARACTERISTIC TIMES (Used for debugging)
c
cdb      tsm=100.
cdb      do 110 i=1, meqn
c      
cdb      if (aqdest(i) .le. 1.e-10) go to 110
cdb      tchar=1./aqdest(i)
c      
cdb      if (tchar .lt. 0.01)then
cdb      write(80,*)tmin,i,yaq(i),yaqprime(i),tchar
cdb      write(6,*) i,yaq(i),yprime(i)
cdb      endif
c      
cdb      if (tchar .lt. tsm) then
cdb      tsm=tchar
cdb      endif
c     
cdb110   continue
c
c     MASS BALANCE FOR SULFUR
c
      sulfrate=yaqprime(3)/gmol(1)+yaqprime(6)/wmol(1)+
     & yaqprime(11)/wmol(2)
      rsrate=sulfrate/(ABS(yaqprime(3))+ABS(yaqprime(11)) +
     &  ABS(yaqprime(6)))
c      if (ABS(rsrate) .ge. 0.01) then
c        write(80, *)'PROBLEM AT ',tmin/60.
c        write(80, *) yaqprime(3),yaqprime(6),yaqprime(11)
c        write(80, *) rsrate
c        write(80, *)'************************'
c      endif
c
c      icount=icount+1
c      if (icount .ge. iprint)then
c      write(6,120)tmin/60.,yaq(11)       !,ph,rsrate,x(1)*1.e12,yaq(13)
c120   format(1x,2(1x,f8.4))
c       write(79,*)tmin/60.,ph
c
c       PRINTING OF ALL REACTION RATES FOR DEBUGGING
csp     write(3,*)tmin/60.,'****(uM/hr)*******'
csp      do i=1,109
csp      write(3,*)i,rr(i)*1.e6*3600.
cs]      enddo
c      icount=0
c      endif
c
      return
      end
