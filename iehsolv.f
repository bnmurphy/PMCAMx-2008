      subroutine iehsolv(ierxn,ierate,iejac,ieslow,dtin,
     &                   H2O,atm,O2,CH4,H2,conc,cncrad,ldark,tcell,
     &                   nirrrxn,rrxn_irr,ldoirr)
c
c-----CAMx v4.02 030709
c
c     IEHSOLV is the driver for the implicit-explicit hybrid chemistry solver.
c     Steady state for some radicals, LSODE for fast species, first order for
c     slow state species. 
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c
c        ierxn             name of the routine which computes the
c                          reaction rates
c        ierate            name of the routine which computes the
c                          species rates for lsode
c        iejac             name of the routine which computes the
c                          Jacobian for lsode
c        ieslow            name of the routine which computes the
c                          slow species rates
c        dtin              time duration to be integrated (hr)
c        H2O               water vapor concentration (ppm)
c        atm               total gas concentration (M) in ppm
c        O2                oxygen concentration (ppm)
c        CH4               methane concentration (ppm)
c        H2                hydrogen concentration (ppm)
c        conc              state species concentrations (ppm)
c        cncrad            radical species concentrations (ppm)
c        ldark             darkness flag, true = dark
c        tcell             temperature of this cell (K)
c        nirrrxn           number of reactions for IPRM
c        ldoirr            flag to calculate IRR data
c
c     Output arguments:
c        conc              state species concentrations (ppm)
c        cncrad            radical species concentrations (ppm)
c        rrxn_irr          array for reaction data for IPRM
c
c     Routines called:
c        LSODE
c
c     Called by:
c        CHEMDRIV
c        RADDRIVR
c
      include "camx.prm"
      include "chmstry.com"
      include "filunit.com"
      include "iehchem.com"
c
c========================== Process Analysis Begin ==================================
c
      real rrxn_old(MXRXN), rrxn_irr(MXRXN)
      logical ldoirr
c
c=========================== Process Analysis End ===================================
c
c --- preserve lsode common block between repeat calls
c
      common /ls0001/ rls, ils
c$omp threadprivate(/ls0001/)
c
c --- parameters for LSODE
c
      parameter (LSDIM = MXSPEC + MXRADCL + 1)
      parameter (LIW = LSDIM + 20)
      parameter (LRW = LSDIM*LSDIM + 9*LSDIM + 22)
      parameter (MF=21, ITOLER=1)
      parameter (ITASK=1, IOPT=1)
      parameter (DTIEH = 0.025)
c
      integer iwork(LIW), neq(3)
      integer ils(39)
      real conc(MXSPEC+1), cncrad(MXRADCL)
      real conc0(MXSPEC+1), cncrad0(MXRADCL)
      real y(LSDIM+5), rwork(LRW)
      real rrxn(MXRXN), rate(LSDIM), gain(LSDIM), loss(LSDIM)
      real rls(218), rtol, atol
      logical ldark
c
      external ierxn, ierate, iejac, ieslow
c
c --- Entry point
c     Tolerance Declaration
      rtol = 5.e-5
      atol = 1.e-6
101   continue
c
c
c========================== Process Analysis Begin ==================================
c
      if (ldoirr) then
        do irxn = 1, MXRXN
           rrxn_old(irxn) = 0.0
           rrxn_irr(irxn) = 0.0
        enddo
      endif
c
c=========================== Process Analysis End ===================================
c
c --- zero the working array elements that can pass options to LSODE
c     then set the max iterations to 1000
c
      do i=5,10
        rwork(i) = 0.
        iwork(i) = 0
      enddo
      iwork(6) = 1000
c
c --- conserve nitrogen within the NXOY scheme
c
      scl = conc(kNXOY) / ( cncrad(kNO3) + 2.0*cncrad(kN2O5) )
      cncrad(kNO3)  = scl * cncrad(kNO3)
      cncrad(kN2O5) = scl * cncrad(kN2O5)
c
c --- load concentrations for IEH solver
c
      j=0
      do i=iessrad+1, nrad
        j=j+1
        y(j)=cncrad(i)
        cncrad0(i)=cncrad(i)
      enddo
      do i=1, nspfst
        j=j+1
        y(j)=conc(i)
        conc0(i)=conc(i)
      enddo
      neq(1) = j
      do i=1,iessrad
        j=j+1
        y(j)=cncrad(i)
        cncrad0(i)=cncrad(i)
      enddo
      do i=nspfst+1, ngas
        j=j+1
        y(j)=conc(i)
        conc0(i)=conc(i)
      enddo
c
c --- pass some values to ierate and iejac
c
      neq(2) = nrad+ngas
      neq(3) = nreact
      y(neq(2)+1) = 0.
      y(neq(2)+2) = H2O
      y(neq(2)+3) = atm
      y(neq(2)+4) = O2
      y(neq(2)+5) = CH4 
      y(neq(2)+6) = H2
c
c --- initialize LSODE
c
      istate = 1 
c
c---  initialize clock
c
      t = 0.
      nstep = int((dtin+0.0001)/DTIEH) 
      nstep = max(nstep,1)
      dtuse = dtin/nstep
c
c --- integrate in steps
c
      do istep = 1, nstep
        tout = istep*dtuse
        rwork(1) = tout
c
c========================== Process Analysis Begin ==================================
c
        if (ldoirr) then
          call ierxn(y,neq(2),rrxn,rk,nreact)
          do irxn = 1, nirrrxn
            rrxn_old(irxn) = rrxn(irxn)
          enddo
        endif
c
c=========================== Process Analysis End ===================================
c
c --- solve radicals and fast species
c
        call lsode(ierate,neq,y,t,tout,ITOLER,rtol,atol,ITASK, 
     &             istate,IOPT,rwork,LRW,iwork,LIW,iejac,MF,
     &             rk,rrxn,gain,loss)
c
c --- check for errors from LSODE
c
        if (istate.eq.-1.or.istate.eq.-4.or.istate.eq.-5) then
           rtol = rtol * 0.1
           atol = atol * 0.1
           write(iout,103) istate
103        format(/,'Information message from the IEH chemistry solver'/
     +              'LSODE reports istate = ',i3, ' in cell:')
           write(iout,'(a,4i4)') 'igrd, i, j, k = ', 
     +          igrdchm,ichm,jchm,kchm 
           write(iout,105) rtol*10,rtol,atol*10,atol
105        format('The tolerances for LSODE will be reduced to'/
     +          'Relative Tolerance: ',e11.4,'(old) ',e11.4,'(new)'/
     +          'Absoulte Tolerance: ',e11.4,'(old) ',e11.4,'(new)'/
     +          'and the calculation will be repeated.')
           go to 101 
        endif

        if (istate.eq.-2.or.istate.eq.-3.or.istate.eq.-6) go to 700 
c
c --- update slow species concentrations
c
        do l = 1, nrad+nspfst
          y(l) = amax1(y(l),0.0)
        enddo
        call ierxn(y,neq(2),rrxn,rk,nreact)
        call ieslow(rrxn,rate,gain,loss,nreact,neq(2),
     &              nrad+nspfst+1,nrad+ngas)


        do l = nrad+nspfst+1, nrad+ngas
          tmp = y(l) + dtuse*rate(l)
          y(l) = amax1(tmp, bdnl(l-nrad))
        enddo
c
c========================== Process Analysis Begin ==================================
c
        if (ldoirr) then
          do irxn = 1, nirrrxn
            rrxn_irr(irxn) = rrxn_irr(irxn) +
     &                       0.5*dtuse*(rrxn(irxn)+rrxn_old(irxn))
          enddo
        endif
c
c=========================== Process Analysis End ===================================
c
      enddo
c
c --- save concentrations
c
 700  j=0
      do i=iessrad+1, nrad
        j=j+1
        cncrad(i) = y(j)
      enddo
      do i=1, nspfst
        j=j+1
        conc(i) = y(j)
      enddo
      do i=1,iessrad
        j=j+1
        cncrad(i) = y(j)
      enddo
      do i=nspfst+1, ngas
        j=j+1
        conc(i) = y(j)
      enddo
      if (istate .lt. 0) go to 900
c
c --- load NXOY scheme
c
      conc(kNXOY) = cncrad(kNO3) + 2.0*cncrad(kN2O5)
c
c --- check state species lower bounds
c
      do i=1, ngas
        conc(i) = amax1(conc(i),bdnl(i))
      enddo
c
      return
c
c --- Error handling
c
 900  write(iout,*)
      write(iout,*) ' ERROR reported in IEHSOLV by the LSODE solver'
      write(iout,*) '  LSODE reporting ISTATE = ', istate
      write(iout,*) '  look up this error code in the LSODE source code'
      write(iout,901) rwork(11),rwork(12),rwork(13),dtin
      write(iout,903) iwork(12),iwork(13),iwork(11)
      write(iout,*)
      write(iout,*) 'ldark, temp(K), water = ',
     &            ldark, tcell, H2O
      write(iout,*) 'M, O2, CH4, H2  = ',
     &            atm, O2, CH4, H2
      write(iout,'(a,4i4)') ' igrd, i, j, k = ', igrdchm,ichm,jchm,kchm 
      write(iout,*) 'No  Name    New Conc  Init Conc '
      do l=1,ngas
        write(iout,905) l, spname(l), conc(l), conc0(l)
      enddo
      do l=1,nrad
        write(iout,905) l, nmrad(l), cncrad(l), cncrad0(l)
      enddo
      write(iout,*) 'No   Rate Con  Rxn Rate'
      do l=1,nreact
        write(iout,910) l, rk(l), rrxn(l)
      enddo
      write(iout,*)
c
      call camxerr()
c
 905  format(i3,2x,a7,1p3e10.3)
 910  format(i3,2x,1p3e10.3)
 901  format('   last step size  = ', e12.3,/
     &       '   next step size  = ', e12.3,/
     &       '   current time    = ', e12.3,/
     &       '   total time step = ', e12.3)
 903  format('   # of steps ',i4,' # of f-s ',i4,' # of j-s ',i4) 
c
      end 
