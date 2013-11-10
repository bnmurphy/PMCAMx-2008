c=====================================================================
c     EQILIBRIUM PARTITIONING OF CONDENSING OR EVAPORATING ORGANICS
c     - revised for PMCAMx by bkoo (03/09/03)
c     - created by bkoo (04/18/02)
c
c     Predefined pointers:
c     - ICPO1  : pointer of the first condensable gas in the gas array
c     - KAPO1 : pointer of the first SOA species in the aerosol array
c     - KPOM  : pointer of primary organic aerosol in the aerosol array
c     NOTE:
c     - There is currently no code to check consistency in MW of SOA
c       species (gas & aerosol) in block.f & soapdat.f
c=====================================================================
      subroutine eqparto(t,q)

      include 'dynamic.inc'
      include 'soap.com'

      parameter (itermax = 1000) ! comparable to ITMAXEQ in equaer.inc

      real*8  q(ntotal),dq,accom,frq(nsec),frqtot,rlambda,frt,t
      real    ctot(nsoap),caer(nsoap),cgas(nsoap),cpre,mwpre,tempk
      real    convfac


      integer iter
      logical prevsp

      real*8 xsum(nsec),mnk,mxk ! bkoo (11/14/02)
      real   csatT(nsoap)       ! bkoo (11/14/02)

      logical lppm

      if(nsoap.ne.ngas-4)
     &  stop'ERROR: nsoap .ne. ngas-4 in eqparto.f (check dynamic.inc)'

      prevsp = .false. ! currently SOAP do not utilize PREVSP

      ng = nsp*nsecx2                 ! gases
      prod=rgas*temp/(1.01325D5*pres) ! conversion from ppm to umoles/m3

      do i=1,nsoap
        cgas(i) = SNGL(q(ng+ICPO1-1+i) * gmw(ICPO1-1+i) / prod) ! ug/m3
        caer(i) = 0.0
        do j=1,nsecx2
          caer(i) = caer(i) + SNGL(q((j-1)*nsp+KAPO1-1+i))
        enddo
        ctot(i) = cgas(i) + caer(i)
      enddo

      cpre = 0.0
      do i=1,nsecx2
        cpre = cpre + SNGL(q((i-1)*nsp+KPOM)) ! primary OC
      enddo
      mwpre = SNGL(emw(KPOM))
      tempk = SNGL(temp)

      convfac = SNGL(1./prod) ! umol/m3 = ppm * convfac
      lppm = .false. ! gases in ppm if true, ugm3 if false

      call get_param(igrdchm,ichm,jchm,kchm,iout,idiag)

      call soap(nsoap,caer,cgas,tempk,convfac,
     &          iout,igrdchm,ichm,jchm,kchm,lppm,
     &          cpre,mwpre,csatT)                          ! bkoo (11/14/02)

      

      rlambda=0.065d0 ! this should be consistent with that in eqpart.f

      do isec = 1,nsecx2                                   ! bkoo (11/14/02)
        indx=(isec-1)*nsp                                  ! bkoo (11/14/02)
        xsum(isec)=0.0
        if(pflag.eq.1) xsum(isec)=q(indx+KPOM)/emw(KPOM)   ! bkoo (11/14/02)
        do isp = 1,nsoap                                   ! bkoo (11/14/02)
          xsum(isec)=xsum(isec)+q(indx+KAPO1-1+isp)/emw(KAPO1-1+isp) ! bkoo (11/14/02)
        enddo                                              ! bkoo (11/14/02)
        ! give any non-zero value to xsum if it is zero in the section
        if(xsum(isec).eq.0.0) xsum(isec) = tinys
      enddo                                                ! bkoo (11/14/02)

      do isp = 1,nsoap
        ! calculate DQ & ACCOM
        dq = q(ng+ICPO1-1+isp)*gmw(ICPO1-1+isp)/prod - DBLE(cgas(isp)) ! ug/m3
        accom = delta(ICPO1-1+isp)
        ! calculate frq
        frqtot = 0.d0
        mnk = 0.d0 ! bkoo (11/14/02)
        mxk = 0.d0 ! bkoo (11/14/02)
        do isec = 1,nsecx2
          frq(isec) = qn(isec)*
     &        ( q(ng+ICPO1-1+isp)*gmw(ICPO1-1+isp)/prod                   ! bkoo (11/14/02)
     &        - q((isec-1)*nsp+KAPO1-1+isp)/emw(KAPO1-1+isp)/xsum(isec) ! bkoo (11/14/02)
     &        * csatT(isp) )*                                       ! bkoo (11/14/02)
     &        dsec(isec)/(1.d0+rlambda/(accom*dsec(isec)))
c          frqtot = frqtot + frq(isec)
          mnk = min(mnk,frq(isec)) ! bkoo (11/14/02)
          mxk = max(mxk,frq(isec)) ! bkoo (11/14/02)
        enddo

        if(dq.gt.0. .and. mnk.lt.0. .and. mxk.gt.0.) then     ! bkoo (11/14/02)
          do isec = 1,nsecx2                                  ! bkoo (11/14/02)
            frq(isec)=max(frq(isec)-mnk,0.d0)                 ! bkoo (11/14/02)
          enddo                                               ! bkoo (11/14/02)
        elseif(dq.lt.0. .and. mxk.gt.0. .and. mnk.lt.0.) then ! bkoo (11/14/02)
          do isec = 1,nsecx2                                  ! bkoo (11/14/02)
            frq(isec)=min(frq(isec)-mxk,0.d0)                 ! bkoo (11/14/02)
          enddo                                               ! bkoo (11/14/02)
        endif                                                 ! bkoo (11/14/02)
        do isec = 1,nsecx2                                    ! bkoo (11/14/02)
          frqtot = frqtot + frq(isec)                         ! bkoo (11/14/02)
        enddo                                                 ! bkoo (11/14/02)
        ! normalize frq
        do isec = 1,nsecx2
          frq(isec) = frq(isec) / frqtot
        enddo
        ! condense all condensing species
        if(dq.gt.0.d0) then
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            q(indx+KAPO1-1+isp) = q(indx+KAPO1-1+isp) + dq * frq(isec)
          enddo
          ! set the gas species concentration to the equilibrium value
          q(ng+ICPO1-1+isp) = DBLE(cgas(isp)) * prod / gmw(ICPO1-1+isp) ! ppm
        ! evaporate all evaporating species
        elseif(dq.lt.0.d0) then
          iter = 0
 100      frt = 1.d0
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            if(frq(isec).gt.0.d0) frt = MAX( MIN(q(indx+KAPO1-1+isp) ! bkoo (10/07/03)
     &                                  /(-dq*frq(isec)),frt), 0.d0)
          enddo
          frqtot = 0.d0
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            q(indx+KAPO1-1+isp) = 
     &               MAX(q(indx+KAPO1-1+isp)+frt*dq*frq(isec),0.d0)
            if(q(indx+KAPO1-1+isp).lt.tinys) frq(isec) = 0.d0
            frqtot = frqtot + frq(isec)
          enddo
          ! check if we should evaporate more
          dq = (1.d0 - frt) * dq
          if(dq.lt.-1.d-7) then
            if(frqtot.gt.tinys) then ! we have sections which are not empty
              if(iter.le.itermax) then ! check infinite loop
                iter = iter + 1
                do isec = 1,nsecx2
                  frq(isec) = frq(isec) / frqtot
                enddo
                goto 100
              endif
            endif
            ! we need to evaporate more to achieve equilibrium
            ! but we completely evaporate the species in all sections
            ! or exceed itermax
cbk            call errprt('EQUIO-F:DQ    ',t,KAPO1-1+isp,-1,dq,0)
            if(dq.lt.-1.d-3) then                                        ! bkoo_dbg
              write(6,'(A9,E15.5,4I4)')'EQUIO-F: ',dq,isp,ichm,jchm,kchm ! bkoo_dbg
            endif                                                        ! bkoo_dbg
          endif
          ! now set the gas species concentration conservatively
          q(ng+ICPO1-1+isp) = DBLE(ctot(isp)) ! ug/m3
          do isec = 1,nsecx2
            indx=(isec-1)*nsp
            q(ng+ICPO1-1+isp) = q(ng+ICPO1-1+isp) - q(indx+KAPO1-1+isp)
          enddo
          q(ng+ICPO1-1+isp) = q(ng+ICPO1-1+isp) * prod / gmw(ICPO1-1+isp) ! ppm
        endif
      enddo


      return
      end


c=====================================================================
c     This subroutine calculates the organic water absorption using
c     fits from UNIFAC for representative SOA species.
c     - bkoo (04/18/02)
c     Currently not used in PMCAMx
c=====================================================================
      subroutine orgwtr(orgtotal,orgwater)

      include 'dynamic.inc'

      real*8 orgtotal(7),orgwater ! bkoo (09/29/02)
      real*8        soawtr(7,70) ! bkoo (09/29/02)
      integer       iswtr ! bkoo (09/29/02)
      common / soawtr / soawtr,iswtr ! bkoo (09/29/02)
      real*8 wtr
      integer idrh

      if(iswtr.ne.1) then ! bkoo (09/29/02)
        orgwater = 0.d0
        return
      endif

c     modified - bkoo (09/29/02)
      idrh = IDINT(rh*100.d0)
      orgwater = 0.d0
      do i=1,7
        if(idrh.lt.30) then
          wtr = soawtr(i,30) * rh / 0.3d0
        elseif(idrh.ge.99) then
          wtr = soawtr(i,70)
        else
          wtr = soawtr(i,idrh-29)
     &        + (rh*100.d0 - DFLOAT(idrh))
     &        * (soawtr(i,idrh-28) - soawtr(i,idrh-29))
        endif
        orgwater = orgwater + wtr * orgtotal(i) * 1.d-9 ! kg/m3
      enddo

      return
      end


c=====================================================================
c     This subroutine calculates partial pressures of organics at the 
c     surface of the particles.
c     To calculate partial pressure of organics at equilibrium, we use the
c     pseudo-ideal solution theory by Bowman et al. (1997) with parameters
c     taken from Strader et al. (1999) for saturation concentrations.
c     We use the Clausius Clapeyron Equation here to calculate
c     the temperature dependence of the saturation concentration.
c     We divide the code into two subroutines so that the temperature
c     dependence calculations are done only once before DIFFUND.
c     - bkoo (04/18/02)
c     - revised for PMCAMx by bkoo (03/09/03)
c=====================================================================
      subroutine orgps1(q)

      include 'dynamic.inc'
      include 'soap.com'

      real*8  q(ntotal),ctotp(nsoap),ctot,rsum
      real*8  csatT(nsoap)
      real*8  conmin
      integer icont
      integer nsol,idx(nsoap)
      logical ispre
      common / csatx / csatT,ispre,nsol,idx

      if(nsoap.ne.ngas-4)
     &  stop'ERROR: nsoap .ne. ngas-4 in eqparto.f (check dynamic.inc)'

      ispre = .FALSE.
      if(pflag.eq.1) ispre = .TRUE.

      do i=1,nsoap
        csatT(i) = csat(i) * rgas * temp / mwsoap(i) * 1.d-6 ! Pa
        csatT(i) = csatT(i)
     &           * exp((deltah(i)/8.314)*(1/cstemp(i)-1/temp))
      enddo

      ng = nsp*nsecx                  ! gases
      prod=rgas*temp/(1.01325D5*pres) ! conversion from ppm to umoles/m3
      conmin = 0.0001d0               ! should be consistent with that of SOAP
      icont = 0
      rsum  = 0.d0
      do i=1,nsoap
        if(flagsoap(i).eq.1) then ! solution-forming species
          ctot = q(ng+ICPO1-1+i) * gmw(ICPO1-1+i) / prod ! ug/m3
          do k=1,nsecx
            ctot = ctot + q((k-1)*nsp+KAPO1-1+i)
          enddo
          ctotp(i) = ctot * rgas * temp / mwsoap(i) * 1.d-6 ! Pa
          ! as in SOAP
          if(ctot.ge.conmin) then
            icont=icont+1
            idx(icont) = i
            rsum = rsum + ctotp(i) / csatT(i)
          else ! for ignored species flux should be zero
            do k=1,nsecx
              ps(k,ICPO1-1+i) = q(ng+ICPO1-1+i)*1.01325d-1*pres ! Pa
            enddo
          endif
        else ! non-solution-forming species
          do k=1,nsecx
            ps(k,ICPO1-1+i) = csatT(i)
          enddo
        endif
      enddo
      nsol = icont
      if(nsol.eq.0) return
      if(rsum.lt.1.d0 .and. pflag.ne.1) then
        do i=1,nsol
          do k=1,nsecx
            ps(k,ICPO1-1+idx(i)) = ctotp(idx(i))
          enddo
        enddo
        nsol = 0
      endif

      return
      end

      subroutine orgps2(q)

      include 'dynamic.inc'

      real*8  q(ntotal),orgtot
      real*8  csatT(ngas-4)
      integer nsol,idx(ngas-4)
      logical ispre
      common / csatx / csatT,ispre,nsol,idx

      if(nsol.eq.0) return

      do i = 1,nsecx
        orgtot = 0.d0
        if(ispre) orgtot=q((i-1)*nsp+KPOM)/emw(KPOM) ! primary OC
        do k = 1,nsol
          orgtot=orgtot+q((i-1)*nsp+KAPO1-1+idx(k))/emw(KAPO1-1+idx(k))
        enddo
        if(orgtot.lt.tinys) then
          do k = 1,nsol
            ps(i,ICPO1-1+idx(k)) = 0.d0
          enddo
        else
          do k = 1,nsol
            ps(i,ICPO1-1+idx(k)) = dmax1( ( q((i-1)*nsp+KAPO1-1+idx(k))
     &                / emw(KAPO1-1+idx(k)) ) / orgtot * csatT(idx(k)) ! Pa
     &                , 0.d0) ! bkoo (09/29/02)
          enddo
        endif
      enddo

      return
      end

