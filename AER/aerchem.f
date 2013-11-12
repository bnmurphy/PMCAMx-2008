c=====================================================================+
c MAY 2000: modified by bkoo                                          +
c           - common test program for EQUI, MADM, and HYBR            +
c           - modified initialization step                            +
c 1/25/02 (tmg)  use wet or dry basis for diameter as needed 
c 1/31/02 (tmg)  call step before returning
c 4/15/02 (tmg)  do not call newdist if soap module will be called next
c 10/22/02 (tmg) added parameter to subroutine step (# of sections)
c 03/07/03 (bkoo) commented out initial call to step (not needed)
c 03/09/03 (bkoo) call newdist regardless of lsoap
c                 (SOAP has been merged with inorganic aerosol module)
c=====================================================================+
      subroutine aerchem(chaero,q,t0,t1,lfrst,ierr)
c
c JUNE 1999 - modified by bkoo (linear interpolation)
c
      include 'dynamic.inc'
      include 'dbg.inc'

      real*8 q(ntotal),qtot1(nexti),qtot2(nexti)
      character*4 chaero
      logical lfrst
      real*4 qmass(nsp)

      real*8 save1(ntotal), save2(nsec), save3(nsec), save4(nexti)  ! bkoo

      common / counters / icount, ifail  ! bkoo (01/08/01)
      common / debug / icnt

      icount = 0  ! tmg (01/21/02)
      icnt = 1
      ierr = 0
c   save initial total ion concentrations for ion balance
c	if ( lfrst ) then ! always do since each cell different (tmg, 01/24/02)
      prod=rgas*temp/(1.01325d5*pres)	! pres (bkoo, 6/9/00)
      qtot0(kH2O)= 0.d0			!gas phase mass
      qtot0(kSO4)= q(naer+ih2so4)/(prod/gmw(ih2so4))  
      qtot0(kNa) = 0.d0  
      qtot0(kNO3)= q(naer+ihno3) /(prod/gmw(ihno3))  
      qtot0(kNH4)= q(naer+inh3)  /(prod/gmw(inh3))  
      qtot0(kCl) = q(naer+ihcl)  /(prod/gmw(ihcl))  
      do ki=1,nexti                 
        do i=1,nsec
          qtot0(ki)=qtot0(ki)+q((i-1)*nsp+ki)	!aerosol mass
        enddo
        save4(ki)=qtot0(ki) ! bkoo
      enddo

      do i=1,nsec
        qt(i)=0.d0
        do ki=2,nsp ! use dry basis first, tmg (01/23/02)
          qt(i)=qt(i)+q((i-1)*nsp+ki) ! total mass per section(ug/m3)
        enddo  
        qn(i)=qt(i)/dsec(i)**3 ! calculate 0th moment(number of particles)
        save2(i)=dsec(i)
        save3(i)=qn(i) ! bkoo
      enddo ! if density = 1g/cm3 units are particles/cm3

      do  i=1,nsp
        qmass(i)=0.
        do k=1,nsec
          qmass(i)=qmass(i)+q((k-1)*nsp+i)
        enddo
      enddo

c      Make sure that the ratio of sulfate to sodium is at least that of 
c      sea salt (0.25 wt ratio from Seinfeld and Pandis, 1998 table 7.8)
c      since presumably the sodium comes from sea-salt.
c      This may not hold for soil or rock derived aerosol.
c
cbk          call step(nsec,q)   ! tmg (10/22/02) bkoo (03/07/03)
c     endif ! do every time
c
      do i=1,ntotal
        save1(i)=q(i)
      enddo

      do  i=1,nsp
        qmass(i)=0.
        do k=1,nsec
          qmass(i)=qmass(i)+q((k-1)*nsp+i)
        enddo
      enddo
c
      aerm = chaero ! tmg 01/11/02

 100  ifail = 0 ! bkoo, 01/08/01

      if(aerm.eq.'MADM') then
         ntotalx  = ntotal
         nsecx    = nsec
         ntotalx2 = 0 ! not used
         nsecx2   = 0
         call madm(t0,t1,q)
      elseif(aerm.eq.'EQUI') then
         ntotalx  = 0 ! not used
         nsecx    = 0 ! not used
         ntotalx2 = ntotal
         nsecx2   = nsec
         call eqpart(t1,q)
      elseif(aerm.eq.'HYBR') then
         ntotalx  = ntotald
         nsecx    = nsecd
         ntotalx2 = ntotale
         nsecx2   = neqsec
         call heqdyn(t0,t1,q)
      else
         write(6,*)'ERROR: AEROSOL MODULE should be MADM, EQUI, or HYBR.'
         stop
      endif

c     if HYBR failed try EQUI after restore q, dsec, qn & qtot0
      if(ifail.eq.1) then
        write(6,*)'Now try EQUI for this step...'
        do j=1,ntotal
          q(j)=save1(j)
        enddo
        do j=1,nsec
          dsec(j)=save2(j)
          qn(j)=save3(j)
        enddo
        do j=1,nexti
          qtot0(j)=save4(j)
        enddo
        aerm='EQUI'
        goto 100
      endif
CDEBUG
c	print *,'AERCHEM: Before Interpolation. con(NO2)='

c
c     Linear Interpolation Routine
      do i=1,nexti
        qtot1(i)=0.0d0
      enddo
      prod=rgas*temp/(1.01325d5*pres) ! pres (bkoo, 6/9/00)
      qtot1(kH2O)= 0.d0               ! gas phase mass
      qtot1(kSO4)= q(naer+ih2so4)/(prod/gmw(ih2so4))  
      qtot1(kNa) = 0.d0  
      qtot1(kNO3)= q(naer+ihno3) /(prod/gmw(ihno3))  
      qtot1(kNH4)= q(naer+inh3)  /(prod/gmw(inh3))  
      qtot1(kCl) = q(naer+ihcl)  /(prod/gmw(ihcl))  
      do ki=1,nexti                 
        do i=1,nsec
          qtot1(ki)=qtot1(ki)+q((i-1)*nsp+ki) ! aerosol mass
        enddo
      enddo
c
c     get dry diameter before calling newdist tmg (01/25/02))
c     (only needed for MADM)
      if (aerm.eq.'MADM') call ddiameter(q)

c
cbk      if (.not.lsoap) call newdist(tt1,q) 
      call newdist(t1,q) ! bkoo (03/09/03)
c
      do  i=1,nsp
        qmass(i)=0.
        do k=1,nsec
          qmass(i)=qmass(i)+q((k-1)*nsp+i)
        enddo
      enddo

      do i=1,nexti
        qtot2(i)=0.0d0
      enddo
      prod=rgas*temp/(1.01325d5*pres) ! pres (bkoo, 6/9/00)
      qtot2(kH2O)= 0.d0	              ! gas phase mass
      qtot2(kSO4)= q(naer+ih2so4)/(prod/gmw(ih2so4))  
      qtot2(kNa) = 0.d0  
      qtot2(kNO3)= q(naer+ihno3) /(prod/gmw(ihno3))  
      qtot2(kNH4)= q(naer+inh3)  /(prod/gmw(inh3))  
      qtot2(kCl) = q(naer+ihcl)  /(prod/gmw(ihcl))  
      do ki=1,nexti                 
        do i=1,nsec
          qtot2(ki)=qtot2(ki)+q((i-1)*nsp+ki) ! aerosol mass
        enddo
      enddo

      call step(nsec,q) ! tmg (01/31/02)
c
c      ----- End Main loop -----
      return
      end
