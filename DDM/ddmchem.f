      subroutine ddmchem(sddm,nfam,nsen,nrad,ngas,nr,nf,ns,nrf,nfs,nrfs,
     &                   njac,nrxn,y,prod,stmp,ipvt,djac,amx,bmx,cmx,
     &                   lzero,dt,wt,cncrad,conc,conc1,bdnl,rr,ddmjac,
     &                   knxoy,kno3,kn2o5,ko1d,ko,ierr,lnoo)
c
c-----CAMx v4.02 030709
c
c     Revised 8/04/00
c
c     DDMCHEM advances the sensitivty coefficients one time step dt.
c     nr, nf and ns are dimensions of key matrix sub-blocks.
c     At night, O atoms are not solved and nr has been reduced by 2.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Input arguments:
c
c        sddm              sensitivity matrix sddm(nfam,nsen)
c        nfam              number of sensitivity families
c        nsen              number of coefficients in each family
c        nrad              number of radicals for this mechanism
c        ngas              number of gas spec for this mechanism
c        nr                number of radicals being solved
c        nf                number of fast state species being solved
c        ns                number of slow state species being solved
c        nrf               nr+nf
c        nfs               nf+ns
c        nrfs              nr+nf+ns
c        njac              mumber of species in the full Jacobian
c        nrxn              number of reactions
c        lzero             vector indicating species to zero out
c        dt                integration time interval (hr)
c        wt                weighting factor in fast species integration
c        cncrad            radical species concentrations (ppm)
c        conc              fast and slow species concentrations (ppm)
c                          fast species =wt*conc1 + (1-wt)*conc0
c                          slow species = conc0  where
c                          conc1 is concentration at end of timestep (t1)
c                          & conc0 is concentration at start of timestep (t0)
c        conc1             fast species concentrations (ppm) at time t1
c        bdnl              lower bounds for fast and slow species (ppm)
c        rr                reaction rates (per hour)
c        ddmjac            name of subroutine to calc full Jacobian
c        knxoy             pointer to species NXOY
c        kno3              pointer to radical NO3
c        kn2o5             pointer to radical N2O5
c        ko1d              pointer to radical O(1D)
c        ko                pointer to radical O(3P)
c        ierr              error flag from LU decomposition
c        lnoo              O atoms not solved here if true
c
c        The following arguments are passed solely to allow their size 
c        to be adjustable: y,prod,stmp,ipvt,djac,amx,bmx,cmx.
c
c     Output arguments:
c
c        sddm              updated sensitivity coefficients
c        ierr              error flag from LU decomposition
c
c     Routines called:
c
c        ddmjac
c        sgefa
c        sgesl
c
      implicit none
c
      integer nfam,nsen,nrad,ngas,nr,nf,ns,nrf,nfs,nrfs,njac,nrxn,
     &        knxoy,kno3,kn2o5,ko1d,ko,mo,
     &        ipvt(nrf),
     &        i,j,k,ierr
      real conc(ngas), cncrad(nrad), conc1(nfs), bdnl(nfs), rr(nrxn), 
     &     dt, wt, fuzz, sddm(nfam,nsen), y(njac+1),
     &     prod(nrf), stmp(nrfs), djac(njac+1,njac+1),
     &     amx(nrf,nfs), bmx(nrf,nrf), cmx(ns,nrfs)
      parameter (FUZZ = 10.0)
      logical lnoo, lzero(nfs)
      external ddmjac
c
c---  Entry point
c
      ierr=0
c
c---  Load all the concentrations into Y
c     Protect against potential divide by zero in DDMJAC
c
      mo=0
      if (lnoo) then
        mo=ko-ko1d+1
        do i=1,mo
          y(i)=1.0
        enddo
      endif
      do i=1,nr
        y(i+mo)=cncrad(i+mo)
      enddo
      do i=1,nfs
        y(i+nr+mo)=conc(i)
      enddo
      y(njac+1)=1.
c
c---  Get the full Jacobian for this mechanism
c
      call ddmjac(njac,nrxn,y,djac,rr)
c
c---  Remove terms for species not in solution
c
      do j=1,nfs
        if (lzero(j)) then
          do k=1,nfam
            sddm(k,j)=0.
          enddo
          do i=1,njac
            djac(j+nr+mo,i)=0.
            djac(i,j+nr+mo)=0.
          enddo
        endif
      enddo
c
c---  Fill the matrices A, B and C block by block
c     Note that ddmjac returns the negative Jacobian
c     djac(i,j) = -dfi/dyj
c
c---  A[-Jrs]
c
      do i=1,nr
        do j=1,ns
          amx(i,j+nf)=djac(i+mo,j+nrf+mo)
        enddo
      enddo
c
c---  A[dt*Jfs]
c
      do i=1,nf
        do j=1,ns
          amx(i+nr,j+nf)=-dt*djac(i+nr+mo,j+nrf+mo)
        enddo
      enddo
c
c---  A[-Jrf*(1-wt)] and B[Jrf*wt]
c
      do i=1,nr
        do j=1,nf
          amx(i,j)=djac(i+mo,j+nr+mo)*(1.0-wt)
          bmx(i,j+nr)=-djac(i+mo,j+nr+mo)*wt
        enddo
      enddo
c
c---  A[I+dt*(1-wt)*Jff] and B[I-dt*wt*Jff]
c
      do i=1,nf
        do j=1,nf
          amx(i+nr,j)=-dt*djac(i+nr+mo,j+nr+mo)*(1.0-wt)
          bmx(i+nr,j+nr)=dt*djac(i+nr+mo,j+nr+mo)*wt
        enddo
      enddo
      do k=1,nf
        amx(k+nr,k)=amx(k+nr,k)+1.
        bmx(k+nr,k+nr)=bmx(k+nr,k+nr)+1.
      enddo
c
c---  B[Jrr]
c
      do i=1,nr
        do j=1,nr
          bmx(i,j)=-djac(i+mo,j+mo)
        enddo
      enddo
c
c---  B[-dt*Jfr]
c
      do i=1,nf
        do j=1,nr
          bmx(i+nr,j)=dt*djac(i+nr+mo,j+mo)
        enddo
      enddo
c
c---  C[dt*Jsr, dt*Jsf, dt*Jss]
c
      do i=1,ns
        do j=1,nr
          cmx(i,j)=-dt*djac(i+nrf+mo,j+mo)
        enddo
        do j=1,nfs
          cmx(i,j+nr)=-dt*djac(i+nrf+mo,j+nr+mo)
        enddo
      enddo
c
c---  Decompose matrix B for later use by SGESL
c     Check for zero determinant
c
      call sgefa(bmx,nrf,nrf,ipvt,ierr)

      if (ierr.ne.0) then
        write(*,*) ' Zero determinant in SGEFA in DDMCHM at ', ierr
        write(*,*) ' # of radicals  = ', nr
        write(*,*) ' # of fast spec = ', nf
        write(*,*)
        return
      endif
c
c---  Update the sensitivities one family at a time
c     1. Begin with (Sf0,Ss0)
c     2. Calculate the product PROD = A*(Sf0,Ss0)
c     3. Solve B*(Sr1,Sf1) = PROD using saved decomposition of B
c     4. Calculate Sfbar = Sf1*wt + Sf0*(1-wt)
c     5. Calculate Ss1 = Ss0 + C*(Sr1,Sfbar,Ss0)
c     6. Save the updated sensitivities (Sf1,Ss1)
c     7. Set S(NXOY) = S(NO3) + 2*S(N2O5)
c
      do k=1,nfam
c
        do i=1,nrf
          prod(i)=0.
        enddo
        do i=1,nrf
          do j=1,nfs
            prod(i)=prod(i)+(amx(i,j)*sddm(k,j))
          enddo
        enddo
c
        call sgesl(bmx,nrf,nrf,ipvt,prod,0)
c
        do i=1,nr
          stmp(i)=prod(i)
        enddo
c
c---  Reset sensitivities Sf1 to 0 if concentrations for fast species
c     were reset to lower bound at end of timestep in trap.
c
        do i=1,nf
          if(conc1(i) .le. FUZZ*bdnl(i)) prod(i+nr) = 0.
          stmp(i+nr)=(sddm(k,i)*(1.-wt))+(prod(i+nr)*wt)
          sddm(k,i)=prod(i+nr)
        enddo
c
        do i=1,ns
          stmp(i+nrf)=sddm(k,i+nf)
        enddo
c
        do i=1,ns
          do j=1,nrfs
            sddm(k,i+nf)=sddm(k,i+nf)+(cmx(i,j)*stmp(j))
          enddo
        enddo
c
        sddm(k,knxoy)=prod(kno3-mo) + 2.*prod(kn2o5-mo)
c
      enddo 
c
      return
c
      end 
