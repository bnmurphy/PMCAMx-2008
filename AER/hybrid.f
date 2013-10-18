c==========================================================================
c   12/17/02: tmg
c             - switched order so equilibrium is first again
c   10/22/02: tmg
c             - added parameter to subroutine step because not all sections
c               are always used
c   01/25/02: tmg
c             - use dry and wet basis for section diameter as appropriate
c   01/27/01: bkoo
c             - change the order of equilibrium and dynamic calculations
c             - change the iteration number
c   01/08/01: bkoo
c             - HYBR-EQUI shift
c   10/11/00: bkoo
c             - establish electroneutrality in the first dynamic section
c   09/20/00: bkoo
c             - combine eqparth & eqpart
c             - change the order of step & diameter at the end of heqdyn
c   06/08/00: bkoo
c             - add coagulation & nucleation
c   MAY 2000: bkoo
c             - correct negative gas conc. after EQPARTh
c             - EQPARTh : EQPART for HYBRID
c
cgy supress error messages to units 39 and 90 (6/12/00)
c
c=========================================================================+

c--------------------------------------------------------------------------
c   HEQDYN  - HYBRID EQUILIBRIUM AND DYNAMIC GAS/AERSOL PARTITIONING
C   ROUTINE FOR USE OF BOTH EQUILIBRIUM AND DYNAMIC GAS/AERSOL PARTITIONING
C
C   THE FIRST 'NEQSEC' SECTIONS USE EQUILIBRIUM PARTITIONING AND THE 
C   REMAINING AEROSOL SECTIONS USE DYNAMIC PARTITIONING
C
C   CONSIDERATIONS FOR WHERE TO CUT THE SECTIONS
C   1.  THE EQUILIBRIUM SECTIONS SHOULD BE SMALL ENOUGH TO REACH EQUILIBRIUM
C   IN A SINGLE TIME STEP.  DP OF ABOUT 1 UM IS PROBABLY SAFE FOR A 10 MIN
C   TIME STEP (SEE WEXLER AND SEINFELD, AE-V24, 1990).
C   2.  ONE SHOULD ALSO CONSIDER THE CHEMICAL HOMOGENEITY OF THE SECTIONS TO 
C   BE PLACED IN EQUILIBRIUM.  IF THE INITIAL CHEMICAL COMPOSITIONS OF 
C   THE EQUILIBRIUM SIZE SECTIONS ARE TOO DIFFERENT ERRORS WILL BE INTRODUCED
C
      SUBROUTINE HEQDYN(T0,T1,Q)
      include 'dynamic.inc'
c     master conc. array: aerosols (ug/m3) then gas (ppm)
      real*8 q(ntotal)
c     conc. array for equilibrium: aerosols (ug/m3) then gas (ppm)
      real*8 qe(ntotale)
c     conc. array for dynamics: aerosols (ug/m3) then gas (ppm)
      real*8 qd(ntotald)

      common / counters / icount, ifail ! bkoo (01/08/01)

c     gas phase "source" due to equilibrium with small sections (ppm/s)
c      real*8  sdgas(ngas),sgasrate(ngas),residual(ngas)
c      common / eqgas / sdgas,sgasrate,residual,t0c
c     Equation for transition regime (Dahneke 1983)
      f(yKnD, alfa)=(1.0d0+yKnD)/(1.0d0+2.0d0*yKnD*(1.0d0+yKnD)/alfa)
c
c     STEP 1/3: CALCULATE NUCLEATION RATE FOR THE WHOLE STEP
c
      call nucl(q)
c
c     STEP 2/3: CALCULATE COAGULATION RATE FOR THE WHOLE STEP
c
      call coagul(q)

c      call step(nsec,q) ! tmg (10/22/02)
c
c      t0c=t0
c     set useful indices
      ng =nsec*nsp
      nge=neqsec*nsp
      ngd=(nsec-neqsec)*nsp
c
c     We iterate to minimize bias towards large sections caused by emmissions
c     only occuring in the dynamic portion of the code
c     beware however, that too many iterations will bias towards the small
c     sections because instantaneous mass transfer will be a bad assumption
c
      iterations=2
      tstart=t0
      tend=tstart+(t1-t0)/float(iterations)
      do iterate=1,iterations
c
c     transfer concentrations to equilibrium array
       do isec=1,neqsec
         indx=(isec-1)*nsp
         do isp=1,nsp
           qe(indx+isp)=q(indx+isp)
         enddo
       enddo
       do isp=1,ngas
         qe(nge+isp)=q(ng+isp) ! bkoo (01/27/01)
       enddo
c     call equilibrium partitioning code
       call eqpart(tstart,qe)
c
c     transfer concentrations to dynamics array
       do isec=neqsec+1,nsec
         indx=(isec-1)*nsp
         indxd=(isec-neqsec-1)*nsp
         do isp=1,nsp
           qd(indxd+isp)=q(indx+isp)
         enddo
       enddo
c     this step dictates sequential changes to the gas phase -
c     first equilibrium then dynamic sections changes
c     Now, first dynamic then equilibrium: bkoo (01/27/01)
c     Switched back to equilibrium first: tmg (12/17/02)
       do isp=1,ngas
         qd(ngd+isp)=qe(nge+isp) ! bkoo (01/27/01) 
       enddo

c     call dynamics partitioning code
       call madm(tstart,tend,qd)
c     if HYBR failed return (bkoo, 01/08/01)
       if(ifail.eq.1) return 
c     set the gas "source" term parameters to represent the change in the gas
c     phase concentrations due to exchange with the equilibrium sections.
c     We assume that the change in the gas phase conc. decays exponentially
c     with a rate equal to the sum of the sectional rates - then NH4 is set
c     equal the rate of negative ions exchanging with the gas phase.  The
c     residual delta Cgas after a time step of exponential decay will be
c     linearized and added to the exponential "source" to preserve mass balance.
c     The decay rate given here is a maximum rate.  If a sepecies completely
c     evaporates in a section within the time step, the decay rate will be
c     slower.  This is neglected here for simplicity.
c
c       do k=1,ngas
c        sgasrate(k)=0.0d0
c        sdgas(k)=qe(nge+k)-q(ng+k)
c        do i=1,neqsec
c         indx=(i-1)*nsp
c         if(sdgas(k).gt.0.0d0.and.                  ! evaporation but
c     &	   ((k.eq.ihcl.and.q(indx+kcl).lt.tinys).or. ! species can't evap.
c     &	   (k.eq.ihno3.and.q(indx+kno3).lt.tinys))) goto 100	   
c         dp = 1.0d-6 * dsec(i)                      ! (moving sections)
c         KnD=2.0d0*lamda(k)/dp
c         sgasrate(k)=sgasrate(k)+2.0d0*pi*dp*diffus(k)*
c     &         f(KnD, delta(k))*qn(i)*1.0d6         ! in sec-1
c 100    enddo
c       enddo
c       sgasrate(inh3)=2*sgasrate(ih2so4)+sgasrate(ihno3)+sgasrate(ihcl)
c       write(50,*)(sgasrate(k),k=1,ngas)
c       do k=1,ngas
c        if(sgasrate(k)*(t1-t0).gt.20) then
c         residual(k)=0.0d0                          ! exp(-20)=~1e-9
c        else
c         residual(k)=sdgas(k)*exp(-sgasrate(k)*(t1-t0))/(t1-t0)
c        endif
c       enddo
c
c      EQPARTh can give negative gas concentrations  ! done in negchk tmg ! (01/18/02)
c       do isp=1,ngas
c         if(qe(nge+isp).lt.0.0d0) then
cgy	   write(90,*)'NEG_GAS: gas(',isp,')=',qe(nge+isp),' t=',tstart
c          qe(nge+isp)=0.0d0
c         endif
c       enddo
c
c     return new concentrations to the original array
       do isec=1,neqsec                             ! from equilibrium sections
         indx=(isec-1)*nsp
         do isp=1,nsp
           q(indx+isp)=qe(indx+isp)
         enddo
       enddo
       do isec=neqsec+1,nsec                        ! from dynamic sections
         indx=(isec-1)*nsp
         indxd=(isec-neqsec-1)*nsp
         do isp=1,nsp
           q(indx+isp)=qd(indxd+isp)
         enddo
       enddo
       do isp=1,ngas
         q(ng+isp)=qd(ngd+isp) ! bkoo (01/27/01)
       enddo

       call negchk(tend,q,nsec)

c     calculate water
c       call step(nsec,q)  ! tmg (10/22/02)

c     calculate new diameters
       call ddiameter(q)
       tstart=tend
       tend=tend+(t1-t0)/float(iterations)
      enddo

      return
      end
