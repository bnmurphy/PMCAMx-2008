c=========================================================================+
c     03/09/03: bkoo - added call to orgps1 to calculate csatT
c     June 8 2000: bkoo                                                   +
c                  - add coagulation & nucleation in madm                 +
c                                                                         +
c     MAY, 2000: modified by bkoo                                         +
c                - combine MADM with HYBRID                               +
c                                                                         +
cgy supress writes to unit 6 (6/12/00)
c
c=========================================================================+

C     Last change:  KC   17 Feb 98    5:18 pm
c	dummy jacobian subroutine (needed by lsode)
c *************************************************************************
c               WRITTEN BY DR. CHRISTODOULOS PILINIS
c                          April 1997
c *************************************************************************
      subroutine madm(t0,t1,q)
c
c     This is the main routine of the new generation
c     multicomponent aerosol dynamics model (MADM)
c *************************************************************************
c               WRITTEN BY DR. CHRISTODOULOS PILINIS
c                          April-May 1997
c *************************************************************************
c
c     t0=initial time (seconds)
c     t1=final time (seconds)
c
      include 'dynamic.inc'

      real*8 q(ntotal)
      dimension pgas(ngas)
c     concetration derivatives of aerosol species in ugr/m3/sec and 
c     gases in ppm/sec.
      common / counters / icount, ifail ! bkoo (01/08/01)

      if(aerm.eq.'MADM') then
c
c     STEP 1/3: CALCULATE NUCLEATION RATE FOR THE WHOLE STEP
c
        call nucl(q)
c
c     STEP 2/3: CALCULATE COAGULATION RATE FOR THE WHOLE STEP
c
        call coagul(q)
c
      endif

c
      call CALCDIF              ! calculate the diffusion coefficients

cbk   calculate the saturation concentrations for organics - bkoo (04/18/02) (03/09/03)
      call orgps1(q)    

      maxstep=200000            ! bkoo (01/08/01)

      neq=ntotalx               ! number of equations
      ng=neq-ngas
      t=t0                      ! initialize time for lsode
      tout=t1                   ! initialize time for lsode
c
c     dt=t1-t0
c
c     ********************* FIRST STEP *********************
c
c     find the change in the size composition distribution
c     due to primary emissions and transport for all sections
c
ck    do i=1,naer
ck       q(i)=q(i)+dt*qsource(i) ! ugr/m3/sec		
ck    end do
c
c     qtot0(kso4)=qtot0(kso4)+4.05e-4*dt   ! add source term to mass balance
c
c     ********************* SECOND STEP *********************
c
c     Call the equilibrium codes to calclulate partial pressures
c     of volatile organic and inorganic species in various sections
c     as well as the physical state of the aerosol (dry/liquid)
c
c     call step(q) ! done in DIFFUND - bkoo (MAY, 2000)
c      call step(q) ! If step in DIFFUND is located at the end of DIFFUND
c                    then this step should be called here.
c                    (bkoo, MAY 22, 2000)
c-----------------------------------------------------------------------------
c
c     ********************* THIRD STEP *********************
c
c     Calculate the mean speed of condensible gas molecules 
c     and the mean free path of condensible gas molecules
c     ( for Kn calculations).
c
      do kx=1,ngas
         vel(kx)=sqrt(8.0d0*rgas*temp/pi/(gmw(kx)*1.0d-3))     ! in m/sec
c
c     calculate lamda when Dahneke equation is used ( TABLE 8.4 )
c
         lamda(kx)= 2.0d0*diffus(kx)/vel(kx)
      end do
c
c     Calculate the evolution of the size composition distribution
c     due to the condensation of all condensible species
c
c
      do while(t.lt.tout)
         tlast=t
         isfirst=.TRUE.   ! bkoo (03/05/02)

         call diffund(neq,t,q,tout) ! tmg (07/02/01) 

         do i=1,nsecx
          ims(i)=0        ! restore IMS's - bkoo (03/05/02)
         enddo

c
c     if failed return with ifail=1 (bkoo, 01/08/01)
c         icount=icount+iwork(11) ! (tmg,  01/21/02)
         if(icount.ge.maxstep .and. t.lt.tout-1.d-5) goto 999

c     put into diffund.f instead ! tmg  (07/02/01)
c         if(aerm.eq.'HYBR') then ! bkoo (02/09/01)
c           call negchk(t,q,nsecd)
c         else
c           call negchk(t,q,nsec)
c         endif

cgy	 write(6,202)t,t-tlast,temp,rh,icount
      enddo
 202  format(2g15.8,2g10.4,' number of steps=',i8)

      return
      
 999  write(*,*)"DLSODE failed here"
      ifail = 1
      return      
      
      end
