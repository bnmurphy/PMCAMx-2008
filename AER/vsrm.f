      subroutine vsrm(gas, aerosol, total, tstart, tend, deltat,
     & temp, iaq, p, rh, chtype)   !kchm, height, chtype) 
ckf     & temp, iaq, p, rh, kchm, height)
ckf     & temp, iaq, p, rh)
     
c  Regarding vsrm argument list:
c  gas, aerosol, total are inputs to the subroutine decisions;
c  iend, istart, and deltat are used to calc. nsteps and tinit (an
c  input to aqoperator);
c  rh, p, iaq, and temp are also inputs to aqoperator.    
      
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'dropcom.inc'
      real*4 gas(ngas_aq), aerosol(nsect,naers)    
      real*4 lwc(2), temp, total, p, rh
      integer chtype

      REAL HU
      real*4 s1
      
      INTEGER NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
    
      common/istep/ istep                      ! used for debugging
      COMMON /SVOD02/ HU, NCFN, NETF, NFE, NJE, NLU, NNI, NQU, NST
      common/time2/tfin

c DETERMINE WHETHER TO EXECUTE THE BULK OR SR CALCULATIONS

CKF      call decisions(aerosol,gas,total,chtype)
      call decisions(aerosol,gas,total,chtype,p,temp)
      
      if (chtype .eq. 2) then 
      lwc(1) = frac1*total                  ! LWC in g/m3 section 1
      lwc(2) = frac2*total                  ! LWC in g/m3 section 2      
      else
      lwc(1) = total
      endif

c     PREPARATION FOR THE RUN

c     SIMULATE CLOUD/FOG

      nsteps=nint((tend-tstart)/deltat)
      do istep = 1, nsteps
c      
      tinit=float(istep-1)*deltat
      tfin=float(istep)*deltat

c     IF CHTYPE = 1, BULK.  IF CHTYPE = 2, 2-SECTION SIZE-RESOLVED
      
      if (chtype .eq. 1) then
        call aqoperator1(tinit, deltat, gas, aerosol,
     & lwc, rh, temp, p, iaq)
        else 
       call aqoperator2(tinit, deltat, gas, aerosol,
     & lwc, rh, temp, p, iaq)
        endif
     
       enddo
 
     
      s1 = 0.
      do i=1, nsect
        s1=s1+aerosol(i,ksvi)*32./96.      ! total sulfur in sulfate (ug/m3)
      enddo
      return
      end
