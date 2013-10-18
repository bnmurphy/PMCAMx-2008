c****************************************************************************
c THIS IS THE ROUTINE THAT CALCULATES THE CHANGES TO THE AEROSOL
c MASS DISTRIBUTION BECAUSE OF COAGULATION FOR Dt.
c
c A HIGH RESOLUTION DISTRIBUTION IS USED FOR THE COAGULATION CALCULATIONS
c BY DIVIDING EACH SECTION OF THE ORIGINAL DISTRIBUTION TO TEN SECTIONS
c THE NEW DISTRIBUTION IS TREATED AS A POPULATION OF MONODISPERSE
c AEROSOL GROUPS. FOR 10 SECTIONS, THE RESULTING CALCULATION USES
c 100 AEROSOL GROUP SOMETHING THAT HAS PROVEN TO BE ACCURATE FOR
c SUCH A CALCULATION
c
c****************************************************************************
c
c PARAMETERS
c   Q(I) : VECTOR OF AEROSOL AND GAS-PHASE CONCENTRATIONS FOR THE
c          COMPUTATIONAL CELL (UNITS: ug/m3 for aerosol)
c   DIAMETROS(i): MEAN DIAMETER OF HIGH RESOLUTION SECTION (um)
c   KCG(I,J): COAGULATION COEFFICIENT BETWEEN PARTICLES OF SECTIONS I AND J
c             THE UNITS ARE cm3/hr
c   CGNEWP(I,J): THE SECTION TO WHICH THE AGGREGATE PARTICLE RESULTING
c                FROM THE COLLISIONS OF PARTICLES FROM SECTION I AND J
c                BELONGS
c   PN(I): PARTICLE NUMBER CONCENTRATION HIGH RESOLUTION (UNITS: cm-3)
c   M(I) : MASS OF ONE PARTICLE IN HIGH RESOLUTION SECTION I (UNITS: ug)
c   DENS : PARTICLE DENSITY (UNITS: ug/m3)
c   DT   : OPERATOR TIMESTEP IN SECONDS
c



      subroutine coagul(q)
c
      include 'dynamic.inc'
      parameter (ndcoag = nsec*nres)
      parameter (dens = 1.0d12)          ! in ug/m3
c
      real*8 q(ntotal), pmtotal(nsec), pm(ndcoag)
      real*8 diametros(ndcoag), pn(ndcoag), m(ndcoag)
      real*8 kcg(ndcoag, ndcoag)
      real*8 dmdt(ndcoag, nsp), dqdt(ntotal)
      real*8 basi, sint, fraction1, fraction2, change
      integer index, index1, index2, isec1, isec2
      integer isec, jres, i ,j, k
      integer cgnewp(ndcoag,ndcoag)
c
c     RETURN IF NO COAGULATION IS REQUIRED
c
      if (icoag .eq. 0) return
c
c CALCULATION OF THE DIAMETERS OF THE HIGH RESOLUTION DISTRIBUTION
c
      basi = dlog10(dsecf(1))
      sint = (dlog10(dsecf(2)) - dlog10(dsecf(1)))/FLOAT(nres)
      basi = basi + sint/2.0d0            ! mid-point of section
      do isec=1, nsec
       do jres=1, nres
        diametros((isec-1)*nres+jres) = 10.0d0**basi
        basi = basi + sint
       enddo
      enddo
c
c CALCULATION OF THE COAGULATION COEFFICIENT AND COLLISION ARRAYS
c
      call kcoag(diametros, kcg, cgnewp)
c
c CALCULATION OF THE NUMBER DISTRIBUTION
c
c* 1. Zeroing of the total mass of each section
c
      do isec = 1, nsec
       pmtotal(isec) = 0.0d0
      enddo
c
c* 2. Calculation of the total mass of each section
c
      do isec = 1, nsec
       do j =1 , nsp
        index = (isec-1)*nsp + j
        pmtotal(isec) = pmtotal(isec) + q(index)
       enddo
      enddo
c
c* 3. Calculation of total mass of each high resolution section
c
      do isec = 1, nsec
       do j = 1, nres
        index = (isec-1)*nres+j
        pm(index) = pmtotal(isec)/FLOAT(nres)
       enddo
      enddo
c
c* 4. Calculation of the number and single particle mass of each high resolution section
c
      do i=1, ndcoag
       m(i) = pi*dens*diametros(i)**3/6.0d18         ! in ug
       pn(i)= 1.0d-6*pm(i)/m(i)                       ! in part/cm3
      enddo
c
c ZEROING OF THE COAGULATION RATE VECTORS
c
      do i=1, ndcoag
       do k=1, nsp
        dmdt(i, k) = 0.0d0
       enddo
      enddo 
c
c CALCULATION OF THE COAGULATION RATES
c
      do i=1, ndcoag
       do j=1, i
         cgrate = kcg(i,j)*pn(i)*pn(j)        ! UNITS cm-3 hr-1
         if ( i .eq. j) cgrate=cgrate/2.0d0   ! to avoid double-counting of i,i collisions
c
       do k=1, nsp
       index1 = ((i-1)/nres)*nsp+ k              ! species position in Q vector
       index2 = ((j-1)/nres)*nsp+ k
       isec1 = ((i-1)/nres) + 1                  ! section in low resolution
       isec2 = ((j-1)/nres) + 1
c
       if (pmtotal(isec1) .le. 0.0d0) then
        fraction1 = 0.0d0
        else
        fraction1 = q(index1) / pmtotal(isec1)    ! fraction of section mass
       endif
c
       if (pmtotal(isec2) .le. 0.0d0) then
        fraction2=0.0d0
        else
        fraction2 = q(index2) / pmtotal(isec2)
       endif
c
       dmdt(i,k)=dmdt(i,k)- cgrate*m(i)*1.d6*fraction1   ! UNITS ug/m3 hr
       dmdt(j,k)=dmdt(j,k)- cgrate*m(j)*1.d6*fraction2   ! UNITS ug/m3 hr
       dmdt(cgnewp(i,j), k)=dmdt(cgnewp(i,j), k) +
     &   cgrate*(m(i)*fraction1+m(j)*fraction2)*1.d6
       enddo
       enddo
      enddo
c
c MASS BALANCE CHECK 1: ** TO BE REMOVED **
c
csp       change = 0.0
csp       do i = 1, ndcoag
csp       do j=1, nsp
csp       change = change + dmdt(i,j)
csp       enddo
csp       enddo
csp       write(3,*)tcom, change
c
c     CALCULATION OF THE COAGULATION RATES FOR THE DIFFERENT AEROSOL SPECIES
c     TRANSITION FROM THE HIGH RESOLUTION TO THE LOW RESOLUTION DISTRIBUTION
c
       do i=1, ntotal
       dqdt(i)=0.0d0
       enddo
c
       do i = 1, ndcoag
       do k = 1, nsp
       index = ((i-1)/nres)*nsp + k
       dqdt(index) = dqdt(index) + dmdt(i,k)        ! in ug m3 hr
       enddo
       enddo
c
c  MASS BALANCE CHECK 2: ** TO BE REMOVED **
c
csp       change = 0.0
csp       do i=1, ntotal
csp       change = change + dqdt(i)
csp       enddo
csp       write(4,*)tcom, change

c
c UPDATING OF THE AEROSOL SIZE/COMPOSITION DISTRIBUTION AFTER COAGULATION
c
c NOTE: IT IS POSSIBLE THAT A GIVEN SECTION MAY DISSAPPEAR COMPLETELY
c       BECAUSE OF OUR ASSUMPTION OF A CONSTANT COAGULATION RATE FOR Dt.
c       THIS CAN RESULT IN A SMALL GAIN OF MASS AND NEGATIVE CONCENTRATIONS
c
c       TO CORRECT THIS THE CONCENTRATION OF THE SPECIES IS SET TO ZERO
c       AND TO CONSERVE MASS THE CONCENTRATION OF THE SAME SPECIES
c       IN THE LARGER SECTION IS REDUCED BY THE SAME AMOUNT
c
       do i = 1, ntotal
       q(i) = q(i) + dqdt(i)*(dt/3600.0d0)
c
        if (q(i) .lt. 0.0d0) then
csp        write(6,*) 'COAGULATION RESULTED IN NEGATIVE CONCENTRATION'
csp        write(6,*) 'MASS GAIN =', -q(i), 'ug/m3'
cgy the syntax stop 'string' doesn't work well on all platforms
cgy        if (i+nsp.gt.naer) stop'COAGUL: beyond aerosol section' ! bkoo (6/8/00)
        if (i+nsp.gt.naer) then
          write(6,*)'COAGUL: beyond aerosol section'
          stop
        endif
	q(i+nsp)=q(i+nsp)+q(i)
        q(i) =0.0d0
        endif
       enddo
c
c END OF CALCULATION
c
      return
      end




      subroutine kcoag(dp, kcg,cgnewp)
c
c     This routine calculates the coagulation coefficient (K) for Brownian
c     diffusion using the formulas of Fuchs (Seinfeld, 1986, Table 10.1) 
c     and puts them in the KCG array
c     It also determines in what section the new particle will exist and 
c     puts the section number in the CGNEWP array
c
c
      include 'dynamic.inc'
c
      integer ndcoag
      real*8 rho, airl, visc
c
      parameter (ndcoag = nsec*nres)
      parameter (rho = 1.0d0)      ! density of particles in g/cm3
      parameter (airl = 0.0651d0)  ! air mean free path in um
      parameter (visc= 1.83d-4)    ! air viscosity in poise
c
      real*8 kcg(ndcoag,ndcoag),difc(ndcoag),cc(ndcoag),cbar
      real*8 gbar,g(ndcoag),l(ndcoag),mass(ndcoag),knudsen(ndcoag)
      real*8 difsum,dpsum,vol(ndcoag),newvol,volbrk
      real*8 dp(ndcoag)
c
      integer cgnewp(ndcoag,ndcoag)
      integer i,j,k
c     
      do i=1,ndcoag      
       vol(i)=dp(i)**3*pi/6.d0                       ! in um3
       mass(i)=vol(i)*rho*1.0d-12                   ! in grams
       cc(i)=DSQRT(8.d0*temp*1.381d-20/pi/mass(i))   ! in m/s
       knudsen(i)=2.0d0*airl/dp(i)
       difc(i)=1.381d-20*temp/3.d0/pi/visc/dp(i)*1.0d4*
     & ((5.d0+4.d0*knudsen(i)+6.d0*knudsen(i)**2+18.d0*knudsen(i)**3)/
     & (5.0d0-knudsen(i)+(8.0d0+pi)*knudsen(i)**2))    !in m2/s
       l(i)=8.0d0*difc(i)/pi/cc(i)*1.d6               !in um
       g(i)=1.0d0/(3.0d0*dp(i)*l(i))*((dp(i)+l(i))**3-
     &      DSQRT((dp(i)**2+l(i)**2)**3))-dp(i)        !in um
      enddo
c
       do i=1,ndcoag
        do j=1,ndcoag
         cbar=DSQRT(cc(i)**2+cc(j)**2)             ! in m/s
         gbar=DSQRT(g(i)**2+g(j)**2)             ! in um
         difsum = difc(i)+difc(j)            ! in m2/s
         dpsum = dp(i)+dp(j)                     ! in um
         kcg(i,j)=2.0d0*pi*difsum*dpsum*3600.0d0/
     &   (dpsum/(dpsum+2.0d0*gbar)+8.0d0*difsum/cbar/dpsum*1.d6) ! in cm3/hr
c
         newvol=vol(i)+vol(j)
         k=max(i,j)
         cgnewp(i,j)=0
              do while (cgnewp(i,j) .eq. 0)
                if (k .lt. ndcoag) then
                   volbrk = 0.5d0*(vol(k)+vol(k+1))
		 else
                   volbrk=2.0d0*vol(ndcoag)+1.0d0
		endif
c
                if(newvol .lt. volbrk) then
                   cgnewp(i,j)=k
		 else
		   k = k+1
		endif
	      enddo
c
       enddo
      enddo
c
      return
      end
