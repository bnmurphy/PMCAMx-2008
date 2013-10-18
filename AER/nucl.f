c****************************************************************************
c THIS IS THE SUBROUTINE FOR H2SO4/H2O NUCLEATION CALCULATIONS.
c
c IT ASSUMES A LINEAR SULFURIC ACID VAPOR CONCENTRATION
c VARIATION FOR THE PERIOD OF THE CALCULATIONS BASED ON
c THE INITIALLY AVAILABLE SULFURIC ACID (IT WILL BE ZERO
c AT THE END OF THE STEP) AND ALSO ASSIGNS ALL THE NUCLEATED
c MASS TO THE FIRST SECTION OF THE DISTRIBUTION
c
c THIS IS APPROPRIATE FOR AN ESTIMATE OF THE EFFECT OF NUCLEATION
c ON THE MASS DISTRIBUTION BUT NOT ON THE NUMBER DISTRIBUTION UNLESS
c THE LOWER LIMIT OF THE DISTRIBUTION HAS BEEN CHOSEN TO BE A FEW NM
c
c THE NUCLEATION RATE USED IS BASED ON THE PARAMETRIZATION OF
c RUSSELL/PANDIS/SEINFELD (JGR, 99, P. 20989-21003, 1994) OF
c THE WORK OF JAECKER/VOIROL AND MIRABEL (ATM. ENVIRON, 23, 2053-2957, 1989)
c
c A SULFURIC ACID NUCLEUS SIZE OF Dnucl IS ASSUMED
*****************************************************************************
c

      subroutine nucl(q)
      include 'dynamic.inc'
      double precision q(ntotal), sulf
      double precision sulf0,tuner, ratelog, rate, jnucl, dtime
      double precision dnucl, dens, dmass, rnucl, conv, totrate
      integer nstep
c
c     RETURN IF NO NUCLEATION IS REQUIRED (bkoo, 6/8/00)
c
      if (inucl .eq. 0) return
c
c If the sulfuric acid concentration is zero, the nucleation rate is zero.
      if(q(naer+ih2so4).le.0.0d0) return	! bkoo, 6/8/00
c
c EMPIRICAL NUCLEATION SCALING FACTOR   
c  (A number of investigators have suggested that the nucleation
c   rate in the atmosphere is a lot higher than what the classical
c   nucleation theory predicts. They scale therefore the nucleation
c   rate by a factor of "tuner". Values ranging from 1 to 10^7
c   have been used. See Russell et al., 1994 for a discussion)
c
      tuner = 1.d0
c
c CALCULATION OF THE H2SO4(g) IN MOLECULES/CM3
c
      sulf0 = q(naer+ih2so4)                 ! in ppm
      sulf0 = sulf0*1.d-9*pres/(0.08206*temp) ! in moles/cm3
      sulf0 = sulf0*6.023e23                  ! in molecules/cm3
c
c CALCULATION OF THE NUCLEATION RATE (in particles cm-3 s-1)
c
c
      nstep = 100
      totrate = 0.0d0
      dtime = dt/FLOAT(nstep)
c
c      do istep = 0, nstep
      do istep = 0, nstep-1			! bkoo, 6/8/00
      sulf = sulf0 - sulf0*FLOAT(istep)/FLOAT(nstep)
      ratelog = DLOG10(tuner) - (64.24+4.7*rh) +
     &   (6.13+1.95*rh)*DLOG10(sulf)
      rate = (10.d0**ratelog)*dtime          ! in particles cm-3
      totrate = totrate + rate               ! in particles cm-3
      enddo
c
      jnucl = totrate*1.d6              ! in particles m-3 
c
c CALCULATION OF THE MASS NUCLEATION RATE (in ug m-3 s-1)
c IT IS ASSUMED THAT THE EQUIVALENT DIAMETER OF SULFURIC ACID
c IS 3 nm. THE ACTUAL DIAMETER OF THE PARTICLE WILL BE A FEW
c nm HIGHER BECAUSE OF THE ADDITION OF WATER
c
      dnucl = 3.0d-9                                 ! nuclei diameter in m
      dens  = 1.0d12                                 ! nuclei density in ug/m3
      dmass = (1.0d0/6.0d0)*pi*dens*dnucl**3         ! nuclei mass in ug/particle
      rnucl = jnucl*dmass*FLOAT(inucl)               ! in ug SULF m-3 (per step)     
c
c ADJUSTMENT OF THE GAS H2SO4 CONCENTRATION FOR RNUCL
c
      conv = (rgas*temp)*rnucl/(pres*1.01325d5*98.d0)! conversion from ug SULF/m3 to ppm 
c
c CHECKING IF WE NUCLEATED MORE THAN WHAT IS AVAILABLE
c
      if (conv .gt. q(naer+ih2so4) ) then
      conv = q(naer+ih2so4)
      q(naer+ih2so4) = 0.0d0
      q(kso4) = q(kso4) + (conv*pres*1.01325d5*98.d0)/(rgas*temp)
      return
      endif
c
      q(naer+ih2so4) = q(naer+ih2so4) - conv
c
c ADJUSTMENT OF THE SULFATE MASS IN THE FIRST SECTION FOR RNUCL
c
      q(kso4) = q(kso4)+rnucl             !sulfate change accounting for MW
c
c
c FOR DEBUGGING ONLY (PLEASE REMOVE)
c H2SO4 (in ppt), rate (part per cm3 per step) , mass nucl (ug/m3)
c
c      write(77,*)q(naer+ih2so4)*1.d6,totrate, rnucl


      return
      end

