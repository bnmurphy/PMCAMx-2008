c******************************************************************
c                        AEROSOL PARAMETERS 
c******************************************************************
c 
c USEFUL CONSTANTS
c
      include 'section_aq.inc'
c
      parameter (pi1 = 3.14159)
      parameter (pi6 = pi1/6.0)
c
c number and distribution of sections
c   section distribution:
c
c      |       |       |       | ...  |           |          |
c        dp(1)   dp(2)   dp(3)   ...   dp(nsect-1)  dp(nsect)
c       =dpmin                                       =dpmax
c
c dpmin and dpmax are the geometric average diameter of particles in the
c   smallest and largest bins.  Cut points are given by sqrt(dp(i)*dp(i+1)).
c
c USE 10 NM AS LOWER BOUND TO AVOID EVAPORATION PROBLEMS IN AQUEOUS MODULE      
c
c      parameter (nsect = 10)
      parameter (dpmax = 10.0e-6)    ! 10.0 um max. diameter [m] 
      parameter (dpmin = 0.05e-6)     ! 0.1 um min. diam. [m]  NEW VALUE 
      parameter (rho = 1.4e12)       ! particle density [ug/m^3]
c
c AEROSOL COMPONENTS IN THE AEROSOL CONCENTRATION VECTOR
c
      parameter (nas =  1)           ! sodium
      parameter (nah =  2)           ! hydrogen
      parameter (naa =  3)           ! ammonium
      parameter (nan =  4)           ! nitrate
      parameter (nac =  5)           ! chloride
      parameter (na4 =  6)           ! sulfate
      parameter (naw =  7)           ! water
      parameter (nae =  8)           ! elemental carbon
      parameter (nao =  9)           ! organics
      parameter (nar = 10)           ! crustal
      parameter (nahso5 = 11)        ! hso5-
      parameter (nahmsa = 12)        ! hmsa
      parameter (naspec = 12)        ! number of aerosol species
c
c condensible gas-phase components in local arrays
c
      parameter (ngca =  1)          ! ammonia
      parameter (ngcn =  2)          ! nitric acid
      parameter (ngcc =  3)          ! hydrochloric acid
      parameter (ngc4 =  4)          ! gas-phase sulfate
      parameter (ngco =  5)          ! gas-phase organics
      parameter (ngcspec = 44)        ! number of condensible gas-phase species
c
c condensible gas-phase components in global gas-phase array
c
c this must be customized to have the correct addresses
c
      parameter (nga =  1)           ! ammonia
      parameter (ngn =  2)           ! nitric acid
      parameter (ngc =  3)           ! hydrochloric acid
      parameter (ng4 =  4)           ! gas-phase sulfate
      parameter (ngo =  5)           ! gas-phase organics
      parameter (ngspec = 44)         ! number of condensible gas-phase species
c
c total number of gas phase species so we know where the aerosol starts
c
      parameter (ngtotal = 100)
      parameter (ngas_aq=ngtotal)
      parameter (naers=naspec)






