C=======================================================================
C
C *** MADM CODE
C *** SUBROUTINE CALCDIF
C *** THIS FUNCTION CALCULATES THE diffusUSION COEFFICIENTS OF GASES 
C     INTO AIR.
C
C ======================== ARGUMENTS / USAGE ===========================
C
C  INPUT:
C  1. [T] 
C     DOUBLE PRECISION variable, the temperature (K). 
C
C  2. [P] 
C     DOUBLE PRECISION variable, the pressure (atm). 
C
C  OUTPUT:
C  1. [diffus]
C     DOUBLE PRECISION array, of length 4, the diffusion coeffs (m2/s):
C     diffus(1) = D(NH3  -AIR)
C     diffus(2) = D(HNO3 -AIR)
C     diffus(3) = D(H2SO4-AIR)
C     diffus(4) = D(HCL  -AIR)
C 
C=======================================================================
C
      SUBROUTINE CALCDIF 
C
	include 'dynamic.inc'
C 
      DOUBLE PRECISION KTE
      COMMON /CINT/ KTE(81), OMGD(81), SIGM(5), PARM(5)
C
	t=temp					! temperature (Kelvin)
        p=pres					! atmospheric pressure (atm)
						! pres (bkoo, 6/9/00)
	amw=29.2d0				! molecular weight of air


      COEF1 = (1.8583D-3/P)*T**(1.5D0)

      DO 10 I=1,4
        SIG     = 0.5*(SIGM(I)+SIGM(5))       	! SIGMA_{12}
        AKTE    = SQRT(T*T/PARM(I)/PARM(5))   	! kT/e_{12}
        OMEGA   = COLINT(AKTE)                	! Collision integral
        COEF2   = SQRT((GMW(I)+AMW)/(GMW(I)*AMW)) 
        diffus(I) = COEF1*COEF2/(SIG*SIG*OMEGA) ! Diff.coef in cm2/s
        diffus(I) = diffus(I)*1.0D-4		! Convert to m2/s
10    CONTINUE
c
c     Diffusion coefficient of all organics, approximation based on
c     eq 11-4.1 in Reid & Prausnitz Edition 4, Page 587
c     Wilke, C. R. and C. Y. Lee: Ind. Eng. Chem., 47: 1253 (1955)
c     We could more elaborate this... ???
c     - bkoo (04/18/02)
c     - revised for PMCAMx by bkoo (03/09/03)
      do i=1,ngas-4
        diffus(ICA41-1+i) = 0.0791 * 1.d-4
      enddo 
C
      RETURN
      END


C=======================================================================
C
C *** MADM CODE
C *** DOUBLE PRECISION FUNCTION COLINT
C *** THIS FUNCTION CALCULATES THE COLLISION INTEGRAL {OMEGA_D} NEEDED 
C     FOR ESTIMATING DIFFUSION COEFFICIENTS.
C
C ======================== ARGUMENTS / USAGE ===========================
C
C  INPUT:
C  1. [AKTE] 
C     DOUBLE PRECISION variable, the parameter kT/e needed for obtaining
C     the collision integral. 
C
C  OUTPUT:
C  1. [FUNCTION COLINT]
C     DOUBLE PRECISION variable, the value of the collision integral. 
C 
C=======================================================================
C
      DOUBLE PRECISION FUNCTION COLINT (AKTE)
C 
      DOUBLE PRECISION KTE, OMGD, SIGM, PARM, AKTE, CF
      COMMON /CINT/ KTE(81), OMGD(81), SIGM(5), PARM(5)
C
      IF (AKTE.LE.KTE(1)) THEN       ! kT/e < 0.3 => kT/e assigned 0.3
         COLINT = OMGD(1)
      ELSE IF (AKTE.GT.KTE(81)) THEN ! kT/e > 400 => kT/e assigned 400
         COLINT = OMGD(81)
      ELSE                           ! 0.3 > kT/e > 400 => use tables.
         DO 10 I=2,81
            IF (KTE(I-1).LT.AKTE .AND. AKTE.LE.KTE(I)) THEN
               CF     = (OMGD(I)-OMGD(I-1))/(KTE(I)-KTE(I-1))
               COLINT = OMGD(I-1) + CF*(AKTE-KTE(I-1))
            ENDIF
10       CONTINUE
      ENDIF
C
      RETURN     
      END

C=======================================================================
C
C *** MADM CODE
C *** BLOCK DATA
C *** PROVIDES THE INTEGRAL COLLISION TABLE NEEDED IN FUNCTION COLINT.
C     DATA WAS OBTAINED FROM TABLE 11-1, PAGE 524 OF:
C     REID R.C., SHERWOOD T.K., "THE PROPERTIES OF GASES AND LIQUIDS",
C     2ND EDITION, MCGRAW-HILL, NEW YORK, 1966.
C     
C     ALSO PROVIDES OTHER NECESSARY DATA FOR DIFFUSION COEFFICIENT
C     CALCULATIONS
C
C=======================================================================
C
      BLOCK DATA CLNT
C 
      include 'dynamic.inc'
      DOUBLE PRECISION KTE
      COMMON /CINT/ KTE(81), OMGD(81), SIGM(5), PARM(5)
C
      DATA KTE /0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65,
     &          0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05,
     &          1.10, 1.15, 1.20, 1.25, 1.30, 1.35, 1.40, 1.45,
     &          1.50, 1.55, 1.60, 1.65, 1.70, 1.75, 1.80, 1.85,
     &          1.90, 1.95, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50,
     &          2.60, 2.70, 2.80, 2.90, 3.00, 3.10, 3.20, 3.30,
     &          3.40, 3.50, 3.60, 3.70, 3.80, 3.90, 4.00, 4.10,
     &          4.20, 4.30, 4.40, 4.50, 4.60, 4.70, 4.80, 4.90,
     &          5.00, 6.00, 7.00, 8.00, 9.00, 10.0, 20.0, 30.0,
     &          40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100., 200.,
     &          400./
C
      DATA OMGD/2.662, 2.476, 2.318, 2.184, 2.066, 1.966, 1.877, 
     &          1.798, 1.729, 1.667, 1.612, 1.562, 1.517, 1.476, 
     &          1.439, 1.406, 1.375, 1.346, 1.320, 1.296, 1.273,
     &          1.253, 1.233, 1.215, 1.198, 1.182, 1.167, 1.153,
     &          1.140, 1.128, 1.116, 1.105, 1.094, 1.084, 1.075,
     &          1.057, 1.041, 1.026, 1.012, .9996, .9878, .9770,
     &          .9672, .9576, .9490, .9406, .9328, .9256, .9186,
     &          .9120, .9058, .8998, .8942, .8888, .8836, .8788,
     &          .8740, .8694, .8652, .8610, .8568, .8530, .8492,
     &          .8456, .8422, .8124, .7896, .7712, .7556, .7424,
     &          .6640, .6232, .5960, .5756, .5596, .5464, .5352,
     &          .5256, .5130, .4644, .4170/
C
C
      DATA SIGM/2.900D0, 3.300D0, 5.500D0, 3.339D0, 3.617D0/
C
      DATA PARM/558.3D0,  475.9D0, 77.3D0, 344.7D0,  97.0D0/
C
      END
