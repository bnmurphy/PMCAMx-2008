      REAL FUNCTION HLCONST ( NAME, TEMP )

C-----------------------------------------------------------------------
C
C  FUNCTION: return the Henry's law constant for the specified substance
C            at the given temperature
C
C  revision history:
C    who        when           what
C  ---------  --------  -------------------------------------
C  S.Roselle  08/15/97  code written for Models-3
C  J.Gipson   06/18/01  added Henry's Law constants 50-55 for saprc99
C  W.Hutzell  07/03/01  added Henry's Law constants 56-57 for Atrazine
C                       and the daughter products from Atrazine and OH
C                       reactions.
C-----------------------------------------------------------------------

      IMPLICIT NONE

C...........INCLUDES and their descriptions

c      INCLUDE SUBST_IODECL    ! I/O definitions and declarations

C...........PARAMETERS and their descriptions:

      INTEGER       MXSPCS              ! Number of substances
      PARAMETER   ( MXSPCS = 58 )

C...........ARGUMENTS and their descriptions

      CHARACTER*(*) NAME                ! name of substance
      REAL          TEMP                ! temperature (K)

C...........SCRATCH LOCAL VARIABLES and their descriptions:

      CHARACTER*16  SUBNAME( MXSPCS )   ! list of substance names
      SAVE          SUBNAME

      INTEGER       SPC                 ! species index

      REAL          A( MXSPCS )         ! Henry's law constants at 298.15K (M/atm)
      SAVE          A
      REAL          E( MXSPCS )         ! enthalpy (like activation energy) (K)
      SAVE          E

C...........EXTERNAL FUNCTIONS and their descriptions:

      INTEGER       HLINDEX

c***********************************************************************
      DATA SUBNAME(  1), A(  1), E(  1) 
     &     / 'O3              ', 1.2E-02, 2.7E+03 /  ! Chameides 1984
      DATA SUBNAME(  2), A(  2), E(  2) 
     &     / 'HO2             ', 9.0E+03, 0.0E+00 /  ! Chameides 1984
      DATA SUBNAME(  3), A(  3), E(  3) 
     &     / 'H2O2            ', 9.7E+04, 6.6E+03 /  ! Chameides 1984
      DATA SUBNAME(  4), A(  4), E(  4) 
     &     / 'NH3             ', 5.8E+01, 4.1E+03 /  ! Chameides 1984
      DATA SUBNAME(  5), A(  5), E(  5) 
     &     / 'NO              ', 1.9E-03, 1.4E+03 /  ! Lide and Frederikse 1995
      DATA SUBNAME(  6), A(  6), E(  6) 
     &     / 'NO2             ', 1.2E-02, 2.5E+03 /  ! Chameides 1984
      DATA SUBNAME(  7), A(  7), E(  7) 
     &     / 'NO3             ', 1.2E+01, 1.9E+03 /  ! Chameides 1984
      DATA SUBNAME(  8), A(  8), E(  8) 
     &     / 'N2O5            ', 1.0E+30, 0.0E+00 /  ! "inf" Sander and Crutzen 1996
      DATA SUBNAME(  9), A(  9), E(  9) 
     &     / 'HNO2            ', 4.9E+01, 4.8E+03 /  ! Chameides 1984
      DATA SUBNAME( 10), A( 10), E( 10) 
     &     / 'HNO3            ', 2.6E+06, 8.7E+03 /  ! Chameides 1984
      DATA SUBNAME( 11), A( 11), E( 11) 
     &     / 'HNO4            ', 2.0E+04, 0.0E+00 /  ! Jacob et al. 1989
      DATA SUBNAME( 12), A( 12), E( 12) 
     &     / 'SO2             ', 1.2E+00, 3.1E+03 /  ! Chameides 1984
      DATA SUBNAME( 13), A( 13), E( 13) 
     &     / 'H2SO4           ', 1.0E+30, 0.0E+00 /  ! infinity
      DATA SUBNAME( 14), A( 14), E( 14) 
     &     / 'METHANE         ', 1.3E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 15), A( 15), E( 15) 
     &     / 'ETHANE          ', 2.0E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 16), A( 16), E( 16) 
     &     / 'PROPANE         ', 1.4E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 17), A( 17), E( 17) 
     &     / 'BUTANE          ', 1.1E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 18), A( 18), E( 18) 
     &     / 'PENTANE         ', 8.1E-04, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 19), A( 19), E( 19) 
     &     / 'HEXANE          ', 6.0E-04, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 20), A( 20), E( 20) 
     &     / 'OCTANE          ', 3.4E-04, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 21), A( 21), E( 21) 
     &     / 'NONANE          ', 2.0E-04, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 22), A( 22), E( 22) 
     &     / 'DECANE          ', 1.4E-04, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 23), A( 23), E( 23) 
     &     / 'ETHENE          ', 4.7E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 24), A( 24), E( 24) 
     &     / 'PROPENE         ', 4.8E-03, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 25), A( 25), E( 25) 
     &     / 'ISOPRENE        ', 1.3E-02, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 26), A( 26), E( 26) 
     &     / 'ACETYLENE       ', 4.1E-02, 1.8E+03 /  ! Wilhelm et al. 1977
      DATA SUBNAME( 27), A( 27), E( 27) 
     &     / 'BENZENE         ', 1.8E-01, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 28), A( 28), E( 28) 
     &     / 'TOLUENE         ', 1.5E-01, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 29), A( 29), E( 29) 
     &     / 'O-XYLENE        ', 2.0E-01, 0.0E+00 /  ! Mackay and Shiu 1981
      DATA SUBNAME( 30), A( 30), E( 30) 
     &     / 'METHANOL        ', 2.2E+02, 0.0E+00 /  ! Snider and Dawson 1985
      DATA SUBNAME( 31), A( 31), E( 31) 
     &     / 'ETHANOL         ', 1.6E+02, 0.0E+00 /  ! Betterton 1992
      DATA SUBNAME( 32), A( 32), E( 32) 
     &     / '2-CRESOL        ', 8.2E+02, 0.0E+00 /  ! Betterton 1992
      DATA SUBNAME( 33), A( 33), E( 33) 
     &     / '4-CRESOL        ', 1.3E+02, 0.0E+00 /  ! Betterton 1992
      DATA SUBNAME( 34), A( 34), E( 34) 
     &     / 'METHYLHYDROPEROX', 3.1E+02, 5.2E+03 /  ! O'Sullivan et al. 1996
      DATA SUBNAME( 35), A( 35), E( 35) 
     &     / 'FORMALDEHYDE    ', 7.0E+03, 6.4E+03 /  ! Chameides 1984
      DATA SUBNAME( 36), A( 36), E( 36) 
     &     / 'ACETALDEHYDE    ', 1.3E+01, 5.7E+03 /  ! Benkelberg et al. 1995
      DATA SUBNAME( 37), A( 37), E( 37) 
     &     / 'GENERIC_ALDEHYDE', 4.2E+03, 0.0E+00 /  ! Graedel and Goldberg 1983
      DATA SUBNAME( 38), A( 38), E( 38) 
     &     / 'GLYOXAL         ', 3.6E+05, 0.0E+00 /  ! Zhou and Mopper 1990
      DATA SUBNAME( 39), A( 39), E( 39) 
     &     / 'ACETONE         ', 3.2E+01, 5.8E+03 /  ! Betterton 1991
      DATA SUBNAME( 40), A( 40), E( 40) 
     &     / 'FORMIC_ACID     ', 3.7E+03, 5.7E+03 /  ! Chameides 1984
      DATA SUBNAME( 41), A( 41), E( 41) 
     &     / 'ACETIC_ACID     ', 5.2E+03, 0.0E+00 /  ! Johnson et al. 1996
      DATA SUBNAME( 42), A( 42), E( 42) 
     &     / 'METHYL_GLYOXAL  ', 3.2E+04, 0.0E+00 /  ! Zhou and Mopper 1990
      DATA SUBNAME( 43), A( 43), E( 43) 
     &     / 'CO              ', 9.5E-04, 1.3E+03 /  ! Wilhelm et al. 1977
      DATA SUBNAME( 44), A( 44), E( 44) 
     &     / 'CO2             ', 3.1E-02, 2.4E+03 /  ! Chameides 1984
      DATA SUBNAME( 45), A( 45), E( 45) 
     &     / 'PAN             ', 2.9E+00, 5.9E+03 /  ! Pandis and Seinfeld 1989
      DATA SUBNAME( 46), A( 46), E( 46) 
     &     / 'MPAN            ', 1.7E+00, 0.0E+00 /  ! Kames and Schurath 1995
      DATA SUBNAME( 47), A( 47), E( 47) 
     &     / 'OH              ', 3.0E+01, 4.5E+03 /  ! Hanson et al. 1992
      DATA SUBNAME( 48), A( 48), E( 48) 
     &     / 'METHYLPEROXY_RAD', 2.0E+03, 6.6E+03 /  ! Lelieveld and Crutzen 1991
      DATA SUBNAME( 49), A( 49), E( 49) 
     &     / 'PEROXYACETIC_ACI', 8.4E+02, 5.3E+03 /  ! O'Sullivan et al. 1996
      DATA SUBNAME( 50), A( 50), E( 50) 
     &     / 'PROPANOIC_ACID  ', 5.7E+03, 0.0E+00 /  ! Kahn et al. 1995
      DATA SUBNAME( 51), A( 51), E( 51) 
     &     / '2-NITROPHENOL   ', 7.0E+01, 4.6E+03 /  ! USEPA 1982
      DATA SUBNAME( 52), A( 52), E( 52) 
     &     / 'PHENOL          ', 1.9E+03, 7.3E+03 /  ! USEPA 1982
      DATA SUBNAME( 53), A( 53), E( 53) 
     &     / 'BIACETYL        ', 7.4E+01, 5.7E+03 /  ! Betteron 1991
      DATA SUBNAME( 54), A( 54), E( 54) 
     &     / 'BENZALDEHYDE    ', 4.2E+01, 4.6E+03 /  ! Zhou and Mopper 1990
      DATA SUBNAME( 55), A( 55), E( 55) 
     &     / 'PINENE          ', 4.9E-02, 0.0E+00 /  ! Karl and Lindinger 1997
      DATA SUBNAME( 56), A( 56), E( 56) 
     &     / 'ATRA            ', 4.1E+05, 6.0E+03 /  ! CIBA Corp (1989) and Scholtz (1999)
      DATA SUBNAME( 57), A( 57), E( 57) 
     &     / 'DATRA           ', 4.1E+05, 6.0E+03 /  ! assumed same as Atrazine
      DATA SUBNAME( 58), A( 58), E( 58) 
     &     / 'ADIPIC_ACID     ', 2.0E+08, 0.0E+00 /  ! Saxena and Hildemann (1996)


C-----------------------------------------------------------------------
C  begin body of subroutine HLCONST

      SPC = HLINDEX( NAME, MXSPCS, SUBNAME )

      IF ( SPC .GT. 0 ) THEN
        HLCONST = A( SPC ) * EXP( E( SPC )
     &          * ( ( 298.0 - TEMP) / ( 298.0 * TEMP) ) )
      ELSE
        HLCONST = 0.0
      END IF

      RETURN
      END              
