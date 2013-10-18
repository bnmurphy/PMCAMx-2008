      INTEGER FUNCTION HLINDEX (NAME, N, NLIST)

C***********************************************************************
C
C  FUNCTION:
C
C    Searches for NAME in list NLIST and returns the subscript
C    (1...N) at which it is found, or returns 0 when NAME not
C    found in NLIST
C
C  PRECONDITIONS REQUIRED:  none
C
C  SUBROUTINES AND FUNCTIONS CALLED:  none
C
C  REVISION HISTORY:
C
C    5/88   Modified for ROMNET
C    9/94   Modified for Models-3 by CJC
C
C***********************************************************************

      IMPLICIT NONE
 
C.......   Arguments and their descriptions:

      CHARACTER*(*) NAME        !  Character string being searched for
      INTEGER       N           !  Length of array to be searched
      CHARACTER*(*) NLIST(*)    !  array to be searched

C.......   Local variable:

      INTEGER       I   !  loop counter

C.....................................................................
C.......   begin body of INDEX1()

      DO 100 I = 1, N

          IF ( NAME .EQ. NLIST( I ) ) THEN    ! Found NAME in NLIST
              HLINDEX = I
              RETURN
          ENDIF

100   CONTINUE

      HLINDEX = 0        !  not found
      RETURN

      END
