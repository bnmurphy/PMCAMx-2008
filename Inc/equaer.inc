C=======================================================================
C
C *** THIS FILE CONTAINS THE DECLARATIONS OF THE GLOBAL CONSTANTS
C     AND VARIABLES. 
C
C
C=======================================================================
C
      PARAMETER (NPAIR=13,NIONS=7,NCOMP=5,NGASAQ=3)
      parameter (itmaxeq=1000)	! maximum iteration in equilibrium codes
C
C *** VARIABLES FOR LIQUID AEROSOL PHASE *******************************
C
      DOUBLE PRECISION MOLAL, MOLALR, M0
      REAL IONIC
      LOGICAL CALAOU, CALAIN, FRST, DRYF
      COMMON /IONS/ MOLAL(NIONS), MOLALR(NPAIR), GAMA(NPAIR), ZZ(NPAIR),
     &              Z(NIONS),     GAMOU(NPAIR),  GAMIN(NPAIR),M0(NPAIR),
     &              GASAQ(NGASAQ),
     &              EPSACT,       COH,           CHNO3,       CHCL,         
     &              WATER,        IONIC,         IACALC,      
     &              FRST,         CALAIN,        CALAOU,      DRYF
C
C *** VARIABLES FOR SOLID AEROSOL PHASE ********************************
C
      COMMON /SALT/ CH2SO4,  CNH42S4, CNH4HS4, CNACL,   CNA2SO4, 
     &              CNANO3,  CNH4NO3, CNH4CL,  CNAHSO4, CLC
C
C
C *** VARIABLES FOR GAS PHASE ******************************************
C
      COMMON /GAS / GNH3, GHNO3, GHCL 

C
C *** EQUILIBRIUM CONSTANTS ********************************************
C
      COMMON /EQUK/ XK1, XK2, XK3, XK4, XK5, XK6, XK7, XK8, XK9, XK10,
     &              XK11,XK12,XK13,XK14,XKW, XK21,XK22,XK31,XK32,XK41,
     &              XK42
C
c     added by bkoo (03/05/02)
      integer metstbl
      common /inpt/ realdum(2*NCOMP+2), intdum, metstbl
C
C *** END OF INCLUDE FILE **********************************************
C
