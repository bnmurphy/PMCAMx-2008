c*** LOADDM
c
      subroutine loaddm(filflg,ncol,nrow,nlay,nddm,grsens,icl,jcl,kcl,
     &                  nsen,nspc,sens,convfac)

c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine loads the sensitivies which are stored in a
c     4-D array from/into a 2-D array in conclusion/preperation 
c     of the DDM chemistry routine.  The 2-D array contains the 
c     family of sensitivities for one cell and all modeled species.
c     The flag "filflg" determines wether the values in the 4-D
c     array are loaded into the 2-D array or vice-versa.  Sensitivities
c     are converted from umol/m3 to ppm for chemistry and back.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c       filflg  L  flag for determining which direction to fill
c                  .TRUE.  = put 2-D values into 4-D gridded array
c                  .FALSE. = put 4-D values into 2-D array
c       ncol    I  number of cells in X direction
c       nrow    I  number of cells in Y direction
c       nlay    I  number of layers 
c       nddm    I  number of total DDM species
c       grsens  R  4-D array of DDM sensitivities
c       icl     I  the X grid location of current cell
c       jcl     I  the Y grid location of current cell
c       kcl     I  the vertical grid location of current layer
c       nsen    I  numbder of DDM families
c       nspc    I  number of modeled species
c       sens    R  2-D array for this cell and species
c       convfac R  conversion from ppm to umol/m3
c
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      logical filflg
      integer ncol
      integer nrow
      integer nlay
      integer nddm
      real    grsens(ncol,nrow,nlay,nddm)
      integer icl
      integer jcl
      integer kcl
      integer nsen
      integer nspc
      real    sens(nsen,nspc)
      real    convfac
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   isen, ispc, iddm
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- task 1: load 4-D array
c
      if( filflg ) then
c
c  --- loop over the modeled species ---
c
        do ispc=1,nspc
c
c  --- loop over the number of DDM families ---
c
          do isen=1,nsen
c
c  --- calculate the index into the DDM species list ---
c
            iddm = iptddm(ispc)+isen-1
            grsens(icl,jcl,kcl,iddm) = sens(isen,ispc)*convfac
c
c  --- next family ---
c
          enddo
c
c  --- next species ---
c
        enddo
c
c  --- task 2: load 2-D array
c
      else
c
c  --- loop over the modeled species ---
c
        do ispc=1,nspc
c
c  --- loop over the number of DDM families ---
c
          do isen=1,nsen
c
c  --- calculate the index into the DDM species list ---
c
            iddm = iptddm(ispc)+isen-1
            sens(isen,ispc) = grsens(icl,jcl,kcl,iddm)/convfac
c
c  --- next family ---
c
          enddo
c
c  --- next species ---
c
        enddo
c
c  --- done ---
c
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
