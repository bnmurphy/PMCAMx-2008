c*** AVGWAL
c
      subroutine avgwal(igrd,ncol,nrow,nlay,nspec,nspsa,dtime,
     &                    deltax,deltay,depth,tempk,press,saconc,conc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine calculates the average concentrations for the
c     WALL OF CELLS types of receptors.  It is not called unless
c     at least one receptor of this type is specified.  The average
c     is taken over all cells in the wall (including layers).  The 
c     concentraions in each cell are weighted by cell volume.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c         nrow   I  number of X cells in grid
c         ncol   I  number of Y cells in grid
c         nlay   I  number of layers
c         nspsa  I  number of tracer species
c         nspec  I  number of species in regular model
c         dtime  R  time step for present concs
c         deltax R  cell width in X direction        
c         deltay R  cell width in Y direction        
c         depth  R  layer depths
c         tempk  R  temperature
c         press  R  pressure
c         saconc R  tracer concentrations 
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     1.  Fixed bug by eliminating unused argument (avtime)
c     2.  Added grid number to recptors defined by cell index
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'bndary.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   nrow
      integer   ncol
      integer   nlay
      integer   nspsa
      integer   nspec
      real      dtime
      real      deltax(nrow)
      real      deltay
      real      depth(ncol,nrow,nlay)
      real      tempk(ncol,nrow,nlay)
      real      press(ncol,nrow,nlay)
      real      saconc(ncol,nrow,nlay,nspsa)
      real      conc(ncol,nrow,nlay,nspec)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   icl, jcl, kcl, ircp, ispc
      real      volume, volsum, sasum(MXTRSP), cncsum(MXSPEC)
      real      cncnow, cnvfac, timewt, cncvoc, cncnox, cnco3
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- calculte the time average weighting for running average
c
      timewt = dtime / dtout
c
c  --- loop over all receptors, except the hourly peak (recptor 1) ---
c
      do 10 ircp=2,nrecep
c
c  --- if receptor is not a WALL OF CELLS type, skip it ---
c
          if( idrcp(ircp) .NE. IDWAL ) goto 10
c
c   --- skip if receptor grid is not this grid ---
c
          if( igrdrcp(ircp) .NE. igrd ) goto 10
c
c  --- initialize the sums to zero ---
c
          volsum = 0.
          do 20 ispc=1,nspsa
             sasum(ispc) = 0.
   20     continue
          do 25 ispc=1,nspec
             cncsum(ispc) = 0.
   25     continue
c
c  --- loop over layers ---
c
          do 30 kcl=kwalbg(ircp),kwalnd(ircp)
c
c  ---- loop over all cells ---
c
              do 40 jcl=jwalbg(ircp),jwalnd(ircp)
                 do 50 icl=iwalbg(ircp),iwalnd(ircp)
c
c  --- slip the boundary cells ---
c
                    if( ibeg(jcl) .EQ. -999 ) goto 50
                    if( icl .LT. ibeg(jcl) .OR. 
     &                                icl .GT. iend(jcl) ) goto 50
c
c  ---- calculate the total volume for this cell and the 
c       conversion factor to ppm ---
c
                    volume = deltax(jcl)/1000. * deltay/1000. * 
     &                                              depth(icl,jcl,kcl)
                    volsum = volsum + volume
                    cnvfac = densfac*( 273./tempk(icl,jcl,kcl)*
     &                                         press(icl,jcl,kcl)/1013) 
c
c   --- loop over tracer species --
c
                    do 60 ispc=1,nspsa
                       sasum(ispc) = sasum(ispc) + volume * 
     &                                saconc(icl,jcl,kcl,ispc) / cnvfac
   60               continue
c
c   --- loop over regular mode species ----
c
                    do 70 ispc=1,nspec
                       cncsum(ispc) = cncsum(ispc) + volume * 
     &                                conc(icl,jcl,kcl,ispc) / cnvfac
   70               continue
c
c   --- next cell ---
c
   50            continue
   40         continue
c
c   --- next layer ---
c
   30     continue
c
c  ---- add concs to the running average ---
c
          do 80 ispc=1,nspsa
              cncnow = 0.
              if( volsum .GT. 0. ) cncnow = sasum(ispc) / volsum
              conrcp(ispc,ircp) = cncnow * timewt + conrcp(ispc,ircp)
   80     continue
          cncnox = 0.
          cncvoc = 0.
          cnco3 = 0.
          do 90 ispc=1,nspec
              cncnow = 0.
              if( volsum .GT. 0. ) cncnow = cncsum(ispc) / volsum
              if( lvocsp(ispc) ) then
                  cncvoc = cncvoc + cncnow  * crbnum(ispc)
              else if( lnoxsp(ispc) ) then
                  cncnox = cncnox + cncnow
              else if( lo3sp(ispc) ) then
                  cnco3 = cnco3 + cncnow
              endif
   90     continue
          rcpvoc(ircp) = cncvoc * timewt + rcpvoc(ircp)
          rcpnox(ircp) = cncnox * timewt + rcpnox(ircp)
          rcpo3(ircp) = cnco3 * timewt +  rcpo3(ircp)
c
c  --- next receptor ----
c
   10 continue
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
c
