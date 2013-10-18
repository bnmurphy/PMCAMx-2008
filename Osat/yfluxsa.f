c*** YFLUXSA
c
      subroutine yfluxsa(ncol,nrow,nlay,nspc,saconc,
     &                    icell,jclbeg,jclend,kcell,dx,dz,mscl,
     &                    fluxo3,cnco3,fluxnox,cncnox,fluxvoc,cncvoc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine updates the tracer concentrations arrays in the
c     common block variables by applying the appropriate flux value
c     calculated from the regulat model transport routine.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       ncol     I  number of columns in the grid
c       nrow     I  number of rows in the grid
c       nlay     I  number of layers in the grid
c       nspc     I  number of species
c       saconc   R  3-D concentration arrays
c       icell    I  cell index in the X direction
c       jclbeg   I  beginning cell in the Y direction
c       jclend   I  ending cell in the Y direction
c       kcell    I  the vertical grid location of current layer
c       dx       R  width of cell in X direction
c       dz       R  depth of layer
c       mscl     R  map scale factor
c       fluxo3   R  flux for the O3 species
c       cnco3    R  regular model concentration for the O3 species
c       fluxnox  R  flux for the NOx species
c       cncnox   R  regular model concentration for the NOx species
c       fluxvoc  R  flux for the VOC species
c       cncvoc   R  regular model concentration for the VOC species
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c     10/02/98   --gwilson--  Orignal development
c     10/30/01   --cemery--   Added map scale factor to flux div
c                             calculation 
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
      integer     ncol
      integer     nrow
      integer     nlay
      integer     nspc
      real        saconc(ncol,nrow,nlay,nspc)
      integer     icell
      integer     jclbeg
      integer     jclend
      integer     kcell
      real        dx(nrow)
      real        dz(ncol,nrow,nlay)
      real        mscl(ncol,nrow)
      real        fluxo3(MX1D)
      real        cnco3(MX1D)
      real        fluxnox(MX1D)
      real        cncnox(MX1D)
      real        fluxvoc(MX1D)
      real        cncvoc(MX1D)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   jcl, i1d , ispc
      real      flux(MX1D), cnctmp(MX1D)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- NOx species ---
c
      do ispc=iptnox,iptvoc-1
c
c   --- normalize the flux on the boundary ---
c
        i1d = 1
        if( fluxnox(i1d) .GE. 0 ) then
           flux(i1d) = saconc(icell,jclbeg-1,kcell,ispc) * 
     &                                       fluxnox(i1d) / cncnox(i1d)
        else
           flux(i1d) = saconc(icell,jclbeg,kcell,ispc) * 
     &                                     fluxnox(i1d) / cncnox(i1d+1)
        endif
c
c   --- loop over all interior cells in this row ---
c
        do jcl=jclbeg,jclend
          i1d = i1d + 1 
c
c   ---  apply flux to the NOx tracer species ---
c
           if ( fluxnox(i1d) .GE. 0. ) then
                flux(i1d) = saconc(icell,jcl,kcell,ispc) *
     &                                    fluxnox(i1d) / cncnox(i1d)
           else
                flux(i1d) = saconc(icell,jcl+1,kcell,ispc) *
     &                                     fluxnox(i1d) / cncnox(i1d+1)
           endif
           cnctmp(i1d) = saconc(icell,jcl,kcell,ispc) - 
     &       (mscl(icell,jcl)*mscl(icell,jcl)*
     &        (flux(i1d) - flux(i1d-1)) / dx(jcl) / dz(icell,jcl,kcell))
        enddo
c
c  --- put new concentrations back into 3-D array ---
c
        i1d = 1
        do jcl=jclbeg,jclend
          i1d = i1d + 1
          saconc(icell,jcl,kcell,ispc) = AMAX1(cnctmp(i1d),BNDLPT)
        enddo
      enddo
c
c   --- VOC species ---
c
      do ispc = iptvoc,ipto3n-1
c
c   --- normalize the flux on the boundary ---
c
        i1d = 1
        if( fluxvoc(i1d) .GE. 0 ) then
            flux(i1d) = saconc(icell,jclbeg-1,kcell,ispc) * 
     &                                      fluxvoc(i1d) / cncvoc(i1d)
        else
            flux(i1d) = saconc(icell,jclbeg,kcell,ispc) * 
     &                                    fluxvoc(i1d) / cncvoc(i1d+1)
        endif
c
c   --- loop over all interior cells in this row ---
c
        do jcl=jclbeg,jclend
           i1d = i1d + 1 
c
c   ---  apply flux to the VOC tracer species ---
c
           if ( fluxvoc(i1d) .GE. 0. ) then
               flux(i1d) = saconc(icell,jcl,kcell,ispc) *
     &                                    fluxvoc(i1d) / cncvoc(i1d)
           else
               flux(i1d) = saconc(icell,jcl+1,kcell,ispc) *
     &                                   fluxvoc(i1d) / cncvoc(i1d+1)
           endif
           cnctmp(i1d) = saconc(icell,jcl,kcell,ispc) - 
     &       (mscl(icell,jcl)*mscl(icell,jcl)*
     &        (flux(i1d) - flux(i1d-1)) / dx(jcl) / dz(icell,jcl,kcell))
        enddo
c
c  --- put new concentrations back into 3-D array ---
c
        i1d = 1
        do jcl=jclbeg,jclend
          i1d = i1d + 1
          saconc(icell,jcl,kcell,ispc) = AMAX1(cnctmp(i1d),BNDLPT)
        enddo
      enddo
c
c   ---  O3 species ---
c
      do ispc = ipto3n,ipttim-1
c
c   --- normalize the flux on the boundary ---
c
        i1d = 1
        if( fluxo3(i1d) .GE. 0 ) then
            flux(i1d) = saconc(icell,jclbeg-1,kcell,ispc) * 
     &                                        fluxo3(i1d) / cnco3(i1d)
        else
            flux(i1d) = saconc(icell,jclbeg,kcell,ispc) * 
     &                                      fluxo3(i1d) / cnco3(i1d+1)
        endif
c
c   --- loop over all interior cells in this row ---
c
        do jcl=jclbeg,jclend
           i1d = i1d + 1 
           if( fluxo3(i1d) .GE. 0. ) then
                flux(i1d) = saconc(icell,jcl,kcell,ispc) *
     &                                    fluxo3(i1d) / cnco3(i1d)
           else
                flux(i1d) = saconc(icell,jcl+1,kcell,ispc) *
     &                                      fluxo3(i1d) / cnco3(i1d+1)
           endif
           cnctmp(i1d) = saconc(icell,jcl,kcell,ispc) - 
     &       (mscl(icell,jcl)*mscl(icell,jcl)*
     &        (flux(i1d) - flux(i1d-1)) / dx(jcl) / dz(icell,jcl,kcell))
        enddo
c
c  --- put new concentrations back into 3-D array ---
c
        i1d = 1
        do jcl=jclbeg,jclend
          i1d = i1d + 1
          saconc(icell,jcl,kcell,ispc) = AMAX1(cnctmp(i1d),BNDLPT)
        enddo
      enddo
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
