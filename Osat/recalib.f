c*** RECALIB
c
      subroutine recalib(ncol,nrow,nlay,nspc,saconc,
     &                                icl,jcl,kcl,modo3,modnox,modvoc)

c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine recalibrates the tracer concentrations back to
c     the levels calculated in the regular model.  This is necessary
c     because the numerical algorithms used in the transport steps
c     cause the tracer species to stray from the regular model.
c     This is due to the sharp gradients at model boundary cells
c     and between source regions.  The regular model does not stratify
c     the concentrations by source region so these gradients do not 
c     appear in the regular model.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c       ncol    number of cells in X direction
c       nrow    number of cells in Y direction
c       nlay    number of layers 
c       nspc    number of tracer species
c       saconc  tracer concentration array
c       icl     I  the X grid location of current cell
c       jcl     I  the Y grid location of current cell
c       kcl     I  the vertical grid location of current layer
c       modo3   R  regular model O3 concentrations 
c       modnox  R  regular model NOx concentrations 
c       modvoc  R  regular model VOC concentrations 
c                  (carbon weighted sum of VOC species)
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
      integer*4 ncol
      integer*4 nrow
      integer*4 nlay
      integer*4 nspc
      real*4    saconc(ncol,nrow,nlay,nspc)
      integer*4 icl
      integer*4 jcl
      integer*4 kcl
      real*4    modo3
      real*4    modnox
      real*4    modvoc
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   i
      real      trvoc, tro3, trnox
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- calculate the sum of the tracer species ----
c
      trnox = 0.
      do 10 i=iptnox,iptvoc-1
         trnox = trnox + saconc(icl,jcl,kcl,i)
   10 continue
      trvoc = 0.
      do 20 i=iptvoc,ipto3n-1
         trvoc = trvoc + saconc(icl,jcl,kcl,i)
   20 continue
      tro3 = 0.
      do 30 i=ipto3n,ipttim-1
         tro3 = tro3 + saconc(icl,jcl,kcl,i)
   30 continue
c
c   --- recalibrate the NOx tracers, if necessary ---
c
      if( trnox .NE. modnox .AND. trnox .GT. 0 ) then
         do 40 i=iptnox,iptvoc-1
            saconc(icl,jcl,kcl,i) = saconc(icl,jcl,kcl,i) + 
     &                 (modnox-trnox) * saconc(icl,jcl,kcl,i) / trnox
            saconc(icl,jcl,kcl,i) = MAX(saconc(icl,jcl,kcl,i),BNDLPT)
   40    continue
      endif
c
c   --- recalibrate the VOC tracers, if necessary ---
c
      if( trvoc .NE. modvoc .AND. trvoc .GT. 0 ) then
         do 50 i=iptvoc,ipto3n-1
            saconc(icl,jcl,kcl,i) = saconc(icl,jcl,kcl,i) + 
     &                 (modvoc-trvoc) * saconc(icl,jcl,kcl,i) / trvoc
            saconc(icl,jcl,kcl,i) = MAX(saconc(icl,jcl,kcl,i),BNDLPT)
   50    continue
      endif
c
c   --- recalibrate the O3 tracers, if necessary ---
c
      if( tro3 .NE. modo3 .AND. tro3 .GT. 0 ) then
         do 60 i=ipto3n,ipttim-1
            saconc(icl,jcl,kcl,i) = saconc(icl,jcl,kcl,i) + 
     &                    (modo3-tro3) * saconc(icl,jcl,kcl,i) / tro3
            saconc(icl,jcl,kcl,i) = MAX(saconc(icl,jcl,kcl,i),BNDLPT)
   60    continue
      endif
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
      return
      end
