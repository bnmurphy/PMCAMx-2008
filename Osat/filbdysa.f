c**** FILBDYSA.F
c
      subroutine filbdysa(igrid,nox,noy,noz,nspec,nspsa,conc,saconc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine fills one hour of boundary conditions and calculates 
c   the NOx and VOC levels.  It then places these concentrations in the 
c   appropriate place in the gridded array used for tracer concentrations.  
c   The O3 concentrations are placed into the concentration arrays
c   for the Ozone tracer species.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c      Argument description:
c           saconc   R  array of tracer concentrations
c       Inputs:
c           igrid    I  grid number of this grid
c           nox      I  number of cells in X direction
c           noy      I  number of cells in Y direction
c           noz      I  number of layers 
c           conc     R  array of regular model concentrations
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     06/06/96   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'bndary.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   igrid
      integer   nox
      integer   noy
      integer   noz
      integer   nspec
      integer   nspsa
      real      conc(nox,noy,noz,nspec)
      real      saconc(nox,noy,noz,nspsa)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      ioff, icl, jcl, izcl, ispc
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- read concentrations for each species ----
c
      do 10 ispc=1,nspec
c
c   --- if the species is a not modelled or not a VOC or NOx or O3
c       species skip it ---
c
          if( .NOT. lvocsp(ispc) .AND. .NOT. lnoxsp(ispc) 
     &                                 .AND. .NOT. lo3sp(ispc) ) goto 10
c
c  --- do the WEST and EAST boundaries ----
c
          do 20 jcl=1,noy
c
c  --- West boundary ---
c
             if( igrid .EQ. 1 ) then
                if( ibeg(jcl) .EQ. -999 ) goto 20
                icl = ibeg(jcl) - 1
             else
                icl = 1
             endif
             do 30 izcl=1,noz
                if( conc(icl,jcl,izcl,ispc) .LE. 0 ) goto 30
c
c   --- if stratifying by boundary, put in seperate position ---
c
                if( lbndry ) then
                   ioff = IDXBWS
                else
                   ioff = 1
                endif
c
c   --- VOC species ---
c
                if( lvocsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptvoc+ioff) = 
     &                      saconc(icl,jcl,izcl,iptvoc+ioff) + 
     &                           conc(icl,jcl,izcl,ispc) * crbnum(ispc)
c
c   --- NOx species ----
c
                 else if( lnoxsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptnox+ioff) = 
     &                      saconc(icl,jcl,izcl,iptnox+ioff) + 
     &                                          conc(icl,jcl,izcl,ispc)
c
c   --- Ozone species, put half into each ozone tracer ---
c
                 else if( lo3sp(ispc) ) then
                    saconc(icl,jcl,izcl,ipto3n+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3n+ioff) + 
     &                                   conc(icl,jcl,izcl,ispc) / 2.0
                    saconc(icl,jcl,izcl,ipto3v+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3v+ioff) + 
     &                                   conc(icl,jcl,izcl,ispc) / 2.0
                 endif
   30        continue
c
c  --- East boundary ---
c
             if( igrid .EQ. 1 ) then
                icl = iend(jcl) + 1
             else
                icl = nox
             endif
             do 40 izcl=1,noz
                if( conc(icl,jcl,izcl,ispc) .LE. 0. ) goto 40
c
c   --- if stratifying by boundary, put in seperate position ---
c
                if( lbndry ) then
                   ioff = IDXBES
                else
                   ioff = 1
                endif
c
c   --- VOC species ---
c
                if( lvocsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptvoc+ioff) = 
     &                      saconc(icl,jcl,izcl,iptvoc+ioff) + 
     &                           conc(icl,jcl,izcl,ispc) * crbnum(ispc)
c
c   --- NOx species ----
c
                 else if( lnoxsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptnox+ioff) = 
     &                      saconc(icl,jcl,izcl,iptnox+ioff) + 
     &                                         conc(icl,jcl,izcl,ispc)
c
c   --- Ozone species, put half into each ozone tracer ---
c
                 else if( lo3sp(ispc) ) then
                    saconc(icl,jcl,izcl,ipto3n+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3n+ioff) + 
     &                                  conc(icl,jcl,izcl,ispc) / 2.0
                    saconc(icl,jcl,izcl,ipto3v+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3v+ioff) + 
     &                                  conc(icl,jcl,izcl,ispc) / 2.0
                 endif
   40        continue
   20     continue
c
c  --- do the SOUTH and NORTH boundaries ----
c
          do 50 icl=1,nox
c
c  --- South boundary ---
c
             if( igrid .EQ. 1 ) then
                 if( jbeg(icl) .EQ. -999 ) goto 50
                 jcl = jbeg(icl) - 1
             else
                 jcl = 1
             endif
             do 60 izcl=1,noz
                if( conc(icl,jcl,izcl,ispc) .LE. 0. ) goto 60
c
c   --- stratifying by boundary, put in seperate position ---
c
                if( lbndry ) then
                    ioff = IDXBST
                else
                    ioff = 1
                endif
c
c   --- VOC species ---
c
                if( lvocsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptvoc+ioff) = 
     &                      saconc(icl,jcl,izcl,iptvoc+ioff) + 
     &                           conc(icl,jcl,izcl,ispc) * crbnum(ispc)
c
c   --- NOx species ----
c
                 else if( lnoxsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptnox+ioff) = 
     &                      saconc(icl,jcl,izcl,iptnox+ioff) + 
     &                                         conc(icl,jcl,izcl,ispc)
c
c   --- Ozone species, put half into each ozone tracer ---
c
                 else if( lo3sp(ispc) ) then
                    saconc(icl,jcl,izcl,ipto3n+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3n+ioff) + 
     &                                  conc(icl,jcl,izcl,ispc) / 2.0
                       saconc(icl,jcl,izcl,ipto3v+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3v+ioff) + 
     &                                  conc(icl,jcl,izcl,ispc) / 2.0
                 endif
   60        continue
c
c  --- North boundary ---
c
             if( igrid .EQ. 1 ) then
                jcl = jend(icl) + 1
             else
                jcl = noy
             endif
             do 70 izcl=1,noz
                if( conc(icl,jcl,izcl,ispc) .LE. 0. ) goto 70
c
c   --- if stratifying by boundary, put in seperate position ---
c
                if( lbndry ) then
                    ioff = IDXBNT
                else
                    ioff = 1
                endif
c
c   --- VOC species ---
c
                if( lvocsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptvoc+ioff) = 
     &                      saconc(icl,jcl,izcl,iptvoc+ioff) + 
     &                           conc(icl,jcl,izcl,ispc) * crbnum(ispc)
c
c   --- NOx species ----
c
                 else if( lnoxsp(ispc) ) then
                    saconc(icl,jcl,izcl,iptnox+ioff) = 
     &                      saconc(icl,jcl,izcl,iptnox+ioff) + 
     &                                         conc(icl,jcl,izcl,ispc)
c
c   --- Ozone species, put half into each ozone tracer ---
c
                 else if( lo3sp(ispc) ) then
                    saconc(icl,jcl,izcl,ipto3n+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3n+ioff) + 
     &                                 conc(icl,jcl,izcl,ispc) / 2.0
                    saconc(icl,jcl,izcl,ipto3v+ioff) = 
     &                      saconc(icl,jcl,izcl,ipto3v+ioff) + 
     &                                 conc(icl,jcl,izcl,ispc) / 2.0
                 endif
   70        continue
   50     continue
c
c  --- next species --
c
   10 continue
c
c  --- return to the calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
