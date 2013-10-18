c**** FILVDSA.F
c
      subroutine filvdsa(nox,noy,noz,nspcs,nspsa,
     &                                     conc,vdep,vdepsa,vdeprt)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine fills the depostion velocity from the regular model
c   depostion velocities.  The tracer vdeps are a concentration
c   weighted average of the regular model vdeps.
c   
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c      Argument description:
c       Outputs:
c           vdepsa   R  depostion velocity for tracer species
c       Inputs:
c           nox      I  number of X cells in the grid
c           noy      I  number of Y cells in the grid
c           noz      I  number of layers in the grid
c           nspcs    I  number of species in the grid
c           nspsa    I  number of tracer species
c           conc     R  regular model concentrations
c           vdep     R  diffusion velocity for regular model species
c           vdeprt   R  diffusion velocities for RTRAC
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     06/06/96   --gwilson--    Original development
c     10/11/97   --gwilson--    Removed unused argument SACONC
c     01/30/02   --gwilson--    Added code for RTRAC probing tool
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'chmstry.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   nox
      integer   noy
      integer   noz
      integer   nspcs
      integer   nspsa
      real      vdep(nox,noy,nspcs)
      real      conc(nox,noy,noz,nspcs)
      real      vdepsa(nox,noy,nspsa)
      real      vdeprt(nox,noy,MXTRSP)
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c     IXNOX   I  index for NOx in local arrays
c     IXVOC   I  index for VOC in local arrays
c
      integer   IXNOX
      integer   IXVOC
c
      parameter( IXNOX = 1 )
      parameter( IXVOC = 2 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      icl, jcl, ispc, idxo3
      real         consum(MXCOLA,MXROWA,IXVOC)
      real         vdsum(MXCOLA,MXROWA,IXVOC)
      real         valnox, valvoc
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- call routine to zero out the arrays ---
c
      call zeros(vdepsa,nox*noy*nspsa)
      call zeros(consum,MXCOLA*MXROWA*IXVOC)
      call zeros(vdsum,MXCOLA*MXROWA*IXVOC)
c
c   --- if doing RTRAC, load from the gridded arrays 
c       filled by DRYDEPRT ---
c
      if( tectyp .EQ. RTRAC ) then
          do ispc=1,nspsa
             do icl=1,nox
                do jcl=1,noy
                  vdepsa(icl,jcl,ispc) = vdeprt(icl,jcl,ispc)
                enddo
             enddo
          enddo
          goto 9999
      endif
c
c  ---- read concentrations for each species ----
c
      do 10 ispc=1,nspcs
c
c   --- if the species is O3, set the index and get next species ---
c
          if( lo3sp(ispc) ) then
             idxo3 = ispc
             goto 10
          endif
c
c   --- if the species is a not modelled or not a VOC or NOx or O3
c       species skip it ---
c
          if( .NOT. lvocsp(ispc) .AND. .NOT. lnoxsp(ispc) ) goto 10
c
c   --- put concentractions in arrays ---
c
           do 20 icl=1,nox
c
c  --- set the beginning and ending of the interior of grid ---
c
c
c   --- loop over cells ---
c
             do 30 jcl=1,noy
c
c   --- VOC species ---
c
                if( lvocsp(ispc) ) then
                   consum(icl,jcl,IXVOC) = consum(icl,jcl,IXVOC) + 
     &                         conc(icl,jcl,1,ispc) * crbnum(ispc)
                   vdsum(icl,jcl,IXVOC) = vdsum(icl,jcl,IXVOC) + 
     &                            conc(icl,jcl,1,ispc) * crbnum(ispc) * 
     &                                               vdep(icl,jcl,ispc)
c
c   --- NOx species ----
c
                else if( lnoxsp(ispc) ) then
                   consum(icl,jcl,IXNOX) = consum(icl,jcl,IXNOX) + 
     &                                              conc(icl,jcl,1,ispc)
                   vdsum(icl,jcl,IXNOX) = vdsum(icl,jcl,IXNOX) + 
     &                         conc(icl,jcl,1,ispc) * vdep(icl,jcl,ispc)
                endif
c
c  --- next cell ---
c
   30         continue
   20     continue
c
c  --- next species --
c
   10 continue
c
c  --- loop over cells and calculate the tracer vdeps ---
c
      do 50 icl=1,nox
         do 60 jcl=1,noy
c
c   --- caclulate NOX vdeps based on regular model conc weighted vdeps ----
c
            valnox = 0.
            if( consum(icl,jcl,IXNOX) .NE. 0. ) valnox =
     &                    vdsum(icl,jcl,IXNOX) /  consum(icl,jcl,IXNOX)
            do 70 ispc=iptnox,iptvoc-1
               vdepsa(icl,jcl,ispc) = valnox
   70       continue
c
c   --- caclulate VOC vdeps based on regular model conc weighted vdeps ----
c
            valvoc = 0.
            if( consum(icl,jcl,IXVOC) .NE. 0. ) valvoc =
     &                    vdsum(icl,jcl,IXVOC) /  consum(icl,jcl,IXVOC)
            do 80 ispc=iptvoc,ipto3n-1
               vdepsa(icl,jcl,ispc) = valvoc
   80       continue
c
c   --- load O3 vdeps from regular model vdeps ----
c
            do 90 ispc=ipto3n,ipttim-1
               vdepsa(icl,jcl,ispc) = vdep(icl,jcl,idxo3)
   90       continue
c
c  ---- next cell ---
c
   60    continue
   50 continue
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
