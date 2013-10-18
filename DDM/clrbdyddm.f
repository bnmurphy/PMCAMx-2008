c**** CLRBDYDDM.F
c
      subroutine clrbdyddm(nox,noy,noz,nspsa,saconc)
c
c-----CAMx 4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine resets the boundary conditions for the DDM species to 
c   zero.  
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c      Argument description:
c       Outputs:
c           saconc   R  array of concentrations to clear
c       Inputs:
c           nox      I  number of cells in X direction
c           noy      I  number of cells in Y direction
c           noz      I  number of vertical cells in grid
c           nspsa    I  number of species 
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
      integer   nspsa
      real      saconc(nox,noy,noz,nspsa)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer nedge, ioff, icl, jcl, izcl, iddm, i
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- set the number of BC edges --
c
      if( lbndry ) then
        nedge = 5
      else
        nedge = 1
      endif
c
c  --- do the WEST and EAST boundaries ----
c
      do 10 jcl=2,noy-1
c
c  --- West boundary ---
c
         if( ibeg(jcl) .EQ. -999 ) goto 10
         icl = ibeg(jcl) - 1
         do 20 izcl=1,noz
c
c   --- if stratifying by boundary, put in seperate position ---
c
            if( lbndry ) then
               ioff = IDXBWS
            else
               ioff = 1
            endif
c
c   --- reset boundary value for each modeled species ---
c
            do i=1,nspec        
              do j=1,nbcddm
                iddm = iptddm(i) + (j-1)*nedge + nicddm + ioff - 1
                saconc(icl,jcl,izcl,iddm) = 0.
              enddo
            enddo
   20    continue
c
c  --- East boundary ---
c
         icl = iend(jcl) + 1
         do 30 izcl=1,noz
c
c   --- if stratifying by boundary, put in seperate position ---
c
            if( lbndry ) then
               ioff = IDXBES
            else
               ioff = 1
            endif
c
c   --- reset boundary value for each modeled species ---
c
            do i=1,nspec        
              do j=1,nbcddm
                iddm = iptddm(i) + (j-1)*nedge + nicddm + ioff - 1
                saconc(icl,jcl,izcl,iddm) = 0.
              enddo
            enddo
   30    continue
   10 continue
c
c  --- do the SOUTH and NORTH boundaries ----
c
      do 40 icl=2,nox-1
c
c  --- South boundary ---
c
         if( jbeg(icl) .EQ. -999 ) goto 40
         jcl = jbeg(icl) - 1
         do 50 izcl=1,noz
c
c   --- stratifying by boundary, put in seperate position ---
c
            if( lbndry ) then
                ioff = IDXBST
            else
                ioff = 1
            endif
c
c   --- reset boundary value for each modeled species ---
c
            do i=1,nspec        
              do j=1,nbcddm
                iddm = iptddm(i) + (j-1)*nedge + nicddm + ioff - 1
                saconc(icl,jcl,izcl,iddm) = 0.
              enddo
            enddo
   50    continue
c
c  --- North boundary ---
c
         jcl = jend(icl) + 1
         do 60 izcl=1,noz
c
c   --- if stratifying by boundary, put in seperate position ---
c
            if( lbndry ) then
                ioff = IDXBNT
            else
                ioff = 1
            endif
c
c   --- reset boundary value for each modeled species ---
c
            do i=1,nspec        
              do j=1,nbcddm
                iddm = iptddm(i) + (j-1)*nedge + nicddm + ioff - 1
                saconc(icl,jcl,izcl,iddm) = 0.
              enddo
            enddo
   60    continue
   40 continue
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
