c**** CLCEWT
c
      subroutine clcewt(jdate,etim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine calculates the weighted reactivity factor for VOC
c   species in each of the source groups.  All of the emissions for the
c   group are read and the emissions are weighted by reactivity factor
c   and summed up.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Argument description:
c        jdate   I   ending date of simulation (YYJJJ)
c        etim    R   ending time of simulation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/04/96   --gwilson--    Original development
c     10/10/96   --gwilson--    Added code to output emissions and 
c                               reactivity for each grouping.
c     12/08/96   --gwilson--    Added code to set the index into tracer
c                               species list for the PiG sources
c     12/12/96   --gwilson--    Fixed bug in reporting total NOx Tons
c     01/09/97   --gwilson--    Fixed (another) bug in reporting total 
c                               NOx Tons
c     01/12/97   --gwilson--    Fixed bug in calculating the source
c                               region in the fine grid
c     11/06/01   --cemery--     Input dates are now Julian
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'filunit.com'
      include 'bndary.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   jdate
      real      etim
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c  CVTTON   R   conversion factor for grams to tons
c
      real   CVTTON
c
      parameter( CVTTON = 907184.7 )
c
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      idx
      integer      i, j, ivoc, ndlast
      integer      inox, ncount, ioff
      real         emssum(MXSPEC,MXTRSP)
      real         emsbas(MXSPEC,MXTRSP), emsoth(MXSPEC,MXTRSP)
      real         emslft(MXCOLA,MXROWA), emstot(MXCOLA,MXROWA)
      real         sumall, sumvoc, sumnox, difmax, diff, ttlast
      real         tonnox(MXTRSP), tonvoc(MXTRSP)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- set the date and times ---
c
      ndlast = jdate 
      ttlast = etim/100.0
      if( ttlast .EQ. 0.0 ) then
         ttlast = 24.0
         ndlast = ndlast - 1 
      endif
c
c   --- initialize the array to zero ---
c
      do 10 j=1,MXTRSP
        tonnox(j) = 0.
        tonvoc(j) = 0.
        do 15 i=1,MXSPEC
          emssum(i,j) = 0.
          emsbas(i,j) = 0.
          emsoth(i,j) = 0.
   15   continue
   10 continue
      do 16 j=1,MXROWA
        do 17 i=1,MXCOLA
           emstot(i,j) = 0.
           emslft(i,j) = 0.
   17   continue
   16 continue
c
c   --- loop over all of the groups ----
c
      call sumgrps(ndlast,ttlast,emstot,emslft,emsbas,emsoth,emssum)
c
c  --- check that all emissions are accounted for (could have some
c      machine fuzz) ----
c
      difmax = -99999999.
      do 80 j=1,nrow(1)
         do 90 i=1,ncol(1)
            if( emstot(i,j) .NE. 0. ) then 
               diff = ABS( (emstot(i,j) - emslft(i,j)) / emstot(i,j) )
               difmax = MAX( difmax,diff )
            endif
   90     continue
   80 continue
      if( .NOT. leftovr .AND. difmax .GE. 0.0001 ) goto 7003
      if( leftovr .AND. difmax .LE. 0.0001 ) goto 7004
c
c  --- all emissions are summed, calculate the weghted fraction ----
c
      do 11 i=1,nsaspc
c
c   --- ignore if this is an initial condition or boundary condition
c       tracer ---
c
         if( ptname(i)(7:8) .EQ. 'IC' ) goto 11
         if( ptname(i)(7:8) .EQ. 'BC' ) goto 11
         sumall = 0.
         sumvoc = 0.
         do 21 idx=1,nspec
            if( lvocsp(idx) .AND. emssum(idx,i) .GT. 0 ) then
                sumall = sumall + emssum(idx,i)
                sumvoc = sumvoc + emssum(idx,i) * reacrt(idx) 
                tonvoc(i) = tonvoc(i) +
     &                     (emssum(idx,i)*crbnum(idx)*16.0)/CVTTON
            endif
            if( lnoxsp(idx) .AND. emssum(idx,i) .GT. 0 ) then
                tonnox(i) = tonnox(i) + (emssum(idx,i)*46.0)/CVTTON
            endif
  21     continue
         if( sumall .GT. 0. ) then
             vocwt(i) = sumvoc / sumall
         else
             vocwt(i) = 0.
         endif
  11  continue
c
c  --- calculate the "leftover" group from lump sums ---
c
      if( leftovr ) then
         do 31 i=1,nregin
            sumall = 0.
            sumvoc = 0.
            do 41 idx=1,nspec
               ivoc = iemvoc - 1 + i + ngroup*nregin
               diff = emsbas(idx,ivoc) - emsoth(idx,ivoc)
               if( lvocsp(idx) .AND. diff .GT. 0. ) then
                   sumall = sumall + diff
                   sumvoc = sumvoc + diff * reacrt(idx)
                   tonvoc(ivoc) = tonvoc(ivoc) + 
     &                               (diff*crbnum(idx)*16.0)/CVTTON
               endif
               inox = iemnox - 1 + i + ngroup*nregin
               diff = emsbas(idx,inox) - emsoth(idx,inox)
               if( lnoxsp(idx) .AND. diff .GT. 0. ) then
                   inox = iemnox - 1 + i + ngroup*nregin
                   tonnox(inox) = tonnox(inox) + (diff*46.0)/CVTTON
               endif
   41        continue
             if( sumall .GT. 0 ) then
                 vocwt(ivoc) = sumvoc / sumall
             else
                  vocwt(ivoc) = 0.
             endif
   31    continue
      endif
c
c  --- echo the data ---
c
      if( ngroup .EQ. 0 ) then
          ioff = 0
          ncount = 0
      else
          if( leftovr ) then
              ncount = ngroup + 1
          else
              ncount = ngroup
          endif
          ioff = 1
      endif
      do 51 i=ioff,ncount
         write(idiag,9004) 'Species   ','     Average ',
     &                                              'Emissions (tons)' 
         write(idiag,9005) ' ','Reactivity','    VOC   ','    NOx   '
         write(idiag,9006) ('-',j=1,60)
         sumnox = 0.
         do 61 j=1,nregin
            if( i .GT. 0 ) then
               inox = iemnox - 1 + j + (i-1)*nregin
            else
               inox = iemnox - 1 + j
            endif
            write(idiag,9007) ptname(inox),0.,0.,tonnox(inox)
            sumnox = sumnox + tonnox(inox)
   61    continue
         write(idiag,9006) ('-',j=1,60)
         write(idiag,9009) 'Total',0.,0.,sumnox
         write(idiag,9004) 'Species   ','     Average ',
     &                                              'Emissions (tons)' 
         write(idiag,9005) ' ','Reactivity','    VOC   ','    NOx   '
         write(idiag,9006) ('-',j=1,60)
         sumvoc = 0.
         do 71 j=1,nregin
            if( i .GT. 0 ) then
               ivoc = iemvoc - 1 + j + (i-1)*nregin
            else
               ivoc = iemvoc - 1 + j
            endif
            write(idiag,9008) ptname(ivoc),vocwt(ivoc),tonvoc(ivoc),0.
            sumvoc = sumvoc + tonvoc(ivoc)
   71    continue
         write(idiag,9006) ('-',j=1,60)
         write(idiag,9010) 'Total',0.,sumvoc,0.
   51 continue
c
c  --- return to the calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in CLCEWT:'
      write(iout,'(/,1X,2A)') 'There is a significant amount ',
     &         'of emissions unaccounted for in source groupings.'
      write(iout,'(1X,2A)') 'You should turn on the "leftover group"',
     &                                          ' flag in job script.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,a)') 'ERROR in CLCEWT:'
      write(iout,'(/,1X,2A)') 'The "leftover" emissions group ',
     &                       'has an insignificant amount of emissions.'
      write(iout,'(1X,2A)') 'You should turn off the ',
     &                            '"leftover group" flag in job script.'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9004 format(/,1X,A10,3X,A13,8X,A)
 9005 format(1X,A10,6X,A10,3X,A10,3X,A10,3X,A10)
 9006 format(100(A1))
 9007 format(1X,A10,6X,F10.1,3X,F10.1,3X,F10.4)
 9008 format(1X,A10,6X,F10.3,3X,F10.4,3X,F10.1)
 9009 format(1X,A10,6X,F10.1,3X,F10.1,3X,F10.2)
 9010 format(1X,A10,6X,F10.3,3X,F10.2,3X,F10.1)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
