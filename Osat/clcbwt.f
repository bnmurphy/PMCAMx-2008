c**** CLCBWT
c
      subroutine clcbwt(idate,btim,jdate,etim,ncols,nrows,nlays)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine calculates the weighted reactivity factor for VOC
c   species for the boundary conditions.  The mass is weighted by layer
c   thickness giving the weighted average for the cell.  The average
c   over all cells is then calculated.  The averages are calculated for
c   the entire bounday and for each boundary seperately.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Argument declarations:
c        idate   I   beginning date of simulation (YYJJJ)
c        btim    R   beginning time of simulation
c        jdate   I   ending date of simulation (YYJJJ)
c        etim    R   ending time of simulation
c        ncols    I   number of columns in coarse grid
c        nrows    I   number of rows in coarse grid
c        nlays    I   number of layers in coarse grid
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/29/96   --gwilson--    Original development
c     11/06/01   --cemery--     Input dates are now Julian
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'chmstry.com'
      include 'bndary.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
      integer   jdate
      real      etim
      integer   ncols
      integer   nrows
      integer   nlays
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 tpspc
      integer      ndate, jdlast, i, j
      integer      idtnow
      real         ttime, ttlast, timnow
      real         sumall, sumvoc 
      real         consum(MXSPEC,0:IDXBTP), ctin, contop(MXSPEC)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- call routine to get the file pointers to the proper place ---
c
      call bndprep(btim,idate,etim,jdate)
      do 15 i=1,MXSPEC
         contop(i) = bdnl(i)
   15 continue
c
c   --- read the top concentration file and store in species array ---
c
      rewind(itopc)
  111 continue
      read(itopc,8000,ERR=7000,END=222) tpspc, ctin
      do 10 i=1,nspec
        if( tpspc .eq. spname(i) ) then
          contop(i) = ctin
          goto 111
        endif
   10 continue
      goto 111
c
c   --- set ending date and time to be consistent with how time is 
c       counted here ----
c      
  222 continue
      ndate = idate
      ttime = btim/100.0
      jdlast = jdate
      ttlast = etim/100.0
      if( ttlast .EQ. 0. ) then
         jdlast = jdlast - 1
         ttlast = 24.0
      endif
c
c  --- initialize the dates and times ---
c
      idtnow = ndate
      timnow = ttime
c
      call rdsumbc(idtnow,timnow,jdlast,ttlast,ncols,nrows,nlays,
     &             contop,consum)
c
c   --- calculate the fractions ---
c
      if( lbndry ) then
         do 91 j=1,IDXBTP
            sumall = 0.
            sumvoc = 0.
            do 12 i=1,nspec
               if( lvocsp(i) .AND. consum(i,j) .GT. 0. ) then
                   sumall = sumall + consum(i,j)
                   sumvoc = sumvoc + consum(i,j) * reacrt(i) 
               endif
   12       continue
            if( sumall .GT. 0. ) then
                vocwt(iptvoc+j) = sumvoc / sumall
            else
                vocwt(iptvoc+j) = 0.
            endif
   91    continue
      else
         sumall = 0.
         sumvoc = 0.
         do 22 i=1,nspec
            if( lvocsp(i) .AND. consum(i,0) .GT. 0. ) then
                sumall = sumall + consum(i,0)
                sumvoc = sumvoc + consum(i,0) * reacrt(i) 
            endif
   22    continue
         if( sumall .GT. 0. ) then
            vocwt(iptvoc+1) = sumvoc / sumall
         else
            vocwt(iptvoc+1) = 0.
         endif
      endif
c
c  --- fill the global arrays of top concentrations for tracer species ---
c
      if( lbndry ) then
         ioff = IDXBTP
      else
         ioff = 1
      endif
      do 32 i=1,nspec
         if( lvocsp(i) ) then
            ptloft(iptvoc+ioff) = ptloft(iptvoc+ioff) + 
     &                             contop(i) * crbnum(i)
         else if( lnoxsp(i) ) then
            ptloft(iptnox+ioff) = ptloft(iptnox+ioff) + contop(i)
         else if( lo3sp(i) ) then
            ptloft(ipto3n+ioff) = ptloft(ipto3n+ioff) + contop(i)/2.0
            ptloft(ipto3v+ioff) = ptloft(ipto3v+ioff) + contop(i)/2.0
         endif
   32 continue
      do 42 i=1,ntotsp
         ptloft(i) = MAX( ptloft(i), BNDLPT )  
   42 continue
c
c  --- return to the calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in CLCBWT:'
      write(iout,'(/,1X,A)') 'Reading top concentrations file.'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 8000 format(A10,F10.0)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
