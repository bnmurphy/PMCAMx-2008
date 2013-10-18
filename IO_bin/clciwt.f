c**** CLCIWT
c
      subroutine clciwt(idate,btim,jdate,etim,ly2k,ncols,nrows,nlays)
c
c-----CAMx v4.03 031205
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine calculates the weighted reactivity factor for VOC
c   species for the initial conditions.  The mass is weighted by layer
c   thickness giving the weighted average for the cell.  The average
c   over all cells is then calculated.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation     
c
c      Argument description:
c        idate   I   beginning date of simulation (YYJJJ)
c        btim    R   beginning time of simulation
c        jdate   I   ending date of simulation (YYJJJ)
c        etim    R   ending time of simulation
c        ly2k    L   Year 2000 flag (T is >=2000)
c        ncols   I   number of columns in coarse grid
c        nrows   I   number of rows in coarse grid
c        nlays   I   number of layers in coarse grid
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/04/96   --gwilson--    Original development
c     11/06/01   --cemery--     Input dates are now Julian
c     11/20/03   --gwilson--    Fixed bug when using the species map
c                               to get the index of modeled species
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'filunit.com'
      include 'chmstry.com'
      include 'bndary.com'
      include 'grid.com'
      include 'flags.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
      integer   jdate
      real      etim
      integer   ncols
      integer   nrows
      integer   nlays
      logical   ly2k
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*4  iname(10)
      integer      idx, iseg, ispc, izcl, i, j, k
      integer      ibegic, iendic, ibgdhp, idtnow, iunit
      real         btimic, etimic, btimhp, timnow
      real         sumall, sumvoc
      real         consum(MXSPEC), congrd(MXCOLA*MXROWA*MXLAYA)
      real         tmp3d(MXCOLA,MXROWA,MXLAYA)
      real         height(MXCOLA,MXROWA,MXLAYA)
      real         depth(MXCOLA,MXROWA,MXLAYA)
      logical      lfound
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- call routine to get the file pointers to the proper place ---
c
      if( .NOT. lrstrt ) then
         iunit = iic
      else
         iunit = irstc
      endif
      ibgdhp = 0
      btimhp = 0.0
c
c   --- call routine to get the file pointers to the proper place,
c       then rewind the files so regular model can write headers again ---
c
      rewind(iunit)
      call cncprep(etim,jdate)
c
c  --- set the current date and time ----
c
      idtnow = idate
      timnow = btim/100.0
      call getdepth(ncols,nrows,nlays,ly2k,ibgdhp,idtnow,btimhp,timnow,
     &              ihtp(1),height,depth)
c
c   --- initialize the array to zero ---
c
      do 50 i=1,MXSPEC
        consum(i) = 0.
   50 continue
c
c   --- read the initial conditions ----
c
  111 continue
      read(iunit,ERR=7001,END=7002) ibegic, btimic, iendic, etimic
      lfound = .TRUE.
      if( ibegic .EQ. idtnow .AND. btimic .LT. timnow) lfound = .FALSE.
      if( ibegic .GT. idtnow ) lfound = .FALSE.
      do 60 ispc=1,nicspc
          do 70 izcl = 1,nlays
              read(iunit,ERR=7001,END=7002) iseg, (iname(i),i=1,10), 
     &                    ((tmp3d(i,j,izcl),i=1,ncols),j=1,nrows)
   70     continue
c
          do k=1,nlays
            do j=1,nrows
              do i=1,ncols
                n3d = i + (j-1)*ncols + (k-1)*ncols*nrows
                congrd(n3d) = tmp3d(i,j,k)
              enddo
            enddo
          enddo
c
c   --- if record does not span this hour, skip it ----
c
         if( .NOT. lfound ) goto 60
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
c
c   --- if the species is a not modelled or not a VOC species skip it ---
c
          idx = 0
          do i=1,nspec
            if( licmap(i,1) .EQ. ispc ) idx = i
          enddo
          if( idx .LE. 0 ) goto 60
          if( .NOT. lvocsp(idx) ) goto 60
c
          call sumicwt(ncols,nrows,nlays,ibeg,iend,deltax,depth,congrd,
     &                 consum(idx))
c

c
c   --- next species ---
c
   60 continue
c
c  --- if hour not found, go back and read it ---
c
      if( .NOT. lfound ) goto 111
c
c  --- all concentrations are summed, calculate the weghted fraction ----
c
      sumall = 0.
      sumvoc = 0.
      do 21 i=1,nspec
         if( lvocsp(i) .AND. consum(i) .GT. 0. ) then
             sumall = sumall + consum(i)
             sumvoc = sumvoc + consum(i) * reacrt(i) 
         endif
  21  continue
      if( sumall .GT. 0. ) then
          vocwt(iptvoc) = sumvoc / sumall
      else
          vocwt(iptvoc) = 0.
      endif
c
c  --- rewind the files to be used by the regular model ----
c
      rewind(iunit)
      rewind(ihtp(1))
c
c  --- return to the calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in CLCIWT:'
      write(iout,'(/,1X,A)') 'Reading initial conditions.'
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in CLCIWT:'
      write(iout,'(/,1X,2A)') 'ERROR: Premature end-of-file reading ',
     &                                            'initial conditions.'
      call camxerr()
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
