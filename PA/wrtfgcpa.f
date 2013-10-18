      subroutine wrtfgcpa(iendat,endtim)
c
c-----CAMx v4.02 030709
c
c     This routine writes the outptut file for the process analysis
c     gridded reaction rate data (CPA data).  This routine is writes 
c     the data for the fine grid.  The format is the same as the 
c     coarse grid average concentration file for the regular model 
c     species.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation 
c
c     Modifications:
c        none
c
c     Input arguments:
c        iendat              current ending date (YYJJJ)
c        endtim              current ending time (HHMM)
c
c     Output arguments:
c        none
c
c     Rountines Called:
c        none
c
c     Called by:
c        CAMx
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'flags.com'
      include 'filunit.com'
      include 'camxfld.com'
      include 'chmstry.com'
      include 'grid.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   iendat
      real      endtim
      integer   ifgptr(MXGRID), ifglvl(MXGRID)
      integer   igrd, idx, i, j, k, l, nlayer(MXGRID)
      real      cnctmp(MXCOLA,MXROWA)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c-----Set the fine grid pointers
c
      do i = 1,ngrid
        ifglvl(i) = 0
      enddo
      do i = 1,ngrid
        do j = 1,nchdrn(i)
          idch = idchdrn(j,i)
          ifgptr(idch) = i - 1 
          ifglvl(idch) = ifglvl(idch) + 1
          do k = 1,nchdrn(idch)
            ifglvl(idchdrn(k,idch)) = ifglvl(idchdrn(k,idch)) + 1
          enddo
        enddo
      enddo
c
c   --- set the number of layers for each grid ---
c
      do igrd = 1,ngrid
        if (l3davg) then
          nlayer(igrd) = nlay(igrd)
        else
          nlayer(igrd) = 1
        endif
      enddo
c
c  --- if this is the first time through, the write the header ---
c      
      if( lfirst ) then
          write(IOWSFC+IDXFIN,ERR=7000) runmsg
          write(IOWSFC+IDXFIN,ERR=7000) ngrid-1, ntotsp
          write(IOWSFC+IDXFIN,ERR=7000) (ptname(i),i=1,ntotsp)
          do igrd = 2,ngrid
             write(IOWSFC+IDXFIN,ERR=7000) inst1(igrd), jnst1(igrd),   
     &        inst2(igrd), jnst2(igrd), meshold(igrd), meshold(igrd),
     &                            ncol(igrd), nrow(igrd), nlayer(igrd), 
     &                                        ifgptr(igrd), ifglvl(igrd)
          enddo
          lfirst = .FALSE.
      endif
c
c-----Write time span
c
      write(IOWSFC+IDXFIN) endtim, iendat
c
c-----Write the concentration data for this hour
c
      do 30 igrd = 2,ngrid
c
        do 40 l = 1,ntotsp
          do 50 k = 1,nlayer(igrd)
            do 60 j = 1,nrow(igrd)
              do 70 i = 1,ncol(igrd)
                idx = i + ncol(igrd)*(j-1) + 
     &                ncol(igrd)*nrow(igrd)*(k-1) +
     &                ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                cnctmp(i,j) = ptconc(ipsa3d(igrd)-1+idx)
   70         continue
   60       continue
            write(IOWSFC+IDXFIN) 
     &             ((cnctmp(i,j),i=1,ncol(igrd)),j=1,nrow(igrd))
   50     continue
   40   continue
c
   30 continue
c
c-----Close files and return to calling routine
c
      call flush(IOWSFC+IDXFIN)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in WRTFGCPA'
      write(iout,9000,ERR=9999) 'Writing gridded output chemical ',        
     &                'process analysis (.cpa) file for coarse grid.'  
      call camxerr()
c     
c-----------------------------------------------------------------------
c    Format statements:   
c-----------------------------------------------------------------------
c      
 9000 format(/,1X,3A)
c      
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
