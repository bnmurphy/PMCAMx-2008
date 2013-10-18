c*** WRFGSA 11/22/96
c
      subroutine wrfgsa(iendat,endtim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine writes the instantaneous file for the tracer
c     species for the fine grids.  
c     If the hour is odd then the first file is written,
c     else the second file is written.  The file is rewound and the
c     header rewritten as well as the hourly data.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c    Argument description:
c     Outputs:
c     Inputs:
c        iendat  I   current date
c        endtim  R   current time
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'grid.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   iendat
      real      endtim
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  ptout(MXTRSP)
      integer       ifgptr(MXGRID), ifglvl(MXGRID)
      integer       iounit, igrd, idx, i, j, k, l, nspcout
      real          cnctmp(MXCOLA,MXROWA), avgtmp(MXCOLA,MXROWA)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- figure out which unit number and file to use --- 
c 
      if( MOD( INT(endtim/100.), 2 ) .EQ. 1 ) then
          iounit = IOWCN1+IDXFIN
          fname = cn1fil(IDXFIN)
      else
          iounit = IOWCN2+IDXFIN
          fname = cn2fil(IDXFIN)
      endif
c
c   --- set the find grid pointers ---
c
      do i=1,ngrid
         ifglvl(i) = 0
      enddo
      do i=1,ngrid
         do j=1,nchdrn(i)
            idch = idchdrn(j,i)
            ifgptr(idch) = i - 1 
            ifglvl(idch) = ifglvl(idch) + 1
            do k=1,nchdrn(idch)
               ifglvl(idchdrn(k,idch)) = ifglvl(idchdrn(k,idch)) + 1
            enddo
         enddo
      enddo
c
c  --- load the species names into local array ---
c
      nspcout = 0
      do i=1,ntotsp
         if( loutsa(i) ) then
             nspcout = nspcout + 1
             ptout(nspcout) = ptname(i)
         endif
      enddo
c
c   --- rewind file and write the file description header ---
c
      rewind(iounit)
      write(iounit,ERR=7000) runmsg
      write(iounit,ERR=7000) ngrid-1, ntotsp
      if( lfirst ) then
          write(IOWSFC+IDXFIN,ERR=7000) runmsg
          write(IOWSFC+IDXFIN,ERR=7000) ngrid-1, nspcout
      endif
c
c  ---- species list ---
c
      write(iounit,ERR=7000) (ptname(i),i=1,ntotsp)
      if( lfirst ) then
          write(IOWSFC+IDXFIN,ERR=7000) (ptout(i),i=1,nspcout)
      endif
c
c  --- read the fine grid descriptions, and check for 
c      consistancy ---
c
      do 20 igrd=2,ngrid
          write(iounit,ERR=7000) inst1(igrd), jnst1(igrd), 
     &                                     inst2(igrd), jnst2(igrd), 
     &             meshold(igrd), meshold(igrd), ncol(igrd), nrow(igrd), 
     &                            nlay(igrd), ifgptr(igrd), ifglvl(igrd)
          if( lfirst ) then
             write(IOWSFC+IDXFIN,ERR=7000) inst1(igrd), jnst1(igrd),   
     &           inst2(igrd), jnst2(igrd), meshold(igrd), meshold(igrd),
     &             ncol(igrd), nrow(igrd), 1, ifgptr(igrd), ifglvl(igrd)
          endif
   20 continue
c
c  --- make sure time span is correct ---
c
      write(iounit,ERR=7000) endtim, iendat
      write(IOWSFC+IDXFIN,ERR=7000) endtim, iendat
c
c   --- read the data for this hour ----
c
      do 50 igrd=2,ngrid
         do 60 l=1,ntotsp
            do 70 k=1,nlay(igrd)
               do 80 j=1,nrow(igrd)
                  do 90 i=1,ncol(igrd)
                     idx =  i + ncol(igrd)*(j-1) + 
     &                      ncol(igrd)*nrow(igrd)*(k-1) +
     &                           ncol(igrd)*nrow(igrd)*nlay(igrd)*(l-1)
                     cnctmp(i,j) = ptconc(ipsa3d(igrd)-1+idx)
                     idx =  i + ncol(igrd)*(j-1) + 
     &                           ncol(igrd)*nrow(igrd)*(l-1)
                     if( loutsa(l) ) avgtmp(i,j) = 
     &                                       ptavrg(ipsa2d(igrd)-1+idx)
   90             continue
   80          continue
               write(iounit) ((cnctmp(i,j),i=1,ncol(igrd)),
     &                                               j=1,nrow(igrd))
   70       continue
            if( loutsa(l) ) write(IOWSFC+IDXFIN) 
     &                   ((avgtmp(i,j),i=1,ncol(igrd)),j=1,nrow(igrd))
   60    continue
   50 continue
c
c  --- close file and return to calling routine ---
c
      lfirst = .FALSE.
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in WRFGSA:'
      write(iout,9000,ERR=9999)'Writing output tracer file: ',
     &                                            fname(:istrln(fname))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,1X,6A)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
