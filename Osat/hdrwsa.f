c*** HDRWSA
c
      subroutine hdrwsa(iounit,fname,ftype,nspcs,nzcell,
     &                                           idate,btim,jdate,etim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine writes the header information to the output tracer
c     instantaneous files and tracer surface concentration files.
c     The data is written to the unit number in the argument list.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c    Argument description:
c      iounit  I  unit number of file to write
c      fname   C  name of file being written
c      ftype   C  file type of file being written 
c      nspcs   I  number of species to write to file
c      nzcell  I  number of layers (will be 1 for surface file)
c      idate   I   beginning date of simulation (YYJJJ)
c      btim    R   beginning time of simulation
c      jdate   I   ending date of simulation (YYJJJ)
c      etim    R   ending time of simulation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c      1/20/99   Grid cell size on file should be meters for all
c                cartesian projections (UTM, LCP, PSP)
c      10/24/01  Removed BSWAP and converted integer strings to character*4
c      11/06/01  Input dates are now Julian
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'grid.com'
      include 'flags.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer      iounit
      character*80 fname
      character*10 ftype
      integer      nspcs
      integer      nzcell       
      integer      idate
      real         btim
      integer      jdate
      real         etim
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*4  ifile(10), note(60), ispec(10,MXTRSP)
      integer      i, j, nseg, ndate, ndlast, izero, nspcout
      real         ttime, ttlast, zero, factr
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data zero /0.0/
      data izero /0/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
      ndate = idate
      ttime = btim/100.0
      ndlast = jdate
      ttlast = etim/100.0
c
c   --- set scaling factor for coordinates ---
c
      if( .NOT.llatlon ) then
         factr = 1000.0
      else
         factr = 1.0
      endif
c
c   ---- put species names into integer array ----
c 
      nspcout = 0
      do 10 j=1,nspcs
          if( loutsa(j) .OR. ftype .EQ. 'AIRQUALITY' ) then
              nspcout = nspcout + 1
              read(ptname(j),'(10A1)') (ispec(i,nspcout),i=1,10)
          endif
   10 continue
c
c   --- file description header ---
c
      read(ftype,'(10A1)') (ifile(i),i=1,10)
      read(runmsg(1:60),'(60A1)') (note(i),i=1,60)
      nseg = 1
      write(iounit,ERR=7000) ifile, note, nseg, nspcout, ndate, ttime,
     &                                                    ndlast, ttlast
c
c   --- regeion description header ----
c
      write(iounit,ERR=7000) zero, zero, iuzon, xorg*factr, yorg*factr, 
     &                       delx*factr, dely*factr, ncol(1), nrow(1), 
     &                           nzcell, izero, izero, zero, zero, zero
c
c  --- segment description header ---
c
      write(iounit,ERR=7000) izero, izero, ncol(1), nrow(1)
c
c  ---- species list ---
c
      write(iounit,ERR=7000) ((ispec(i,j),i=1,10),j=1,nspcout)
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in HDRWSA:'
      write(iout,9000,ERR=9999)'Writing output tracer file: ',
     &                                            fname(:istrln(fname))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,1X,2A)
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
c
