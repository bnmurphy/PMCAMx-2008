c*** INSTSA
c
      subroutine instsa(idum,iendat, endtim, nox, noy, noz, nspsa, 
     &                                                 saconc, saavrg)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine writes the instantaneous file for the tracer
c     species.  If the hour is odd then the first file is written,
c     else the second file is written.  The file is rewound and the
c     header rewritten as well as the hourly data.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c    Argument description:
c      idum    I  not used -- needed for NETCDF version
c      iendat  I  ending date for this hour
c      endtim  R  ending time for this hour
c      nox     I  number of X cells in domain
c      noy     I  number of Y cells in domain
c      noz     I  number of layers in domain
c      nspsa   I  number of species
c      saconc  R  array of concentrations (instantaneous)
c      saavrg  R  array of surface concentrations (average)
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c      10/24/01  Removed BSWAP and converted integer strings to character*4
c      11/06/01  Removed calls to JULDATE and CALDATE (all dates in Julian)
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
       integer   idum
       integer   iendat
       real      endtim
       integer   nox
       integer   noy
       integer   noz
       integer   nspsa
       real      saconc(nox,noy,noz,nspsa)
       real      saavrg(nox,noy,nspsa)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      integer       idate, jdate
      integer       iounit, isegmt, i, j, k, l, nspcout
      character*4   ispec1(10,MXTRSP), ispec2(10,MXTRSP)
      real          btim, etim
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- set up the dates and times for instantaneous files ---
c
      jdate = iendat
      btim = endtim 
      etim = endtim + 10.0
      idate = jdate
      if( etim .EQ. 24.0 ) then
         jdate = jdate + 1
         etim = 0.0
         if( MOD(jdate,1000) .GT. 365 ) then
            if( MOD(INT(jdate/1000),4) .EQ. 0 ) then
               if( MOD(jdate,1000) .EQ. 367 )
     &                     jdate = (INT(jdate/1000)+1)*1000 + 1
            else
               jdate = (INT(jdate/1000)+1)*1000 + 1
            endif
         endif
      endif
c
c   ---- put species names into integer array ----
c 
      nspcout = 0
      do 10 j=1,ntotsp
          read(ptname(j),'(10A1)') (ispec1(i,j),i=1,10)
          if( loutsa(j) ) then
             nspcout = nspcout + 1
             read(ptname(j),'(10A1)') (ispec2(i,j),i=1,10)
          endif
   10 continue
c
c  --- figure out which unit number and file to use --- 
c 
      if( MOD( INT(btim/100.0), 2 ) .EQ. 1 ) then
          iounit = IOWCN1+IDXCRS
          fname = cn1fil(IDXCRS)
      else
          iounit = IOWCN2+IDXCRS
          fname = cn2fil(IDXCRS)
      endif
c
c  --- rewind file and call routine to write the header ---
c
      rewind(iounit)
      call hdrwsa( iounit, fname, 'AIRQUALITY', nsaspc, noz,
     &                                        idate, btim, jdate, etim )
c
c   --- write the data for this hour ----
c
      isegmt = 1
      btim = btim / 100.
      etim = etim / 100.
      write(iounit,ERR=7000) idate, btim, jdate, etim 
      do 20 l=1,nsaspc
         do 30 k=1,noz
            write(iounit) isegmt, (ispec1(i,l),i=1,10), 
     &                              ((saconc(i,j,k,l),i=1,nox),j=1,noy)
   30    continue
   20 continue
c
c  --- write the surface concentration file if requested ----
c
      if( lsfcfl(IDXCRS) ) then
c
c   --- set up the dates and times for instantaneous files ---
c
          jdate = iendat
          etim = ANINT(endtim)/100.
          idate = jdate
          if( dtout .GE. 60.0 ) then
              btim = ANINT( 1000*(etim - ANINT(dtout)/60.) )/1000.
          else
              btim = ANINT( 1000*(etim - ANINT(dtout)/100.) )/1000.
          endif
          if( btim .LT. 0. ) then
              btim = btim + 24.0
              idate = idate - 1
          endif
          write(IOWSFC+IDXCRS,ERR=7000) idate, btim, jdate, etim
          do 40 l=1,nsaspc
            if( loutsa(l) ) write(IOWSFC+IDXCRS) isegmt, 
     &           (ispec2(i,l),i=1,10),((saavrg(i,j,l),i=1,nox),j=1,noy)
   40     continue
          if( ntrtim .GT. 0 ) then
             do 50 l=nsaspc+1,ntotsp
               write(IOWSFC+IDXCRS) isegmt, (ispec2(i,l),i=1,10), 
     &                                        ((0.0,i=1,nox),j=1,noy)
   50        continue
          endif
      endif
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
      write(iout,'(//,A)') 'ERROR in INSTSA:'
      write(iout,9000,ERR=9999)'Writing output tracer file: ',
     &                                          fname(:istrln(fname))
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
