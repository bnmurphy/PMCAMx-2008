c*** HDRCPRT
c
      subroutine hdrcprt(idate,btim,jdate,etim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine writes the header information for the receptor
c     decay rates file for RTRAC.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Argument declarations:
c        idate   I   beginning date of simulation (YYJJJ)
c        btim    R   beginning time of simulation
c        jdate   I   ending date of simulation (YYJJJ)
c        etim    R   ending time of simulation
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c    04/14/97  --gwilson--  Changed the way number of source groups
c                           is specified
c    11/16/99  --gwilson--  Fixed a bug in writing the header to the 
c                           receptor file.  There was an extra line being
c                           written that messed up the browser.
c    11/06/01  --cemery--   Input dates are now Julian
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'tracer.com'
      include 'rtracchm.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
      integer   jdate
      real      etim
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*40 cdate
      integer      ndate, ndlast, i
      real         ttime, ttlast
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- if no receptors requested, just bail ---
c 
      if( nrcprt .LE. 0 ) goto 9999
c
c  --- set up the date and time ---
c
      ndate = idate
      ttime = btim/100.0
      ndlast = jdate
      ttlast = etim/100.0
      if( ttlast .EQ. 0. ) then
          ttlast = 24.0
          ndlast = ndlast - 1
      endif
c
c  --- get the date and time and write header records ---
c
      call fdate( cdate )
      write(IOWAVG,9010,ERR=7000) 'CAMx', 
     &                      runmsg(:istrln(runmsg)), verson
      write(IOWAVG,'(A)',ERR=7000) cdate(1:24)
c
c  --- write the dates and times of simuation ---
c
      write(IOWAVG,9000,ERR=7000)
      write(IOWAVG,9001,ERR=7000) ndate, ttime, ndlast, ttlast
      write(IOWAVG,9008,ERR=7000)'Average Interval ',ANINT(dtout)/60.
c
      write(IOWAVG,9002,ERR=7000) 'Number of gas tracers groupings ',
     &                                                          nrtgas
c
c   ---- write the tracer species names ----
c
      write(IOWAVG,9000,ERR=7000)
      write(IOWAVG,9000,ERR=7000) 'Tracer Names'
      write(IOWAVG,9003,ERR=7000) ' ',' ',' ',(ptname(i),i=1,nrtgas)
c
c   --- do the receptors ---
c
      write(IOWAVG,9002,ERR=7000) 'Number of receptors ',nrcprt
      write(IOWAVG,9007,ERR=7000) 'Receptor','Grid','I-cell',
     &                                          'J-cell','K-cell'
      do i=1,nrcprt
        write(IOWAVG,9006,ERR=7000) i, idomrt(i), ircprt(i), 
     &                                     jrcprt(i), krcprt(i)
      enddo
c
c  --- write the header for rest of data ----
c
      write(IOWAVG,9000,ERR=7000)
      write(IOWAVG,9000,ERR=7000) 'Time Varying Tracer Data'
c
c  --- return to calling routine ---
c
      call flush(IOWAVG)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in HDRRCPRT:'
      write(iout,9020,ERR=9999) 'Writing Header of decaying rates ',
     &                                   'file for receptors in RTRAC.'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(2(:,A,','))
 9001 format(' File Duration  ,',2(I10.5,',',F10.2,','))
 9002 format(A,',',I10)
 9003 format(500(:,A9,','))
 9006 format(5(I10,','))
 9007 format(5(A10,','))
 9008 format(A,',',F10.4)
 9010 format(3(A,','))
 9020 format(/,1X,2(A,','))
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
