      subroutine cvtdate(nyr,nmo,ndy,time,ndate,jdate,ly2k) 
c   
c-----CAMx v4.02 030709
c   
c     CVTDATE converts year, month, day to calender (YYMMDD) and
c     julian (YYJJJ) dates
c                             
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c             
c     Modifications:   
c        11/05/01  Added Y2K flag and conversion to Julian date
c 
c     Input arguments:
c        nyr                 year (YYYY)
c        nmo                 month (MM)
c        ndy                 day (DD)
c        time                time (HHMM)
c 
c     Output arguments: 
c        time                time (HHMM)
c        ndate               calender date (YYMMDD)
c        jdate               Julian date (YYJJJ)
c        ly2k                Y2K flag
c               
c     Routines Called:   
c        JULDATE
c               
c     Called by:   
c        STARTUP 
c 
      logical ly2k
      dimension nday(12)
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----Entry point
c
      ly2k = .false.
      if (nyr.ge.2000) ly2k = .true.
      nday(2) = 28
      if (mod(nyr,4).eq.0) then
        nday(2) = 29
      endif
c
c-----Convert hour 2400 to 0000 of the following day
c
      if (time.ge.2400.) then 
        time = time - 2400. 
        ndy = ndy + 1 
        if (ndy.gt.nday(nmo)) then 
          ndy = 1
          nmo = nmo + 1 
          if (nmo.gt.12) then
            nmo = 1
            nyr = nyr + 1
          endif
        endif
      endif
c
c-----Create calendar and Julian date stamps
c
      ndate = 10000*mod(nyr,100) + 100*nmo + ndy
      jdate = ndate
      call juldate(jdate)
c
      return
      end
