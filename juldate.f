      subroutine juldate(idate)
c 
c-----CAMx v4.02 030709
c 
c     JULDATE converts date from calender (YYMMDD) format to Julian
c     (YYJJJ) format
c                           
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none
c 
c     Input arguments: 
c        idate               calender date (YYMMDD) 
c             
c     Output arguments: 
c        idate               julian date (YYJJJ) 
c             
c     Routines Called: 
c        none 
c             
c     Called by: 
c        CNCPREP
c        RDFGCON
c        READAHO
c        READINP
c        STARTUP 
c
      dimension nday(12)
      data nday/31,28,31,30,31,30,31,31,30,31,30,31/
c
c-----Entry point
c
      iyear = idate/10000
      imonth = (idate - iyear*10000)/100
      iday = idate - iyear*10000 - imonth*100
c
      nday(2) = 28
      if (mod(iyear,4).eq.0) nday(2) = 29
      mday = 0
      do 10 n = 1,imonth-1
        mday = mday + nday(n)
 10   continue
      jday = mday + iday
      idate = iyear*1000 + jday
c
      return
      end
