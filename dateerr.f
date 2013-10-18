      subroutine dateerr(ly2k,string)
c   
c-----CAMx v4.02 030709
c   
c     DATEERR writes an error message whenever date mismatches occur
c             when reading input files
c                             
c     Copyright 2001, 2002, 2003
c     ENVIRON International Corporation
c             
c     Modifications:   
c       02/24/3003  --gwilson--  Added string to identify file being read. 
c 
c     Input arguments:
c        ly2k             Year 2000 flag (T is >=2000)
c 
c     Output arguments: 
c        none 
c               
c     Routines Called:   
c        CAMXERR
c               
c     Called by:   
c        READINP
c 
      include "camx.prm"
      include "filunit.com"
c
      logical ly2k
      character*(*) string
c
c-----Entry point
c
      write(iout,*) 
      write(iout,*) string
      write(iout,'(2A)') 'A date error or mismatch has occurred ',
     &                                             'reading the file.'
      write(iout,*) 
      if (ly2k) then
        write(iout,*)'You are modeling year 2000 or later --'
        write(iout,*)'CAMx assumes that all input file dates',
     &               ' are in Julian format (YYJJJ)'
      endif
      write(iout,*) 
c
      call camxerr()
c
      stop
c
      end
