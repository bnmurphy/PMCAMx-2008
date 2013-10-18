      subroutine gresmscl(ngrid,time,date,idiag,pigdump)
c
c-----CAMx v4.02 030709
c
c     GRESMSCL calculates total NO and NO2 mass in all PiGs by grid
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        07/05/02  gwilson  Changed name to GRESMSCL to accommadate
c                           seperate routine for IRON-PIG
c
c     Input arguments:
c        ngrid               number of grids
c        time                model time (HHMM)
c        date                model date (YYJJJ)
c        idiag               Diagnostic file unit number
c        pigdump             PiG dumped mass of NO, NO2, HNO3 by grid (umol)
c
c     Output arguments:
c        none
c
c     Subroutines called:
c        none
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "pigsty.com"
c
      integer   date
      real*8 pigmass(MXSPEC,MXGRID),pigdump(MXSPEC,MXGRID)
c
c-----Entry point
c
      do igrd=1,ngrid
        npigon(igrd) = 0
        do l = 1,3
          pigmass(l,igrd) = 0.
        enddo
      enddo
c
      do n=1,npig
        igrd = ingrd(n)
        lnewt(n)=.false.
        if (igrd.ne.0) then
          npigon(igrd) = npigon(igrd) + 1
          do l = 1,3
            pigmass(l,igrd) =  pigmass(l,igrd) + puffmass(l,n)
          enddo
        endif
      enddo
c
      write(idiag,'(/,a,f10.0,i10.5)')'PiG diagnostics at: ',time,date 
      do igrd = 1, ngrid 
        write(idiag,'(a,i5,a,i5)') '# active PiGs in grid: ',igrd, 
     &                             ' is ',npigon(igrd) 
        write(idiag,'(a,3(1pe10.3))') 
     &             '   Total PiG NO, NO2, HNO3 mass: ', 
     &              pigmass(1,igrd),pigmass(2,igrd),pigmass(3,igrd) 
        write(idiag,'(a,3(1pe10.3))') 
     &              'Total dumped NO, NO2, HNO3 mass: ', 
     &              pigdump(1,igrd),pigdump(2,igrd),pigdump(3,igrd) 
      enddo 
c
      return
      end
