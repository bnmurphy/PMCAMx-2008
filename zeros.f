      subroutine zeros(a,n)
c
c-----CAMx v4.02 030709
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001 , 2002, 2003
c     ENVIRON International Corporation
c
c     ZEROS sets input array to zero
c
      dimension a(n)
c
      do i=1,n
        a(i) = 0.
      enddo
c
      return
      end
