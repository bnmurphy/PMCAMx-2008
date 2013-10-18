      subroutine ironmscl(ngrid,time,date,idiag,pigdump)
c
      integer date
      include 'camx.prm'
      real*8  pigdump(MXSPEC,MXGRID)
c
c-----CAMx v4.02 030709
c
c
      if (ngrid.ne.999) then
        write(*,'(A)') 'The PiG option IRON is not yet supported'
        write(*,'(A)') 'Choose PiG option GREASD or NONE'
        call camxerr()
      endif
c
      return
      end
