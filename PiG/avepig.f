      subroutine avepig(igrd,dt,ncol,nrow,nlay,nlayav,dx,dy,mapscl,
     &                  height,nspav,nspc,lmap,tempk,press,conc,avcnc)
      real mapscl(ncol,nrow)
      dimension lmap(nspc)
c
c-----CAMx v4.02 030709
c
c
      if (igrd.ne.999) then
        write(*,'(A)') 'The PiG option IRON is not yet supported'
        write(*,'(A)') 'Choose PiG option GREASD or NONE'
        call camxerr()
      endif
c
      return
      end
