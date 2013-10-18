      subroutine irondriv(igrd,iptr2d,ncol,nrow,nlay,lrad,itzon,
     &                    dt,dx,dy,mapscl,height,rkv,tempk,press,water,
     &                    windu,windv,cldtrns,fcloud,cellat,cellon,
     &                    conc,cncrad,pigdump)
      dimension cncrad(ncol,nrow,nlay,*)
      real*8 pigdump(*)
      real mapscl(ncol,nrow)
c
      if (igrd.ne.999) then
        write(*,'(A)') 'The PiG option IRON is not yet supported'
        write(*,'(A)') 'Choose PiG option GREASD or NONE'
        call camxerr()
      endif
c
      return
      end
