      subroutine aero_err()
c   
c-----PMCAMx v3.01 020531
c   
c     CAMXERR writes the final message whenever CAMx terminates
c             due to an error
c                             
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001 
c     ENVIRON International Corporation
c             
c     Modifications:   
c        none 
c 
c     Input arguments:
c        none 
c 
c     Output arguments: 
c        none 
c               
c     Routines Called:   
c        none 
c               
c     Called by:   
c        LININT
c 
      include "camx.prm"
      include "chmstry.com"
      include "camxfld.com"
      include "filunit.com"
      include "grid.com"
c
      data densfac /44.9/
c
c-----Entry point
c
c-----Calcutate unit conversion factor-- umol/m^3 -> ppm
c
      convfac =  densfac*(273./298.)
c
c-----Calculate pointers and output temp,press,conc,cncrad,etc...
c
      if ( igrdchm .eq. 1 ) then
	n3d=ichm + ncol(igrdchm)*(jchm-1) + 
     &           ncol(igrdchm)*nrow(igrdchm)*(kchm-1)
      endif

      write(*,*) 
      write(*,*) 
      write(*,*) ' CAMx is stopping because an error has occured in the aerosol routines '
      write(*,*) ' See the .out output file for details'
      write(*,*) ' igrd,i,j,k :',igrdchm,ichm,jchm,kchm
      write(*,*) 
      write(*,*) 
c
      write(*,*) 'Temperature :',tempk(iptr3d(igrdchm)-1+n3d)
      write(*,*) 'Pressure    :',press(iptr3d(igrdchm)-1+n3d)
      write(*,*) 'Water       :',water(iptr3d(igrdchm)-1+n3d)
cbk      write(*,*) 'LWC         :',lwc(iptr3d(igrdchm)-1+n3d)
      write(*,*)
      write(*,*) ' Concentrations are :'
      do l=1,ngas
	n4d = n3d + ncol(igrdchm)*nrow(igrdchm)*nlay(igrdchm)*(l-1)
	write(*,'(a10,1x,e13.6)') spname(l),
     &                       conc(iptr4d(igrdchm)-1+n4d)/convfac
      enddo
      do l=ngas+1,nspec
	n4d = n3d + ncol(igrdchm)*nrow(igrdchm)*nlay(igrdchm)*(l-1)
	write(*,'(a10,1x,e13.6)') spname(l),conc(iptr4d(igrdchm)-1+n4d)
      enddo
c
      write(*,*)
      write(*,*) ' Radicals are :'
      do l=1,nrad
	n4d = n3d + ncol(igrdchm)*nrow(igrdchm)*nlay(igrdchm)*(l-1)
	write(*,'(a10,1x,e13.6)') nmrad(l),cncrad(iptrad(igrdchm)-1+n4d)
      enddo
c
      write(iout,*) 
      write(iout,*) 
      write(iout,*) ' CAMx is stopping because of the error(s) ',
     &              'described above'
      write(iout,*) 
      write(iout,*) 
c

c      do i=1,ncol(igrdchm)
c       do j=1,nrow(igrdchm)
c	n3d=i + ncol(igrdchm)*(j-1) + 
c     &           ncol(igrdchm)*nrow(igrdchm)*(kchm-1)
c        do l=1,nspec
c	  n4d = n3d + ncol(igrdchm)*nrow(igrdchm)*nlay(kchm)*(l-1)
c	  write(*,'(a10,3(1x,i3),1x,e12.6)') spname(l),i,j,kchm,conc(iptr4d(igrdchm)-1+n4d)
c	enddo
c       enddo
c      enddo
c	  

      stop
c
      end

      subroutine get_param(igrdchm_c,ichm_c,jchm_c,kchm_c,
     &                                      iout_c,idiag_c)

      include "camx.prm"
      include "chmstry.com"
      include "filunit.com"

      integer igrdchm_c,ichm_c,jchm_c,kchm_c,iout_c,idiag_c

      igrdchm_c = igrdchm
      ichm_c    = ichm
      jchm_c    = jchm
      kchm_c    = kchm
      iout_c    = iout
      idiag_c   = idiag

      return
      end

