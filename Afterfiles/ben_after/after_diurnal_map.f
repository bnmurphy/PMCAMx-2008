      SUBROUTINE  after_diurnal_map(outdir,scenario,nlay,clay,ilay,ndate,
     &                              cdate,mspc,cspc,nspc,mconc)
c######################################################################
c   BNM
c     This program processes PMCAMx results
c     The parameters are hardcoded for the specific output file
c     Ben Murphy rewrote this post-processor 3-11-08
c       outputs Gas, PM 2.5, and NOT(PM 10) averaged over each day
c
c     Modifications
c	nmol to ppb conversion for CPO,COO,CAS,CBS - (7-1-08)
c
c######################################################################

c
c Define constants and arrays
c

      parameter(ifin  = 10)
      parameter(ifout = 30)

      parameter(ncol = 97)
      parameter(nrow = 90)
      parameter(nhr  = 24)    ! Hours in each day
      parameter(gasfac = 1000)
      parameter(nsite = 4)

      integer MW
      integer ndate, mspc, nspc
      integer isite
      integer nlay, ilay
      character*2   clay(nlay)       !Character for layer

      parameter(R    = 83.145)    !cm^3*bar/mol/K
      parameter(T    = 300)       !K
      parameter(P    = 1.01325)   !bar

      character*198 fname1
      character*198 outdir, indir, indir1
      character*190 fname
      character*198 scenario   !the name of the output folder
      character*198 scenario2  !name of input file

      character*4   name(10),note(60),mspec(10,nspc)
      character*2   chr(nhr)         !Hour
      character*3   cdate(ndate)   !Julian Date
      character*10  spcname(nspc)    !Input species name string
      character*7   cspc(mspc)      !Output species name vector
      character*3   cext(nsite)

      real*4 maxh(mspc)
      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 mconc(ndate,31,ncol,nrow,mspc)  !Conc read from input file
      real*4 sconh(mspc,ncol,nrow,nhr)       !Conc for output (avg, sum, etc)
      real*4 date1, a                        !Dummies for output algorithm
      
      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate,col(nsite),row(nsite)
      integer*4 ihr,msec,icut,line
      integer*4 iday,ix,iy,iz,io

c
c Open Output Files
c
	  print *,'\n','Opening Diurnal Profile Output Files...','\n\n'

	  do io = 1,mspc

	    isblk = INDEX(cspc(io),' ')
          if(isblk.eq.0) isblk = 5
          fname = outdir(1:INDEX(outdir,' ')-1)//
     &		   scenario(1:INDEX(scenario,' ')-1)//'/'//
     &		   cspc(io)(1:isblk-1)//'.diurn.'//clay(ilay)

	    isblk = INDEX(fname,' ')
	    open((ifout+io),file=fname(1:isblk-1),status='UNKNOWN')
	    write(ifout+io, 46), 'Hour', 'X','Y','CONC'

	  end do


c
c  Write to Output Files
c
      print *,'Exporting Data to directory: ',outdir(1:INDEX(outdir,' ')-1)
      do io = 1,mspc
      do ihr = 1,nhr
        do i = 1,ncol
	  do j = 1,nrow
	    a = 0
	    do idate=1,ndate
		a = a + mconc(idate,ihr,i,j,io)/ndate
	    enddo
                write( (ifout+io), 47),ihr,i,j,a
	  enddo
	  end do
      end do  !io
      end do  !idate


c
c  Close Output Files
c
	  print *,'Closing Dirunal Output Files...','\n'
	  do io = 1,mspc
		close (ifout+io)	!Close Output Files
	  end do

  46   format(A4,2x,A1,3x,A1,3x,A4)
  47   format(i2,4x,i2,2x,i2,2x,e12.5)

       write(6,*)'ALL DAYS DONE'

 990   end
