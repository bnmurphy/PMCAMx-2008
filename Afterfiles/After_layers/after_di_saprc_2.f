      SUBROUTINE  after_di_saprc(outdir,scenario,nlay,clay,ilay,ndate,
     &				 cdate,mspc,cspc,nspc,mconc)
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
c      parameter(mspc = 88)    ! Number of output species
      parameter(gasfac = 1000)
      parameter(nsite = 9)

      integer MW
      integer ndate, mspc, nspc
      integer isite
      integer nlay, ilay
      character*2   clay(nlay)       !Character for layer

      parameter(R    = 83.145)    !cm^3*bar/mol/K
      parameter(T    = 300)	  !K
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
      real*4 sconh(mspc,ncol,nrow,nhr)      !Conc for output (avg, sum, etc)
      real*4 date1, a			    !Dummies for output algorithm
      
      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate,col(nsite),row(nsite)
      integer*4 ihr,msec,icut,line
      integer*4 iday,ix,iy,iz,io


c      data cdate /'193','194','195','196','197','198','199','200',
c     &	 	  '201','202','203','204','205','206','207','208',
c     &		  '209'/
     
      data chr /'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &          '21','22','23','24'/


c     Sites = Pittsburgh, St. Louis, New York City, Atlanta, Chicago
c     Lat   = 40.45       38.65,     40.78,         33.65,   41.88
c     Lon   = -80.0,      -90.63,    -73.96,        -84.42,  -87.63

c      data cext /'PIT','STL','NYC','ATL','CHI','DFN','PSP','NE1','NE2','HOU','TOT'/
      data cext /'PIT','STL','NYC','ATL','CHI','DFN','PSP','HOU','TOT'/
c      data cext /'TOT','PIT','ATL','DFN','LMI','NYC','CHI'/
c      data cext /'TOT','LMI','NYO','LER','NYC','CHI'/
c      data cext /'FIN','TOT'/
c      data cext /'TOT'/

c      data col /65, 41, 79, 58, 47, 70, 72, 74, 84, 30, 1/
      data col /65, 41, 79, 58, 47, 70, 72, 30, 1/
c      data col /1, 65, 58, 70, 47, 79, 47/
c      data col /1, 47, 83, 58, 79, 47/
c      data col /132, 1/
c      data col /1/

c      data row /51, 42, 55, 28, 52, 38, 60, 28, 64, 14, 1/
      data row /51, 42, 55, 28, 52, 38, 60, 14, 1/
c      data row /1, 51, 28, 38, 56, 55, 52/
c      data row /1, 56, 55, 53, 55, 52/
c      data row /84, 1/
c      data row /1/

c
c  LOOP THROUGH OUTPUT FILES
c
      print *,'Exporting Diurnal Data to directory: ',outdir(1:INDEX(outdir,' ')-1),'\n'
      do isite = 1,nsite
      do io = 1,mspc

c 	Open Output File
	  print '(A37,A5,A1,A5)','Opening Diurnal Profile Output File ',cspc(io),'-',cext(isite),'...'

	  isblk = INDEX(cspc(io),' ')
          if(isblk.eq.0) isblk = 5
          fname = outdir(1:INDEX(outdir,' ')-1)//
     &		   scenario(1:INDEX(scenario,' ')-1)//'/'//
     &		   cspc(io)(1:isblk-1)//'.'//cext(isite)//'.'//clay(ilay)

	  isblk = INDEX(fname,' ')
	  open(ifout,file=fname(1:isblk-1),status='UNKNOWN')
	  write(ifout, 46), 'Date', 'Hour', 'X','Y','CONC'

        
c 	Write to Output File
        do idate=1,ndate
	  do ihr = 1,nhr
	    if (cext(isite).eq.'TOT') then
	    	a = 0
	    	do i = 1,ncol
		    do j = 1,nrow
		    	a = a + mconc(idate,ihr,i,j,io)/8730
		    enddo
		enddo
	  else 
	    i=col(isite)     ! Column First
	    j=row(isite)     ! Then Row
	    a = mconc(idate,ihr,i,j,io)
	  endif
          write( ifout, 47),cdate(idate), ihr,i,j,a
	    end do
      end do  !idate

c  	Close Output Files
	  print *,'    Closing Diurnal Output File...'
	  close (ifout)	!Close Output Files

      end do  !io
      end do  !isite


c
c
 46    format(A4,2x,A4,2x,A1,3x,A1,3x,A4)
 47    format(A3,3x,i2,4x,i2,2x,i2,2x,e12.5)

       write(6,*)'ALL DAYS DONE'

 990   end
