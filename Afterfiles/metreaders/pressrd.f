      program pressrd
c######################################################################
c
c    Pressrd opens PMCAMx layer height input files and prints the data
c	to be read and plotted by Matlab
c	Output File Name:  hgt.layer
c
c    Written by Ben Murphy - 6-19-2008
c
c#####################################################################
c
      implicit none

      integer ifin,ifouth, ifout, ifout2
      integer ndate, ncol, nrow, nlay, nhr

      parameter (ifin  = 20)
      parameter (ifouth = 29)
      parameter (ifout = 30)
      parameter (ifout2= 31)
      parameter (ndate = 9)
      parameter (ncol  = 97)
      parameter (nrow  = 90)
      parameter (nlay  = 14)
      parameter (nhr   = 25)

      character*3   	cdate(ndate), chr(24)
      character*2       clay(nlay)
      real		hour
      integer		idate, ilay
      character*198 	filename, indir, outdir, fname, fname2, fnout, fnout2, hfile
      real 		height(ncol,nrow,nlay,nhr,ndate), press(ncol,nrow,nlay,nhr,ndate)
      real              tee,h

      integer d, io, p(500000), ihr
      integer k,n,l,i,t,j,ll
      real a, b

c
c  Specify Date Vectors
c
c      data cdate /'193','194','195','196','197','198','199','200','201',
c     &		  '202','203','204','205','206','207','208','209'/
      data cdate /'194','195','196','197','198','199','200','201','202'/
c       data cdate /'001','002','003','004','005','006','007','008','009',
c     &		   '010','011','012','013','014','015','016','017'/
c	data cdate /'356','357','358','359','360','361','362','363',
c     &		    '364','365'/

      data chr /'01','02','03','04','05','06','07','08','09','10','11','12',
     &		'13','14','15','16','17','18','19','20','21','22','23','24'/

      data clay /'01','02','03','04','05','06','07','08','09','10','11','12',
     &		 '13','14'/

c
c  Set Directory Structure
c
c      indir  = '/home2/bnmurphy/PMCAMx_Input/Base_Input/2003/January/pres/'
c      fname  = 'camx.zp.2003'
c      fname2 = '.updt'
c      fnout  = 'jan.height'
c      fnout2 = 'jan.press'

c      indir  = '/home2/bnmurphy/PMCAMx_Input/Base_Input/2003/July/pres/'
c      fname  = 'camx.zp.2003'
c      fname2 = '.updt'
c      fnout  = 'july03.height'
c      fnout2 = 'july03.press'
c      outdir = '/home2/bnmurphy/PMCAMx_Input/mm5camx_Ben_IO/confirm/out/'

      indir  = '/home/BaseInput/July/pres/'
      fname  = 'pres.2001'
      fname2 = '.4rpos.36.14.mm5.ld.camx'
      fnout  = 'july01.height'
      fnout2 = 'july01.press'
      outdir = '/home/mday/processed_output/height/'

c
c Open output avg file
c

	hfile = outdir(1:INDEX(outdir,' ')-1)//'PandH.avg'
        open(unit=ifouth,file=hfile,status='UNKNOWN')
	write(ifouth,*)'	Layer ','	Pressure ','	Height '
c
c  Begin Date Loop
c
      do d=1,ndate
	print *,'Day: ',cdate(d)

c
c      Open Input File
c
         filename = indir(1:INDEX(indir,' ')-1)// fname(1:INDEX(fname,' ')-1)//
     &		   cdate(d)//fname2(1:INDEX(fname2,' ')-1)
         open (UNIT=ifin, FILE=filename, FORM='UNFORMATTED', STATUS='OLD')

c
c      Read Time-Variant Portion
c
	do ihr = 1,25
	    do k = 1,nlay
	    	read(ifin,end=200), hour, idate, ((height(i,j,k,ihr,d),
     &				i=1,ncol),j=1,nrow)
	    	read(ifin), hour, idate, ((press(i,j,k,ihr,d),i=1,ncol),j=1,nrow)
	    enddo
	enddo

 200    continue
        close(ifin)

      enddo !date

c
c  Tally and Print Average Pressure Profile
c
      do k = 1,nlay
        tee = 0
	h = 0
        do d = 1,ndate
          do i = 1,ncol
            do j = 1,nrow
              do ihr = 1,24
                tee = tee + press(i,j,k,ihr,d)/24/ncol/nrow/ndate
                h = h + height(i,j,k,ihr,d)/24/ncol/nrow/ndate
              enddo
            enddo
          enddo
        enddo
	write(ifouth,*)k,tee,h
        print *,'Average Pressure (Layer',k,') = ',tee
       	print *,'Average Height (Layer',k,') = ',h
      enddo


c
c  Output Data
c
c    Begin Layer Loop
      do ilay = 1,nlay

c
c      Open Diurnal Map Output Files
c
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout(1:INDEX(fnout,' ')-1)//'.diurn.'//clay(ilay)
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'HOUR','X','Y','Height(m)'
	 !print *
	 !print *,'Filename: ',filename

         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout2(1:INDEX(fnout2,' ')-1)//'.diurn.'//clay(ilay)
         open (UNIT = ifout2, FILE=filename, STATUS='unknown')
         write(ifout2,900), 'HOUR','X','Y','Pressure(mb)'
	    
c	 Average diurnal trend for all days and write map data
	   do ihr = 1,24
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do idate = 1,ndate
			a = a + (height(i,j,ilay,ihr,idate)+height(i,j,ilay,ihr+1,idate))/2/ndate
			b = b + (press(i,j,ilay,ihr,idate)+press(i,j,ilay,ihr+1,idate))/2/ndate
		    enddo
		    write (ifout,901), chr(ihr), i, j, a
		    write (ifout2,901), chr(ihr), i, j, b
		enddo  !j
	     enddo  !i
	    enddo  !ihr
	close(ifout)
	close(ifout2)

c
c      Open Daily Average Map Output Files
c
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout(1:INDEX(fnout,' ')-1)//'.daily.'//clay(ilay)
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'DATE','X','Y','Height(m)'
	 !print *
	 !print *,'Filename: ',filename

         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout2(1:INDEX(fnout2,' ')-1)//'.daily.'//clay(ilay)
         open (UNIT = ifout2, FILE=filename, STATUS='unknown')
         write(ifout2,900), 'DATE','X','Y','Pressure(mb)'
	    
c	 Average daily trend for all days and write map data
	   do idate = 1,ndate
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do ihr = 1,24
			a = a + (height(i,j,ilay,ihr,idate)+height(i,j,ilay,ihr+1,idate))/2/24
			b = b + (press(i,j,ilay,ihr,idate)+press(i,j,ilay,ihr+1,idate))/2/24
		    enddo
		    write (ifout,901), cdate(idate), i, j, a
		    write (ifout2,901), cdate(idate), i, j, b
		enddo  !j
	     enddo  !i
	    enddo  !idate
	close(ifout)
	close(ifout2)


       enddo !ilay

 300    format(1X,F10.7)

 900  format(A4,2x,A2,2x,A2,2x,A12)
 901  format(A3,3x,I2,4x,I2,4x,E10.4)
 902  format(A3,3x,I2,4x,I2,4x,E10.4)

      end
 
