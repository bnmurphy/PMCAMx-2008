      program temprd
c######################################################################
c
c    Temprd opens PMCAMx layer temperature input files and prints the data
c	to be read and plotted by Matlab
c	Output File Name:  temp.daily.[layer]
c
c    Written by Ben Murphy - 6-19-2008
c
c#####################################################################
c
      implicit none

      integer ifin, ifout, ifout2
      integer ndate, ncol, nrow, nlay, nhr

      parameter (ifin  = 20)
      parameter (ifout = 30)
      parameter (ifout2= 31)
      parameter (ndate = 17)
      parameter (ncol  = 97)
      parameter (nrow  = 90)
      parameter (nlay  = 14)
      parameter (nhr   = 25)

      character*3   	cdate(ndate), chr(24)
      character*2       clay(nlay)
      real		hour
      integer		idate, ilay
      character*198 	filename, indir, outdir, fname, fname2, fnout
      real 		temp(ncol,nrow,nlay,nhr,ndate), temps(ncol,nrow,nhr,ndate)
      real		tee

      integer d, io, p(500000), ihr
      integer k,n,l,i,t,j,ll, icol, irow, lPIT, lDiurn, lSurf
      real a, b

c
c  Specify Date Vectors
c
      data cdate /'193','194','195','196','197','198','199','200','201',
     &		  '202','203','204','205','206','207','208','209'/
c      data dates /'193','194'/

      data chr /'01','02','03','04','05','06','07','08','09','10','11',
     &	       	'12','13','14','15','16','17','18','19','20','21','22',
     &		'23','24'/
      data clay /'01','02','03','04','05','06','07','08','09','10','11',
     &		'12','13','14'/

c
c  Set Directory Structure
c
c temp.2001$DATE.4rpos.36.14.mm5.ld.camx

      indir  = '/home/BaseInput/July/temp/'
      fname  = 'temp.2001'
      fname2 = '.4rpos.36.14.mm5.ld.camx'
      fnout  = 'temp'
      outdir = '/home/mday/processed_output/2001meteorology/' !'/home/mday/Temp/Output/'


c
c  Begin Date Loop
c
      do d=1,ndate
	print *,'Day: ',cdate(d)

c
c      Open Input File
c
         filename = indir(1:INDEX(indir,' ')-1)// 
     &		fname(1:INDEX(fname,' ')-1)//
     &		cdate(d)
     &          // fname2(1:INDEX(fname2,' ')-1)
         open (UNIT=ifin,FILE=filename,FORM='UNFORMATTED',STATUS='OLD')

c
c      Read Time-Variant Portion
c
	do ihr = 1,25
	    read(ifin), hour,idate,((temps(i,j,ihr,d),i=1,ncol),j=1,nrow)
	    do k = 1,nlay
	    	read(ifin,end=200), hour, idate, ((temp(i,j,k,ihr,d),
     &				i=1,ncol),j=1,nrow)
	    enddo
	enddo

 200    continue
        close(ifin)

      enddo !date

c
c  Tally and Print Average Temperature Profile
c
      do k = 1,nlay
	tee = 0
	do d = 1,ndate
	  do i = 1,ncol
	    do j = 1,nrow
	      do ihr = 1,24
		tee = tee + temp(i,j,k,ihr,d)/24/ncol/nrow/ndate
	      enddo
	    enddo
	  enddo
	enddo
	print *,'Average T (Layer',k,') = ',tee
      enddo

c
c  Output Data
c
      lSurf = 0
      if (lSurf.eq.1) then
c    Write Surface Layer Files
          filename = outdir(1:INDEX(outdir,' ')-1)//
     &			'tempsurf.diurn.01'
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'HOUR','X','Y','Temp(K)'
	 print *
	 print *,'Filename: ',filename

c	 Average diurnal trend for all days and write map data
	   do ihr = 1,24
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do idate = 1,ndate
			a = a + (temps(i,j,ihr,idate)+temps(i,j,ihr+1,idate))/2/ndate
		    enddo
		    write (ifout,901), chr(ihr), i, j, a
		enddo  !j
	     enddo  !i
	    enddo  !ihr
	close(ifout)

c      Open Daily Average Map Output Files
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			'tempsurf.daily.01'
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'DATE','X','Y','TEMP(K))'
	 print *
	 print *,'Filename: ',filename

c	 Average diurnal trend for all days and write map data
	   do idate = 1,ndate
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do ihr = 1,24
			a = a + (temps(i,j,ihr,idate)+temps(i,j,ihr+1,idate))/2/24
		    enddo
		    write (ifout,901), cdate(idate), i, j, a
		enddo  !j
	     enddo  !i
	    enddo  !idate
	close(ifout)
      endif  !lSurf

c    Begin Layer Loop
      do ilay = 1,nlay
	lDiurn = 0
	if(lDiurn.eq.1) then
c
c      Open Diurnal Map Output Files
c
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout(1:INDEX(fnout,' ')-1)//'.diurn.'//clay(ilay)
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'HOUR','X','Y','Temp(K)'
	 print *
	 print *,'Filename: ',filename

c	 Average diurnal trend for all days and write map data
	   do ihr = 1,24
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do idate = 1,ndate
			a = a + (temp(i,j,ilay,ihr,idate)+temp(i,j,ilay,ihr+1,idate))/2/ndate
		    enddo
		    write (ifout,901), chr(ihr), i, j, a
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
         write(ifout,900), 'DATE','X','Y','Temp(K)'
	 print *
	 print *,'Filename: ',filename

c	 Average diurnal trend for all days and write map data
	   do idate = 1,ndate
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do ihr = 1,24
			a = a + (temp(i,j,ilay,ihr,idate)+temp(i,j,ilay,ihr+1,idate))/2/24
		    enddo
		    write (ifout,901), cdate(idate), i, j, a
		enddo  !j
	     enddo  !i
	    enddo  !idate
	close(ifout)
	endif   !lDiurn

        lPIT = 1
      if (lPIT.eq.1) then
       filename= outdir(1:INDEX(outdir,' ')-1)// fnout(1:INDEX(fnout,' ')-1)//
     &          '.ATL.'//clay(ilay) !.PIT.
        open (unit = 20, file=filename, status='unknown')
        write(20,903), 'DATE','HOUR','X','Y','Temp (K)'
          do idate = 1,ndate
          do ihr = 1,24
                icol = 58 !65
                irow = 29 !51
                    write (20,904), cdate(idate),chr(ihr), icol, irow, temp(icol,irow,ilay,ihr,idate)
          enddo
          enddo
        close(20)
      endif

       enddo !ilay

 300    format(1X,F10.7)

 900  format(A4,2x,A2,2x,A2,4x,A12)
 901  format(A3,3x,I2,4x,I2,4x,E10.4)
 903  format(A4,2x,A4,2x,A2,2x,A2,2x,A12)
 904  format(A4,2x,A4,2x,I2,2x,I2,2x,E10.4)

      end
 
