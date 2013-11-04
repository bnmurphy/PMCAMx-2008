      program vdifrd07
c######################################################################
c
c    Vdifrd opens PMCAMx layer Kzz input files and prints the data
c	to be read and plotted by Matlab
c	Output File Name:  vdif.daily.[layer]
c
c    Written by Ben Murphy, tweaked by Melissa Day
c
c#####################################################################
c
      implicit none

      integer ifin, ifout, ifout2
      integer ndate, ncol, nrow, nlay, nhr

      parameter (ifin  = 20)
      parameter (ifout = 30)
      parameter (ifout2= 31)
      parameter (ndate = 31)
      parameter (ncol  = 148)
      parameter (nrow  = 112)
      parameter (nlay  = 16)
      parameter (nhr   = 25)

      character*3   	cdate(ndate), chr(24)
      character*2       clay(nlay)
      real		hour
      integer		idate, ilay
      character*198 	filename, indir, outdir, fname, fname2, fnout
      real 		vdif(ncol,nrow,nlay,nhr,ndate), vdifs(ncol,nrow,nhr,ndate)
      real		tee, Kzz

      integer d, io, p(500000), ihr
      integer k,n,l,i,t,j,ll, icol, irow, lPIT, lDiurn, lSurf
      real a, b

c
c  Specify Date Vectors
c
      data cdate /'182','183','184','185','186','187','188','189','190','191',
     &            '192','193','194','195','196','197','198','199','200','201',
     &		  '202','203','204','205','206','207','208','209','210',
     &            '211','212'/
c      data dates /'193','194'/

      data chr /'01','02','03','04','05','06','07','08','09','10','11',
     &	       	'12','13','14','15','16','17','18','19','20','21','22',
     &		'23','24'/

      data clay /'01','02','03','04','05','06','07','08','09','10','11',
     &		'12','13','14','15','16'/

c
c  Set Directory Structure
c
      indir  = '/home/mday/2007_input/input/vdif/' !'/home/BaseInput/July/temp/'
      fname  = 'vdif.2007'
      fname2 = '.conus.36.16.wrf.ld.camx.OB70mesh' !'.4rpos.36.14.mm5.ld.camx'
      fnout  = 'vdif'
      outdir = '/home/mday/2007_input/output/vdif/' !'/home/mday/Temp/Output/'


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
         print *, 'Input file: ',filename

c
c      Read Time-Variant Portion
c
	do ihr = 1,nhr
c	    read(ifin), hour,idate,((temps(i,j,ihr,d),i=1,ncol),j=1,nrow)
	    do k = 1,nlay
	    	read(ifin,end=200), hour, idate, ((vdif(i,j,k,ihr,d),
     &				i=1,ncol),j=1,nrow)
	    enddo
	enddo

 200    continue
        close(ifin)

      enddo !date
      print *,'Reading done..'
c
c  Tally and Print Vertical Dispersion Coeffiicents
c
      do k = 1,nlay
        Kzz = 0
       do ihr = 1,nhr
	tee = 0
	  do i = 1,ncol
	    do j = 1,nrow
	      do d = 1,ndate
		tee = tee + vdif(i,j,k,ihr,d)/ncol/nrow/ndate
	      enddo
	    enddo
	  enddo
	enddo
	print *,'Layer:',k,' Hour:',ihr,' Kzz:', tee
      enddo

c
c  Output Data
c

c    Begin Layer Loop
      do ilay = 1,nlay
	lDiurn = 1
	if(lDiurn.eq.1) then
c
c      Open Diurnal Map Output Files
c
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout(1:INDEX(fnout,' ')-1)//'.diurn.'//clay(ilay)
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'HOUR','X','Y','Kzz(m^2/s)'
	 print *
	 print *,'Filename: ',filename

c	 Average diurnal trend for all days and write map data
	   do ihr = 1,24
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do idate = 1,ndate
			a = a + (vdif(i,j,ilay,ihr,idate)+vdif(i,j,ilay,ihr+1,idate))/2/ndate
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
         write(ifout,900), 'DATE','X','Y','Kzz(m^2/s)'
	 print *
	 print *,'Filename: ',filename

c	 Average diurnal trend for all days and write map data
	   do idate = 1,ndate
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
		    do ihr = 1,24
			a = a + (vdif(i,j,ilay,ihr,idate)+vdif(i,j,ilay,ihr+1,idate))/2/24
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
     &          '.PIT.'//clay(ilay)
        print *,'Filename: ',filename
        open(UNIT = ifout2, FILE=filename, STATUS='unknown')
        write(ifout2,903), 'DATE','HOUR','X','Y','Kzz(m^2/s)'
          do idate = 1,ndate
          do ihr = 1,24
            icol = 65
            irow = 51
              write (ifout2,904), cdate(idate),chr(ihr), icol, irow, vdif(icol,irow,ilay,ihr,idate)
          enddo
          enddo
        close(ifout2)

        filename= outdir(1:INDEX(outdir,' ')-1)// fnout(1:INDEX(fnout,' ')-1)//
     &          '.ATL.'//clay(ilay)
        print *,'Filename: ',filename
        open(UNIT = ifout2, FILE=filename, STATUS='unknown')
        write(ifout2,903), 'DATE','HOUR','X','Y','Kzz(m^2/s)'
          do idate = 1,ndate
          do ihr = 1,24
            icol = 58
            irow = 29
              write (ifout2,904), cdate(idate),chr(ihr), icol, irow, vdif(icol,irow,ilay,ihr,idate)
          enddo
          enddo
        close(ifout2)
      endif  !lPIT

       enddo !ilay

 300    format(1X,F10.7)

 900  format(A4,3x,A2,3x,A2,4x,A12)
 901  format(A3,3x,I3,2x,I3,2x,E10.4)
 903  format(A4,2x,A4,2x,A3,2x,A3,2x,A12)
 904  format(A4,2x,A4,2x,I2,2x,I2,2x,E10.4)

      end
