      program rainrd
c######################################################################
c
c    Rainrd opens PMCAMx layer height input files and prints the data
c	to be read and plotted by Matlab
c	Output File Name:  rain.diurn.[layer]
c
c    Written by Ben Murphy - 4-20-2009
c
c#####################################################################
c
      implicit none

      integer ifin, ifout, ifout2
      integer ndate, ncol, nrow, nlay, nhr, nyr

      parameter (ifin  = 20)
      parameter (ifout = 30)
      parameter (ifout2= 31)
      parameter (ndate = 17) !31
      parameter (ncol  = 97)
      parameter (nrow  = 90)
      parameter (nlay  = 14)
      parameter (nhr   = 25)
      parameter (nyr   = 1)  !20

      character*3   cdate(ndate), chr(24)
      character*2   clay(nlay)
      character*4   cyr(nyr), year
      character*15  cldhdr
      integer       cldcol, cldrow, cldlay, dum
      real	    hour, date
      integer	    idate, ilay
      character*198 filename, indir, outdir, fname, fname2, fnout2
      character*6   fnout(5)
      real 	    cldwtr(ncol,nrow,nlay,nhr,ndate), ranwtr(ncol,nrow,nlay,nhr,ndate)
      real 	    snowtr(ncol,nrow,nlay,nhr,ndate), gplwtr(ncol,nrow,nlay,nhr,ndate)
      real 	    cldopd(ncol,nrow,nlay,nhr,ndate)
      real 	    tee, tee2, tee3, tau, coefcld, ctrns, energy(nlay), cloud(nlay)
      real          fcld, zenang, cldrat, flcoud
      integer       io, p(500000), ihr, icol, irow, iabov, iyr
      integer       k,n,l,i,t,j,ll
      integer       lexplicit, lPIT, lDiurn, lDaily
      real          a,b,c,d,e
      real          e

c      data cyr /'1991','1992','1993','1994','1995','1996','1997','1998','1999','2000'
c     &         ,'2051','2052','2053','2054','2055','2056','2057','2058','2059','2060'/
     data cyr /'2001'/

c      data cdate /'182','183','184','185','186','187','188','189','190','191','192',
c     &            '193','194','195','196','197','198','199','200','201','202','203',
c     &            '204','205','206','207','208','209','210','211','212'/

      data cdate /'193','194','195','196','197','198','199','200','201',
     &		  '202','203','204','205','206','207','208','209'/
c      data cdate /'356','357','358','359','360','361','362','363','364','365'/

      data chr /'01','02','03','04','05','06','07','08','09','10','11','12',
     &		'13','14','15','16','17','18','19','20','21','22','23','24'/
      data clay /'01','02','03','04','05','06','07','08','09','10','11','12','13','14'/
      data fnout /'cldwtr','ranwtr','snowtr','gplwtr','cldopd'/

c  Define year
      do iyr = 1,nyr
      year = cyr(iyr)
      print *, 'Starting ',year

c
c  Set Directory Structure
c
c  /home/BaseInput/July/clranew/clra.2001$DATE.4rpos.36.14.mm5.ld.camx
c  /home/PMCAMx/met_hd/CAMx_Input/$YEAR/camx.cr.$YEAR$DATE.updt

      indir  = '/home/BaseInput/July/clranew/'
      outdir = '/home/mday/processed_output/2001meteorology/'
      fname  = 'clra.2001'
      fname2 = '.4rpos.36.14.mm5.ld.camx'
      fnout2 = '.Jul01'

c      indir  = '/home/mday/PMCAMx/met_hd/CAMx_Input/' // year // '/'
c      outdir = '/home/mday/PMCAMx/met_hd/CAMx_Input/' // year // '/'
c      fname  = 'camx.cr.' // year      !'clra.2001'
c      fname2 = '.updt'                 !'.4rpos.36.14.mm5.ld.camx'
c      fnout2 = '.Jul'//year            !'.Jul01'

c
c  Begin Date Loop
c
      do idate = 1,ndate
	print *,'Day: ',cdate(idate)

c
c      Open Input File
c
         filename = indir(1:INDEX(indir,' ')-1)// fname(1:INDEX(fname,' ')-1)//
     &		   cdate(idate)//fname2(1:INDEX(fname2,' ')-1)
         print *, filename
         open (UNIT=ifin, FILE=filename, FORM='UNFORMATTED', STATUS='OLD')

c
c      Read Header
c
	read(ifin) cldhdr, cldcol, cldrow, cldlay
	print *,'Input file header: ',cldhdr
        print *,'Input file has ',cldcol,' cols, ',cldrow,' rows, and ',cldlay,' lays.'
c
c      Read Time-Variant Portion
c
	do ihr = 1,25
	    read(ifin,end=200), hour, dum
	    do k = 1,nlay
	    	read(ifin),((cldwtr(i,j,k,ihr,idate),i=1,ncol),j=1,nrow) !g m-3
	    	read(ifin),((ranwtr(i,j,k,ihr,idate),i=1,ncol),j=1,nrow) !g m-3
c               Convert ranwtr to rainfall rate
		do icol = 1,ncol
		  do irow = 1,nrow
		    ranwtr(icol,irow,k,ihr,idate) = (ranwtr(icol,irow,k,ihr,idate) /1e6 *1e7)**1.14  !mm hr-1
		  enddo
		enddo
c	    	read(ifin),((snowtr(i,j,k,ihr,idate),i=1,ncol),j=1,nrow) !g m-3
c	    	read(ifin),((gplwtr(i,j,k,ihr,idate),i=1,ncol),j=1,nrow) !g m-3
	    	read(ifin),((cldopd(i,j,k,ihr,idate),i=1,ncol),j=1,nrow) !g m-3
	    enddo
	enddo

 200    continue
        close(ifin)

      enddo !date

c
c  Tally and Print Average Rainwater and Cloud Optical Depth Profile
c
      do k = 1,nlay
        tee = 0
        tee2 = 0
        tee3 = 0
        do idate = 1,ndate
          do i = 1,ncol
            do j = 1,nrow
              do ihr = 1,24
                tee = tee + ranwtr(i,j,k,ihr,idate)/24/ncol/nrow/ndate
                tee2 = tee2 + cldopd(i,j,k,ihr,idate)/24/ncol/nrow/ndate
                tee3 = tee3 + cldwtr(i,j,k,ihr,idate)/24/ncol/nrow/ndate
              enddo
            enddo
          enddo
        enddo
        print *,'Average Rain (Layer',k,') = ',tee
        print *,'Average CldOD (Layer',k,') = ',tee2
        print *,'Average CldWtr (Layer',k,') = ',tee3
      enddo

c
c  Output Data
c
c    Begin Layer Loop
      do ilay = 1,nlay

c
c      Open Diurnal Map Output Files
c
	lDiurn = 0
	if(lDiurn.eq.1) then

	do io = 1,5
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout(io)//fnout2(1:INDEX(fnout2,' ')-1)//'.diurn.'//clay(ilay)
         open (UNIT = ifout+io, FILE=filename, STATUS='unknown')
         write(ifout+io,900), 'HOUR','X','Y',fnout(io)
	 !print *
	 !print *,'Filename: ',filename
	enddo
	    
c	 Average diurnal trend for all days and write map data
	   do ihr = 1,24
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
c		    c = 0
c		    d = 0
 		    e = 0
		    do idate = 1,ndate
			a = a + (cldwtr(i,j,ilay,ihr,idate)+cldwtr(i,j,ilay,ihr+1,idate))/2/ndate
			b = b + (ranwtr(i,j,ilay,ihr,idate)+ranwtr(i,j,ilay,ihr+1,idate))/2/ndate
c			c = c + (snowtr(i,j,ilay,ihr,idate)+snowtr(i,j,ilay,ihr+1,idate))/2/ndate
c			d = d + (gplwtr(i,j,ilay,ihr,idate)+gplwtr(i,j,ilay,ihr+1,idate))/2/ndate
			e = e + (cldopd(i,j,ilay,ihr,idate)+cldopd(i,j,ilay,ihr+1,idate))/2/ndate
		    enddo
		    write (ifout+1,901), chr(ihr), i, j, a
		    write (ifout+2,901), chr(ihr), i, j, b
c		    write (ifout+3,901), chr(ihr), i, j, c
c		    write (ifout+4,901), chr(ihr), i, j, d
		    write (ifout+5,901), chr(ihr), i, j, e
		enddo  !j
	     enddo  !i
	    enddo  !ihr
	do io = 1,5
	  close(ifout+io)
	enddo
	endif

c
c      Open Daily Average Map Output Files
c
	lDaily = 0
	if (lDaily.eq.1) then

	do io = 1,5
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout(io)//fnout2(1:INDEX(fnout2,' ')-1)//'.daily.'//clay(ilay)
         open (UNIT = ifout+io, FILE=filename, STATUS='unknown')
         write(ifout+io,900), 'DATE','X','Y',fnout(io)
c	 print *
	 print *,'Filename: ',filename
	enddo
    
c	 Average diurnal trend for all days and write map data
	   do idate = 1,ndate
	     do i = 1,ncol
	    	do j = 1,nrow
		    a = 0
		    b = 0
c		    c = 0
c		    d = 0
		    e = 0
		    do ihr = 1,24
			a = a + (cldwtr(i,j,ilay,ihr,idate)+cldwtr(i,j,ilay,ihr+1,idate))/2/24
			b = b + (ranwtr(i,j,ilay,ihr,idate)+ranwtr(i,j,ilay,ihr+1,idate))/2/24
c			c = c + (snowtr(i,j,ilay,ihr,idate)+snowtr(i,j,ilay,ihr+1,idate))/2/24
c			d = d + (gplwtr(i,j,ilay,ihr,idate)+gplwtr(i,j,ilay,ihr+1,idate))/2/24
			e = e + (cldopd(i,j,ilay,ihr,idate)+cldopd(i,j,ilay,ihr+1,idate))/2/24
		    enddo
		    write (ifout+1,901), cdate(idate), i, j, a
		    write (ifout+2,901), cdate(idate), i, j, b
c		    write (ifout+3,901), cdate(idate), i, j, c
c		    write (ifout+4,901), cdate(idate), i, j, d
		    write (ifout+5,901), cdate(idate), i, j, e
		enddo  !j
	     enddo  !i
	    enddo  !idate
	do io =1,5
	   close(ifout+io)
	enddo
	endif

       enddo !ilay

c
c   Output explicit (every hour and day) rain fields
c
      lexplicit = 1
      if (lexplicit.eq.1) then
      	do idate = 1,ndate
         filename = outdir(1:INDEX(outdir,' ')-1)//fname(1:INDEX(fname,' ')-1)//
     &              '.rainwtr.'//cdate(idate)//'.01'
	  open (unit = 20, file=filename, status='unknown')
          write(20,900), 'HOUR','X','Y','Rain Rate'
          do ihr = 1,24
	    do icol = 1,ncol
		do irow = 1,nrow
		    write (20,901), chr(ihr), icol, irow, ranwtr(icol,irow,1,ihr,idate)
		enddo
	    enddo
	  enddo
	  close(20)
        enddo
      endif

c
c     Output explicit rain fields for specific sites
c

      lPIT = 1
      if (lPIT.eq.1) then
        filename = outdir(1:INDEX(outdir,' ')-1)//fname(1:INDEX(fname,' ')-1)//
     &             '.cldwtr.ATL.01'  !.PIT.
	open (unit = 20, file=filename, status='unknown')
        write(20,903), 'DATE','HOUR','X','Y','Rain Rate'
	  do idate = 1,ndate
          do ihr = 1,24
		icol = 58 !65
		irow = 29 !51
		    write (20,904), cdate(idate),chr(ihr), icol, irow, ranwtr(icol,irow,1,ihr,idate)
	  enddo
	  enddo
	close(20)
      endif
	  
      !Calculate coefcld and print	
      icol = 58 !65
      irow = 29 !51
      filename = outdir(1:INDEX(outdir,' ')-1)//fname(1:INDEX(fname,' ')-1)//
     &           '.coefcld.ATL.01'  !.PIT.
	  open (unit = 20, file=filename, status='unknown')
      write(20,903), 'DATE','HOUR','X','Z','COEFCLD'
      do idate = 1,ndate
        do ihr = 1,24
          do ilay = 1, nlay
            tau = cldopd(icol,irow,ilay,ihr,idate)
			if (tau.lt.5.) then
              energy(ilay) = 1
			  cloud(ilay) = 0
			else
		      energy(ilay) = (5-exp(-tau))/(4+0.42*tau)
			  cloud(ilay) = 1
		    endif
            if (energy(ilay).ne.1) then
				iabov = 0
				ctrns = energy(ilay)
				fcld = cloud(ilay)
			else
				iabov = 1
				ctrns = energy(1)
				fcld = cloud(1)
		    endif
			zenang = 0.0175 * 30
			if (iabov.eq.1) then
				cldrat = 1 + (1 - ctrns)*cos(zenang)
			else
				cldrat = 1.6*ctrns*cos(zenang)
			endif
			coefcld = 1 + fcld*(cldrat - 1)
		    write (20,904), cdate(idate),chr(ihr), icol, ilay, coefcld
          enddo
        enddo
      enddo

 300    format(1X,F10.7)

 900  format(A4,2x,A2,2x,A2,2x,A12)
 901  format(A3,3x,I2,4x,I2,4x,E10.4)
 902  format(A3,3x,I2,4x,I2,4x,E10.4)
 903  format(A4,2x,A4,2x,A2,2x,A2,2x,A12)
 904  format(A4,2x,A4,2x,I2,2x,I2,2x,E10.4)

      end do
      end
 
