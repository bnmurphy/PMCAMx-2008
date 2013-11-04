      program windrd
c######################################################################
c
c    Windrd opens PMCAMx wind input files and prints the data
c	to be read and plotted by Matlab
c	Output File Name:  wind.layer
c
c    Written by Ben Murphy - 6-19-2008
c
c#####################################################################
c
      implicit none

      integer ifin, ifout, ifout2, ifout3
      integer ndate, ncol, nrow, nlay, nhr, nyr
      parameter (ifin  = 20)
      parameter (ifout = 30)
      parameter (ifout2= 31)
      parameter (ifout3= 32)
      parameter (ndate = 31)
      parameter (ncol  = 97)
      parameter (nrow  = 90)
      parameter (nlay  = 14)
      parameter (nhr   = 25)
      parameter (nyr   = 20)

      character*3   	cdate(ndate), chr(24)
      character*2       clay(nlay)
      character*4       cyr(nyr), year
      real		hour
      character*1       lstagger
      character*198 	filename, indir, outdir, fname, fname2, fname3, fnout
      real 		uwind(ncol,nrow,nlay,nhr,ndate), vwind(ncol,nrow,nlay,nhr,ndate)
      real		dummy, tee, tee2, a, b
      integer           lexplicit, lDaily
      integer           ihr, idate, ilay, iyr
      integer           d, i, j, k, n

c
c  Specify Data Vectors
c
      data cyr /'1991','1992','1993','1994','1995','1996','1997','1998','1999','2000'
     &         ,'2051','2052','2053','2054','2055','2056','2057','2058','2059','2060'/
c      data cyr /'2001'/

      data cdate /'182','183','184','185','186','187','188','189','190','191','192',
     &            '193','194','195','196','197','198','199','200','201','202','203',
     &            '204','205','206','207','208','209','210','211','212'/
c      data cdate /'193','194','195','196','197','198','199','200','201',
c     &		  '202','203','204','205','206','207','208','209'/

      data chr /'01','02','03','04','05','06','07','08','09','10','11','12',
     &          '13','14','15','16','17','18','19','20','21','22','23','24'/
      data clay /'01','02','03','04','05','06','07','08','09','10','11','12','13','14'/

c  Define year
      do iyr = 1,nyr
      year = cyr(iyr)
      print *, 'Starting ',year

c
c  Set Directory Structure
c
c      indir  = '/home/BaseInput/July/wind/'
c      fname  = 'wind.'//year !2001'
c      fname2 = '.4rpos.36.14.mm5.ld.camx'
c      fnout  = '.Jul'//year
c      outdir = '/home/mday/processed_output/2001meteorology/'

      indir  = '/home/mday/PMCAMx/met_hd/CAMx_Input/' // year // '/'
      outdir = outdir
      fname  = 'camx.uv.'//year
      fname2 = '.updt'
      fnout  = '.Jul'//year

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
         open (UNIT=ifin, FILE=filename, FORM='UNFORMATTED', STATUS='OLD')
	 print *,'Filename: ',filename

c
c      Read Time-Variant Portion
c
	do ihr = 1,nhr
	    read (ifin, end=200), hour, d
	    do k = 1,nlay
	    	read(ifin), ((uwind(i,j,k,ihr,idate),i=1,ncol),j=1,nrow) !m/s
	    	read(ifin), ((vwind(i,j,k,ihr,idate),i=1,ncol),j=1,nrow) !m/s
	    enddo
	    read(ifin), dummy
	enddo

 200    continue
        close(ifin)

      enddo  !Next Date

c
c  Tally and Print Average Wind Profile
c
      do idate = 1,ndate
        tee = 0
        do k = 1,nlay
          do i = 1,ncol
            do j = 1,nrow
              do ihr = 1,24
                tee = tee + uwind(i,j,k,ihr,idate)/24/ncol/nrow/ndate
                tee2 = tee2 + vwind(i,j,k,ihr,idate)/24/ncol/nrow/ndate
              enddo
            enddo
          enddo
        enddo
        print *,'Average UWIND (Layer',d,') = ',tee
        print *,'Average VWIND (Layer',d,') = ',tee2
      enddo

c
c      Open Daily Average Map Output Files
c
        lDaily = 1
        if (lDaily.eq.1) then
        do ilay = 1,2 !nlay
         fname = outdir(1:INDEX(outdir,' ')-1)//'uwind'//fnout(1:INDEX(fnout,' ')-1)//'.daily.'//clay(ilay)
         open (UNIT = ifout, FILE=fname, STATUS='unknown')
         write(ifout,900), 'DATE','X','Y','uwind(m/s)'
         fname2 = outdir(1:INDEX(outdir,' ')-1)//'vwind'//fnout(1:INDEX(fnout,' ')-1)//'.daily.'//clay(ilay)
         open (UNIT = ifout2, FILE=fname2, STATUS='unknown')
         write(ifout2,900), 'DATE','X','Y','vwind(m/s)'
         fname3 = outdir(1:INDEX(outdir,' ')-1)//'uvwind'//fnout(1:INDEX(fnout,' ')-1)//'.daily.'//clay(ilay)
         open (UNIT = ifout3, FILE=fname3, STATUS='unknown')
         write(ifout3,902), 'DATE','X','Y','uwind(m/s)','vwind(m/s)'

         print *,'Filenames: ',fname,fname2,fname3

c        Average diurnal trend for all days and write map data
           do idate = 1,ndate
             do i = 1,ncol
                do j = 1,nrow
                    a = 0
                    b = 0
                    do ihr = 1,24
                        a = a + (uwind(i,j,ilay,ihr,idate)+uwind(i,j,ilay,ihr+1,idate))/2/24
                        b = b + (vwind(i,j,ilay,ihr,idate)+vwind(i,j,ilay,ihr+1,idate))/2/24
                    enddo
                    write (ifout,901), cdate(idate), i, j, a
                    write (ifout2,901), cdate(idate), i, j, b
                    write (ifout3,903), cdate(idate), i, j, a, b
                enddo  !j
             enddo  !i
           enddo  !idate
           close(ifout)
           close(ifout2)
           close(ifout3)
        enddo !ilay
        endif

c
c      Output explicit (every hour and day) files
c
       lexplicit = 0
       if (lexplicit.eq.1) then
        do ilay = 1,1
         do idate = 1,ndate
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &                  fnout(1:INDEX(fnout,' ')-1)//cdate(idate)//'.uwind.'//clay(ilay)
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'HOUR','X','Y','Wind(m/s)'
         !print *,'Filename: ',filename

            do ihr = 1,24
              do i = 1,ncol
                do j = 1,nrow
                    write (ifout,901), chr(ihr), i, j, uwind(i,j,ilay,ihr,idate)
                enddo
              enddo 
            enddo         
          close(ifout) 
         enddo  !ndate
        enddo   !nlay
       endif

 300    format(1X,F10.7)

 900  format(A4,2x,A2,2x,A2,4x,A12)
 901  format(A3,3x,I2,4x,I2,4x,E10.4)
!            F5.0,3x,I2,4x,I2,2x,F8.3)
 902  format(A4,2x,A2,2x,A2,4x,A12,A12)
 903  format(A3,3x,I2,4x,I2,4x,E10.4,2x,E10.4)
      end do
      end
