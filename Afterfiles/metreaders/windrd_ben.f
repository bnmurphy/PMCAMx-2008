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

      integer ifin, ifout
      integer ndate, ncol, nrow, nlay

      parameter (ifin  = 20)
      parameter (ifout = 30)
      parameter (ndate = 17)
      parameter (ncol  = 97)
      parameter (nrow  = 90)
      parameter (nlay  = 14)

      character*3   	dates(ndate)
      real		hour
      integer		idate
      character*1       lstagger
      character*198 	filename, indir, outdir, fname, fname2, fnout
      real 		uwind(ncol,nrow,nlay,25,ndate), vwind(ncol,nrow,nlay,25,ndate)
      real		dummy, tee, tee2

      integer d, i, ihr, n, k, j

c
c  Specify Date Vectors
c
      data dates /'193','194','195','196','197','198','199','200','201',
     &		  '202','203','204','205','206','207','208','209'/
c      data dates /'193','194'/
c        data dates /'356','357','358','359','360','361','362','363',
c     &              '364','365'/

c
c  Set Directory Structure
c
c      indir  = '/home2/bnmurphy/PMCAMx_Input/Base_Input/2003/July/wind/'
c      fname  = 'camx.uv.2003'
c      fname2 = '.updt'
c      fnout  = 'wind.'
c      outdir = '/home2/bnmurphy/PMCAMx_Input/mm5camx_Ben_IO/confirm/out/'

      indir  = '/usr/BaseInput/July/wind/'
      fname  = 'wind.2001'
      fname2 = '.4rpos.36.14.mm5.ld.camx'
      fnout  = 'wind.'
      outdir = '/home/mday/processed_output/PMCAMx2008_BASE_sulf'

c
c  Begin Date Loop
c
      do d=1,ndate
	print *,'Day: ',dates(d)

c
c      Open Input File
c
         filename = indir(1:INDEX(indir,' ')-1)// fname(1:INDEX(fname,' ')-1)//
     &		   dates(d)//fname2(1:INDEX(fname2,' ')-1)
         open (UNIT=ifin, FILE=filename, FORM='UNFORMATTED', STATUS='OLD')

c
c      Open Output File
c
         filename = outdir(1:INDEX(outdir,' ')-1)//
     &			fnout(1:INDEX(fnout,' ')-1)//dates(d)//'.layer'
         open (UNIT = ifout, FILE=filename, STATUS='unknown')
         write(ifout,900), 'HOUR','X','Y',(n,n=1,14)
	 !print *
	 !print *,'Filename: ',filename

c
c      Read Time-Variant Portion
c
	do ihr = 1,25
	    read (ifin, end=200), hour, idate
	    do k = 1,nlay
	    	read(ifin), ((uwind(i,j,k,ihr,d),i=1,ncol),j=1,nrow)
	    	read(ifin), ((vwind(i,j,k,ihr,d),i=1,ncol),j=1,nrow)
	    enddo
	    read(ifin), dummy

	    do i = 1,ncol
	    	do j = 1,nrow
		    write (ifout,901), hour, i, j, (uwind(i,j,k,ihr,d),k=1,14)
		enddo
	    enddo

	enddo

 200    continue
        close(ifin)
	close(ifout)

      enddo  !Next Date

c
c  Tally and Print Average Wind Profile
c
      do d = 1,ndate
        tee = 0
        do k = 1,nlay
          do i = 1,ncol
            do j = 1,nrow
              do ihr = 1,24
                tee = tee + uwind(i,j,k,ihr,d)/24/ncol/nrow/ndate
                tee2 = tee2 + vwind(i,j,k,ihr,d)/24/ncol/nrow/ndate
              enddo
            enddo
          enddo
        enddo
        print *,'Average UWIND (Layer',d,') = ',tee
        print *,'Average VWIND (Layer',d,') = ',tee2
      enddo


 300    format(1X,F10.7)

 900  format(A4,2x,A2,2x,A2,14(8x,I2))
 901  format(F5.0,3x,I2,4x,I2,14(2x,F8.3))

      end
 
