      program pointread
c######################################################################
c
c    Pointread opens PMCAMx point emissions files based on 2002 NEI
c	inventory compatible with SAPRC, reads in all information,
c	and writes it back out in a human readable format
c
c    Written by Ben Murphy - 6-9-2008
c
c#####################################################################
c
      parameter (ifin  = 20)
      parameter (ifout = 30)

c      integer begtim=0000
      integer begdate
c      integer endtim=2400
      integer enddate
      integer dxmod
      integer dymod
      integer d, io, p(500000)
      integer k,n,l,i,a,t,j,m
      parameter (ndate=17)
      parameter (ncol = 97)
      parameter (nrow = 90)

      character*3   	dates(ndate)
      character*198 	filename, indir, outdir, fname, fname2
      character*4   	ptspec1(10,200)
      character*4   	name1(10)
      character*10  	name2(1000)
      real 		ptemis1(120000,500)

      character*4 ifile1(10),note1(60)
      integer 	  nseg1,nptspc1,idat11,idat21,izone1,nx1,ny1,nz1,idum1
      real 	  tim11,tim21,orgx1,orgy1,utmx1,utmy1,dx1,dy1
      integer 	  npts1
      real 	  flowrat1(120000)
      real 	  effph1(120000),temp(120000)
      real 	  xloc1(120000),yloc1(120000)
      real 	  hstk1(120000),dstk1(120000)
      real 	  tstk1(120000),vstk1(120000)
      real 	  xstk(120000,1),ystk(120000,1)
      integer 	  map(120000,3), counter, loc
c
      begdate=274
      enddate=304
      dxmod=97
      dymod=90
      counter=0

c
c  Specify Date and Comparison Vectors
c
      data dates /'193','194','195','196','197','198','199','200','201',
     &		  '202','203','204','205','206','207','208','209'/

c
c  Set Directory Structure
c

      indir = '/home/BaseInput/July/SAPRC/prep99/semivol/'
      outdir = '/home/mday/processed_output/sites/'
      fname = 'point.camx.9902.2001'
      fname2 = '.bin'

c
c  Begin Date Loop
c
      do d=1,ndate

c
c      Open Input and Output Files
c
         filename = indir(1:INDEX(indir,' ')-1)// fname(1:INDEX(fname,' ')-1)//
     &		   dates(d)//fname2(1:INDEX(fname2,' ')-1)
         open (UNIT=ifin, FILE=filename, FORM='UNFORMATTED', STATUS='OLD')

         filename = outdir(1:INDEX(outdir,' ')-1)// fname(1:INDEX(fname,' ')-1)//
     &             dates(d)//fname2(1:INDEX(fname2,' ')-1)//'.txt'
         open (UNIT = ifout, FILE=filename, STATUS='unknown')

         write(6,*) 'Date: ', dates(d)

c      Header Line 1 (name,note,ione,nspec,ibdate,btime,iedate,etime)
         read(ifin) ifile1,note1,nseg1,nptspc1,idat11,tim11,idat21,tim21
	 write(6,*), (ifile1(i)(1:1),i=1,10),':  ',(note1(i)(1:1),i=1,60)
	 write(6,*), 'Number of Species: ',nptspc1
	 write(6, fmt='(A16,I4,2x,A16,E8.3)'), 'Beginning Date: ',idat11,'Beginning Time: ',tim11
	 write(6,fmt='(A10,I4,2x,A10,E8.3)'), 'End Date: ',idat21,'End Time: ',tim21

c      Header Line 2 (rdum,rdum,iutm,xorg,yorg,delx,dely,nx,ny,nz,idum,idum,rdum,rdum,rdum)
         read(ifin) orgx1,orgy1,izone1,utmx1,utmy1,dx1,dy1,nx1,ny1,nz1
c	 write(ifout,fmt='(A6,E9.3,2x,A6,E9.3,2x,A10,I3,1x,I3,1x,I2)'), 
c     &			'delx:',dx1,'dely:',dy1, 'nx,ny,nz: ',nx1,ny1,nz1

c      Header Line 3 (ione,ione,nx,ny)
         read(ifin) (idum1,idum1,idum1,idum1,n=1,nseg1)
c	 write(ifout,fmt='(I4,2x,I4)'), nx1, ny1

c      Header Line 4 (mspec(l),l=1,nspec)
         read(ifin) ((ptspec1(n,l),n=1,10),l=1,nptspc1)
c	 do i = 1,nptspc1
c		write(ifout,fmt='(A13,I3,1x,A4,2x,12A1)'), 
c     &			'Species No:',i,'is:',(ptspec1(n,i)(1:1),n=1,10)
c	 end do

c      Header Line 5 (ione, nstk)
         read(ifin) idum1,nptsrc1
c	 write(ifout,*)
c	 write(ifout,fmt='(A26,I6)'),'Number of Point Sources: ',nptsrc1

c      Header Line 6 (xstk, ystk, stack height(m), stack diam(m), 
c		stack Temp(K), stack vel(m/hr) )
         read(ifin) (xloc1(n),yloc1(n),hstk1(n),dstk1(n),tstk1(n),vstk1(n),
     &             n=1,nptsrc1)
c	 write( ifout, fmt='(A9,2x,2(A4,6x),A6,3(6x,A4))' ),
c     &		'PtSource','xloc','yloc','height','diam','temp','vel'
c	 do n = 1,nptsrc1
c		write(ifout,fmt='(I6,5x,E9.3,2x,E9.3,2x,E9.3,2x,E9.3,2x,E9.3,2x,E9.3)' ),
c     &			n,xloc1(n),yloc1(n),hstk1(n),dstk1(n),tstk1(n),vstk1(n)
c	 end do

c
c      Map pt sources to grid
c      
c       write(ifout,*)

       counter = 0
 	write(ifout,*)'Stack # ',' x ',' y ','Plume (m) '
       do n=1,nptsrc1
         xstk(n,1)=xloc1(n)/1000.+900.
         ystk(n,1)=yloc1(n)/1000.+1620.
         map(n,1)=1
         map(n,2)=1+INT(xstk(n,1)/36)
         map(n,3)=1+INT(ystk(n,1)/36)
         if(map(n,2).gt.ncol .or. map(n,3).gt.nrow .or. 
     &	    map(n,2).lt.0    .or. map(n,3).lt.0) then
           map(n,1)=0
	 else
          write(ifout,*)n,map(n,2),map(n,3),hstk1(n)
	 end if
       end do

c         t = 0
c 100     read(ifin,end=200) idat11,tim11,idat21,tim21
c         t = t + 1
cc       Read (ione, number of stacks)
c         read(ifin) idum1,npts1
c         read(ifin) (idum1,idum1,idum1,flowrat1(m),effph1(m),m=1,npts1)
c
c	   write(ifout,fmt='(A10,I6,A7,2x,I2,1x,I2)'), 'Stack No. ',n,' is at',map(n,2),map(n,3)
c	   write(ifout,*)n,map(n,2),map(n,3),effph1(n)
c           counter = counter + 1
c           p(counter) = n
c         endif

c        goto 100

c       enddo

c       write(ifout,*)
c       write(ifout,*), counter,' stacks are in domain.'

c
c      Write Location Representation
c
c	write(ifout,*)
c	write(ifout,*),'Locations'
c	write(ifout,fmt='(4A12)'),'col','row','site id'
c	do i = 1,ncol
c	   do j = 1,nrow
c		do n = 1,nptsrc1
c		   if (map(n,2).eq.i .and. map(n,3).eq.j) then
c			write(ifout,*), i,j,n
c		   end if
c		end do
c	   end do
c	end do

c
c      Read Time-Variant Portion
c

c       Time (ibdate,btime,iedate,etime)
c	 t = 0
c
c 100     read(ifin,end=200) idat11,tim11,idat21,tim21
c	 write(ifout,*)
c	 write(ifout,fmt='(A15,I5,2x,E9.3,2x,A15,I5,2x,E9.3)'), 
c     &		'Beg date/time: ',idat11,tim11,'End date/time: ',idat21,tim21
c
c	 t = t + 1
c	Read (ione, number of stacks)
c         read(ifin) idum1,npts1
c         read (ifin) (idum1,idum1,idum1,flowrat1(n),effph1(n),n=1,npts1)
c	 write(ifout,*)
c	 write(ifout,*), 'Effective Plume Heights'
c	 write(ifout,fmt='(A20,A20,A20)'), 'Stack No.','Flowrate','Plume Height'
c	 write(ifout,fmt='(6x,I8,6x,6x,E8.2,6x,5x,E10.4)'), (p(n), flowrat1(p(n)),effph1(p(n)),n=1,counter)
c	 write(ifout,*)
c
c       Species Concentrations
c
c        do ll = 1,nptspc1
c            read(ifin) idum1,(name1(i),i=1,10),(ptemis1(n,ll),n=1,npts1)
c	    write (name2(ll),*), (name1(i)(1:1),i=1,10)
c	    write(ifout,*)
c	    write(ifout,*),'Species: ',name2(ll)
c	    write(ifout,fmt='(4(4x,A5,I6,A2,E9.3))'),
c     &			('Stack',p(n),': ',ptemis1(p(n),ll),n=1,counter)
c        enddo
c
c        goto 100

 200    continue
        close(ifin)

      enddo  !Date

 300    format(1X,F10.7)

 900  format(A4,2x,A4,2x,A2,2x,A2,2x,A4,2x,A14)
 901  format(A3,3x,I2,4x,I2,2x,I2,2x,I4,2x,E11.3)

      end
 
