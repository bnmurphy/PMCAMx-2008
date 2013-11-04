      PROGRAM  after
c######################################################################
c   BNM
c     This program processes PMCAMx results
c     The parameters are hardcoded for the specific output file
c     Ben Murphy rewrote this post-processor 3-11-08
c       outputs Gas, PM 2.5, and NOT(PM 10) averaged over each day
c
c######################################################################

c
c Define constants and arrays
c

      parameter(ifin  = 10)
      parameter(ifout = 30)

      parameter(ncol = 97)
      parameter(nrow = 90)
      parameter(nlay = 14)     ! Number of layers of interest
      parameter(nspc = 459)    ! Number of input species (SAPRC=459,CB4=158)
      parameter(nhr  = 24)     ! Hours in each day
      parameter(mspc = 28)     ! Number of output species
      parameter(agg  = 9)      ! Number of Aggregate Species
      parameter(simdays = 17)  ! Num of Simulation Days
      parameter(gasfac = 1000)

      character*198 fname1
      character*198 outdir, indir
      character*190 fname
      character*198 scenario   !the name of the output folder
      character*198 scenario2  !name of input file

      character*4   name(10),note(60),mspec(10,nspc)
      character*2   chr(nhr)         !Hour
      character*3   cdate(simdays)   !Julian Date
      character*2   clay(nlay)
      character*10  spcname(nspc)    !Input species name string
      character*7   pname(mspc)      !Output species name vector
      character*7   aggname(agg)     !Aggregate species name vector

      real*4 maxh(mspc+agg)
      real*4 mconc(simdays, 31, ncol, nrow, mspc+agg)
      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 conh(nspc,nhr,ncol,nrow,nlay)  !Conc read from input file
      real*4 sconh(mspc+agg,ncol,nrow)      !Conc for output (avg, sum, etc)
      real*4 date1
      
      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate,col,row
      integer*4 ihr,msec,icut,line
      integer*4 iday,ix,iy,zlay

c      data cdate /'194','195','196','197','198','199','200','201','202'/

      data cdate /'193','194','195','196','197','198','199','200',
     &	 	  '201','202','203','204','205','206','207','208',
     &		  '209'/
     
      data chr /'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &          '21','22','23','24'/

      data clay /'01','02','03','04','05','06','07','08','09','10',
     &                 '11','12','13','14'/

      data pname /'NO','NO2','O3','ISOP','TERP','SESQ','HNO3','HO2H',
     &		  'SO2','SULF','NH3','CPO','COO','CBS','CAS',
     &		  'APO','AOO','ABS','AAS','POC','PEC','CRST','PH2O',
     &		  'PCL','NA','PNH4','PNO3','PSO4'/

      data aggname /'TotNH4','TotNO3','SOA','POA','TotOM',
     &            'cSOA','cPOA','TotOC','PM25'/

c
c Define Input/Output Locations
c
      indir = '/home/mday/Output/'
      outdir = '/home/mday/processed_output/'
      scenario = 'PMCAMx2008_biogenic50_p5_ht'   !output folder
      scenario2 = 'PMCAMx2008_biogenic50_p5_ht'      !input file

c
c Initialize Maximum Concentration Counter
c
      do i = 1,mspc
         maxh(i) = 0.0
      enddo

c
c   Begin Looping Over Layers
c    
      do ilay = 1,nlay

c Output Files
c
	print *,'\n','Opening Output Files...','\n\n'
	do io = 1,(mspc+agg)

	  if (io.le.mspc) then  	!Input Species
	     isblk = INDEX(pname(io),' ')
             if(isblk.eq.0) isblk = 5
             fname = outdir(1:INDEX(outdir,' ')-1)//
     &		   scenario(1:INDEX(scenario,' ')-1)//'/'//
     &		   pname(io)(1:isblk-1)//'.daily'//clay(ilay)
	  else				!Aggregate Species
	     isblk = INDEX(aggname(io-mspc),' ')
	     if (isblk.eq.0) isblk = 7
	     fname = outdir(1:INDEX(outdir,' ')-1)//
     &		   scenario(1:INDEX(scenario,' ')-1)//'/'//
     &		   aggname(io-mspc)(1:isblk-1)//'.daily'//clay(ilay)
	  end if
        
	  isblk = INDEX(fname,' ')
	  open(unit=(ifout+io),file=fname(1:isblk-1),status='UNKNOWN')
      
	end do
c
c Iterate Days (input files)
c
      do idate=1,simdays

c
c   Open Input File
c
	if (nlay.gt.1) then
        	fname1='/4rpos.baseE.'//cdate(idate)//'.'//
     &			scenario2(1:INDEX(scenario2,' ')-1)//'.avrg'//clay(ilay)
	else
                fname1='/4rpos.baseE.'//cdate(idate)//'.'//
     &              scenario2(1:INDEX(scenario2,' ')-1)//'.avrg'//clay(ilay)
	endif
        fname = indir(1:INDEX(indir,' ')-1)//
     &          scenario2(1:INDEX(scenario2,' ')-1)//
     &          fname1(1:INDEX(fname1,' ')-1)
        write(6,*)'INPUT FILE: ',fname
        open(unit=ifin, file=fname, form='UNFORMATTED', status='old')
 
c    
c Clear Temporary Variables
c
       do io = 1,mspc+agg
            do i = 1,ncol
                do j = 1,nrow
                    do k = 1,nhr
                        mconc(idate,k,i,j,io) = 0.0
                    end do
                end do 
            end do
        end do

c
c   Read Input Header
c
      read(ifin,End=990)name,note,ione,nspec,ibdate,btime,iedate,etime
      write(6,*)(name(i)(1:1),i=1,10)
      write(6,*)(note(i)(1:1),i=1,60)
      write(6,*),nspec,ibdate,btime,iedate,etime
      If(nspec .ne. nspc)stop'nspec is not equal to nspc'
      read(ifin,END=990)rdum,rdum,iutm,xorg,yorg,delx,dely,nx,ny,nz
     &                  ,idum,idum,rdum,rdum,rdum
      write(6,*)iutm,xorg,yorg,delx,dely
      write(6,*)nx,ny,nz
      read(ifin, END=990)izero,izero,nx,ny
      write(6,*)izero,nx,ny
      if(nx.ne.ncol)stop'nx is not equal to ncol'
      if(ny.ne.nrow)stop'ny is not equal to nrow'
      if(nz.ne.nlay)stop'nz is not equal to nlay'

      read(ifin,END=990)((mspec(i,j),i=1,10),j=1,nspec)
      do isp = 1,nspc
	write (spcname(isp),*), (mspec(n,isp)(1:1),n=1,10)
      end do

c
c   Read Conc'ns into memory and close input file
c
      print *,'\nGetting Hourly Conc Data From Day ',cdate(idate),'...\n'
      do ihr = 1,nhr
	 read(ifin)ibdate,btime,iedate,etime
	 write(6,*)'ihr = ',ihr-1
	 do isp = 1,nspec
c	  write(*,*)isp,'/459 finished'
	    do zlay = 1,1
c	      write(*,*)'zlay= ',zlay
	       read(ifin)ione,(mspec(n,isp),n=1,10),
     &		  ((conh(isp,ihr,i,j,zlay),i=1,nx),j=1,ny)
c		write(*,*)'read= ',ione
c	       write (spcname(isp),*), (mspec(n,isp)(1:1),n=1,10)
	    end do  !zlay
	 end do  !isp
      end do  !ihr

      close(ifin)  !Close input file

c
c   Initialize Output Species Matrix
c
	print *,'Initializing Output Concentration Arrays...'
	do m=1,mspc
	   do i=1,ncol
	      do j=1,nrow
		 sconh(m,i,j)=0.0
	      end do
	   end do
	end do

c
c   Sum up concentrations
c
      print *,'Summing Hourly Concentrations...'
      do isp = 1,nspec
	isund = INDEX(spcname(isp),'_')
	      
	if (isund.eq.0) then !Gas

c                  ***** Gases *****
	    do kspc = 1,mspc  !Loop through output species
	      if (spcname(isp)(1:INDEX(spcname(isp),' ')-1).eq.
     &		  pname(kspc)(1:INDEX(pname(kspc),' ')-1)) then
            	do khr = 1,nhr
                  do i = 1,nx
            	  do j = 1,ny
		  do k = 1,1
		    sconh(kspc,i,j) = sconh(kspc,i,j) 
     &			+ conh(isp,khr,i,j,k)*gasfac
		  end do
		  end do
		  end do
		end do
	      elseif (spcname(isp)(1:3).eq.pname(kspc)(1:3).and.
     &		       (spcname(isp)(1:3).eq.'CPO'.or.
     &		        spcname(isp)(1:3).eq.'COO'.or.
     & 			spcname(isp)(1:3).eq.'CBS'.or.
     &			spcname(isp)(1:3).eq.'CAS')    ) then
		do khr = 1,nhr
		  do i = 1,nx
		  do j = 1,ny
		  do k = 1,1
		    sconh(kspc,i,j) = sconh(kspc,i,j)
     &			+ conh(isp,khr,i,j,k)*gasfac
		  end do
		  end do
 		  end do
		end do
	      end if
	    end do
	      
	else   !Aerosol

c                 ***** Aerosols *****
c		Doesn't do PM 10 yet!!!
	    read (spcname(isp)(isund+1:isund+3),*), bin
	    if (bin.le.6) then  !PM 2.5
	       do kspc = 1,mspc !Loop through output species
		  if ( spcname(isp)(1:3).eq.pname(kspc)
     &		     .or.spcname(isp)(1:4).eq.pname(kspc)
     &		     .or.(spcname(isp)(1:2).eq.'NA'
     &			  .and.pname(kspc)(1:2).eq.'NA') ) then
		     do khr = 1,nhr
			do i = 1,nx
			do j = 1,ny
			do k = 1,1
			   sconh(kspc,i,j) = sconh(kspc,i,j)
     &		             + conh(isp,khr,i,j,k)
			end do
		     	end do
			end do
		     end do
		  end if
	       end do
	    end if  !PM 2.5
		
	end if  !Gas or Aerosols
      end do	    !isp

c
c  Calculate Aggregate Quantities
c    Species Index MUST match up with vector <pname>
c
	print *,'Calculating Aggregate Species Concentrations...'

	do i = 1,nx
	do j = 1,ny
	
c	TotNH4 = PNH4 + NH3
      	sconh(29,i,j) = sconh(26,i,j)+sconh(11,i,j)*(1.0133*10**2*17/(8.314*300))
c	TotNO3 = PNO3 + HNO3	
	sconh(30,i,j) = sconh(27,i,j)+sconh(7,i,j)*(1.0133*10**2*63/(8.314*300))
c	CIS (currently refers to CPO?)
	  sconh(12,i,j) = sconh(12,i,j)*(1.0133*10**2*250/(8.314*300))
c	CR1 (currently refers to COO?)
          sconh(13,i,j) = sconh(13,i,j)*(1.0133*10**2*250/(8.314*300))
c       SOA = ABS + AAS        
          sconh(31,i,j) = sconh(18,i,j) + sconh(19,i,j)
c	POA = APO + AOO + POC      
	  sconh(32,i,j) = sconh(16,i,j) + sconh(17,i,j) + sconh(20,i,j)
c	TotOM = SOA + POA
	  sconh(33,i,j) = sconh(31,i,j) + sconh(32,i,j)

c       cSOA = ABS + AAS with OM->OC conversion factors      
          sconh(34,i,j) = sconh(18,i,j)/2 + sconh(19,i,j)/2
c       cPOA = APO + AOO + POC with OM->OC conversion factors      
          sconh(35,i,j) = sconh(16,i,j)/1.4 + sconh(17,i,j)/2 + sconh(20,i,j)/2
c       TotOC = cSOA + cPOA
          sconh(36,i,j) = sconh(34,i,j) + sconh(35,i,j)

c	PM2.5 = OM + PEC/2 + CRST + PCL + NA + ~PH2O + PNH4 + PNO3 + PSO4        
	  sconh(37,i,j) = sconh(33,i,j) + sconh(21,i,j)/2 + sconh(22,i,j)
     &                  + sconh(24,i,j) + sconh(25,i,j)
     &		        + sconh(26,i,j) + sconh(27,i,j) + sconh(28,i,j)

	end do
	end do

c
c  Write to Output Files
c
      print *,'Exporting Data to directory: ',outdir(1:INDEX(outdir,' ')-1)
      do io = 1,(mspc+agg)
        
	do j=1,ny     ! ROW LATER
          do i=1,nx   ! COLUMN FIRST
            write( (ifout+io), *),cdate(idate),j,i,sconh(io,i,j)/nhr
            maxh(io) = AMAX1(maxh(io),sconh(io,i,j))
          end do
        end do
c
      end do  !io

c
c  Done with Day -> Start Next One
c
	print *,'Done with Day ',cdate(idate),'\n\n'
      end do !Day -> idate iterates

c
c  Close Output Files
c
	print *,'Closing Output Files...','\n'
	do io = 1,(mspc+agg)
		close ((io+ifout))	!Close Output Files
	end do

      end do !lay -> iterate to next layer


c 41   format(10A,60A,3I10,F17.7,I10,F17.7)
c 42   format(2F17.7,I10,4F17.7,5I10,3F17.7)
c 43   format(4I10)
c 44   format(1580A)
c 45   format(I10,F17.7,I10,F17.7)
c 46   format(I10,10A,8730F17.7)
  47   format(i3,1x,i2,1x,i2,1x,f10.5)

       write(6,*)'ALL DAYS DONE'

 990   end
