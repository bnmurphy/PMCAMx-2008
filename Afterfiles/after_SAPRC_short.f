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
c      implicit none

      parameter(ifin  = 10)
      parameter(ifout = 30)

      parameter(ncol = 97)
      parameter(nrow = 90)
      parameter(nlay = 1)
      parameter(nspc = 459)   ! Number of input species (459)
      parameter(nhr  = 24)    ! Hours in each day
      parameter(mspc = 78)    ! Number of output species
      parameter(agg  = 9)     ! Number of Aggregate Species
      parameter(simdays = 12)  ! Num of Simulation Days
      parameter(gasfac = 1000)

      character*198 fname1
      character*198 outdir, indir
      character*190 fname
      character*198 scenario   !the name of the output folder
      character*198 scenario2  !name of input file

      character*4   name(10),note(60),mspec(10,nspc)
      character*2   chr(nhr)         !Hour
      character*3   cdate(simdays)   !Julian Date
      character*10  spcname(nspc)    !Input species name string
      character*7   pname(mspc)      !Output species name vector
      character*7   aggname(agg)     !Aggregate species name vector

      real*4 maxh(mspc+agg)
      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 conh(nspc,nhr,ncol,nrow,nlay)  !Conc read from input file
      real*4 sconh(mspc+agg,ncol,nrow)      !Conc for output (avg, sum, etc)
      real*4 date1
      
      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate,col,row
      integer*4 ihr,msec,icut,line
      integer*4 iday,ix,iy

      data cdate /'182','183','184','185','186','187','188','189','190','191','192','193'/
c     &            '194','195','196','197','198','199','200',
c     &	 	  '201','202','203','204','205','206','207','208',
c     &		  '209'/
     
      data chr /'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &          '21','22','23','24'/


      data pname /'NO','NO2','O3','PAN','PAN2','MPAN','PBZN','NPHE',
     &		  'RNO3','CRES','DCB2','DCB3','HNO4','BALD','HONO',
     &		  'XN','HCHO','CCHO','RCHO','BACL','PROD','DCB1',
     &		  'PHEN','ISOP','ISPD','MVK','METH','MGLY','GLY',
     &		  'TERP','BPIN','LIMO','MONO','SESQ','HNO3','HO2H',
     &		  'HC2H','CO2H','CO3H','RC2H','RC3H','ACET','MEK',
     &		  'MEOH','COOH','ROOH','CO','ETHE','ALK1','ALK2',
     &		  'ALK3','ALK4','ALK5','ARO1','ARO2','OLE1','OLE2',
     &		  'NXOY','SO2','SULF','NH3','CPO','COO','CBS','CAS',
     &		  'APO','AOO','ABS','AAS','POC','PEC','CRST','PH2O',
     &		  'PCL','NA','PNH4','PNO3','PSO4'/

      data aggname /'TotNH4','TotNO3','SOA','POA','TotOM',
     &            'cSOA','cPOA','TotOC','PM25'/

c
c Define Input/Output Locations
c
      indir = '/home/mday/mday/met_hd/Output/'
      outdir = '/home/mday/mday/met_hd/processed_output/'
      scenario = 'PMCAMx2008_MGN_1994'   !output folder
      scenario2 = 'PMCAMx2008_MGN3_1994'  !input folder

c
c Initialize Maximum Concentration Counter
c
      do i = 1,mspc
         maxh(i) = 0.0
      enddo

c
c Open Output Files
c
	print *,'\n','Opening Output Files...','\n\n'
	do io = 1,(mspc+agg)

	  if (io.le.mspc) then  	!Input Species
	     isblk = INDEX(pname(io),' ')
             if(isblk.eq.0) isblk = 5
             fname = outdir(1:INDEX(outdir,' ')-1)//
     &		   scenario(1:INDEX(scenario,' ')-1)//'/'//
     &		   pname(io)(1:isblk-1)//'.daily'
	  else				!Aggregate Species
	     isblk = INDEX(aggname(io-mspc),' ')
	     if (isblk.eq.0) isblk = 7
	     fname = outdir(1:INDEX(outdir,' ')-1)//
     &		   scenario(1:INDEX(scenario,' ')-1)//'/'//
     &		   aggname(io-mspc)(1:isblk-1)//'.daily'
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
        fname1='/4rpos.baseE.'//cdate(idate)//'.'//
     &		scenario2(1:INDEX(scenario2,' ')-1)//'.avrg01'
        fname = indir(1:INDEX(indir,' ')-1)//
     &          scenario2(1:INDEX(scenario2,' ')-1)//
     &          fname1(1:INDEX(fname1,' ')-1)
        write(6,*)'INPUT FILE: ',fname
        open(unit=ifin, file=fname, form='UNFORMATTED', status='old')
     
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
      do ihr = 1,24
	 read(ifin)ibdate,btime,iedate,etime
	 write(6,*)'ihr = ',ihr-1
	 do isp = 1,nspec
	    do ilay = 1,nz
	       read(ifin,END=990)ione,(mspec(n,isp),n=1,10),
     &		  ((conh(isp,ihr,i,j,ilay),i=1,nx),j=1,ny)
c	       write (spcname(isp),*), (mspec(n,isp)(1:1),n=1,10)
	    end do  !ilay
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
		    sconh(kspc,i,j) = sconh(kspc,i,j) 
     &			+ conh(isp,khr,i,j,1)*gasfac
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
		    sconh(kspc,i,j) = sconh(kspc,i,j)
     &			+ conh(isp,khr,i,j,1)*gasfac
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
			   sconh(kspc,i,j) = sconh(kspc,i,j)
     &		             + conh(isp,khr,i,j,1)
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
      	sconh(79,i,j) = sconh(76,i,j)+sconh(61,i,j)*(1.0133*10**2*17/(8.314*300))
c	TotNO3 = PNO3 + HNO3	
	sconh(80,i,j) = sconh(77,i,j)+sconh(35,i,j)*(1.0133*10**2*63/(8.314*300))
c	CIS (currently refers to CPO?)
	  sconh(62,i,j) = sconh(62,i,j)*(1.0133*10**2*250/(8.314*300))
c	CR1 (currently refers to COO?)
          sconh(63,i,j) = sconh(63,i,j)*(1.0133*10**2*250/(8.314*300))
c       SOA = ABS + AAS        
          sconh(81,i,j) = sconh(68,i,j) + sconh(69,i,j)
c	POA = APO + AOO + POC      
	  sconh(82,i,j) = sconh(66,i,j) + sconh(67,i,j) + sconh(70,i,j)
c	TotOM = SOA + POA
	  sconh(83,i,j) = sconh(81,i,j) + sconh(82,i,j)

c       cSOA = ABS + AAS with OM->OC conversion factors      
          sconh(84,i,j) = sconh(68,i,j)/2 + sconh(69,i,j)/2
c       cPOA = APO + AOO + POC with OM->OC conversion factors      
          sconh(85,i,j) = sconh(66,i,j)/1.4 + sconh(67,i,j)/2 + sconh(70,i,j)/2
c       TotOC = cSOA + cPOA
          sconh(86,i,j) = sconh(84,i,j) + sconh(85,i,j)

c	PM2.5 = OM + PEC/2 + CRST + PCL + NA + ~PH2O + PNH4 + PNO3 + PSO4        
	  sconh(87,i,j) = sconh(83,i,j) + sconh(71,i,j)/2 + sconh(72,i,j)
     &                  + sconh(74,i,j) + sconh(75,i,j)
     &		        + sconh(76,i,j) + sconh(77,i,j) + sconh(78,i,j)

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

c 41   format(10A,60A,3I10,F17.7,I10,F17.7)
c 42   format(2F17.7,I10,4F17.7,5I10,3F17.7)
c 43   format(4I10)
c 44   format(1580A)
c 45   format(I10,F17.7,I10,F17.7)
c 46   format(I10,10A,8730F17.7)
  47   format(i3,1x,i2,1x,i2,1x,f10.5)

       write(6,*)'ALL DAYS DONE'

 990   end
