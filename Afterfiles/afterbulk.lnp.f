      program after

c  This program processes PMCAMx-UF .AVRG files to create hourly files for species.
c  NOTES:
c  io -  a counter within the specified array, ie io = 2
c  for cn_name is the second CN in the array (CN10).

      parameter(ifin = 10)
      parameter(ifout = 30)
      parameter(ncol = 97) ! Number of columns in simulated domain
      parameter(nrow = 90) ! Number of rows in simulated domain
      parameter(nlay = 1) ! Number of vertical layers simulated (1 if no 3D average output) 
      parameter(nspc = 615) ! Number of avg species
      parameter(nhr  = 24) ! Number of hours in simulation day
      parameter(mspc = 22) ! Number of species to be processed
      parameter(simdays = 15) ! Number of simulation days to be processed
      parameter(nsect = 43) ! Number of size bins
      parameter(trsp = 15) ! Number of tracting species
      parameter(cnsp = 5) ! Number of types of CN's

      character*99 fname1
      character*99 fname2
      character*99 fdir1
      character*99 fdir2
      character*198 fname
      character*198 fdir
      character*198 fnamex
      character*4 name(10),note(60),mspec(10,nspc)
      character*10 spcname(nspc)
      character*10 pname(mspc)
      character*3 cdate(simdays)
      character*2 chr(nhr)
      character*10 pname2(trsp)
      character*10 cn_name(cnsp) 

      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 conh(nspc,ncol,nrow,nlay) ! data read from ascii file
      real*4 sconh(mspc,ncol,nrow,nlay) ! saving variables
      real*4 cn(cnsp,ncol,nrow,nlay)

      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate
      integer*4 ihr,msec,icut
      integer*4 iday,ix,iy

c SPECIFY OUTPUT FILE DIRECTORY 
cexample      fdir1 = '/home/lnp/processedoutput/NewEmissNoSO4Emiss/'
cexample      fdir2 = 'Processed_NoSO4/After_Bulk/All_Locations/'

      fdir1 = '/home/lnp/processedoutput/Base_NoNatGas_NucON_SA.07312012/'
      fdir2 = 'Processed_Base/After_Bulk/All_Locations/'

      isblk2=INDEX(fdir1,' ')
      fdir=fdir1(1:isblk2-1)//fdir2

c SPECIFY SIMULATED DAYS TO PROCESS
      data cdate /'195', '196', '197', '198', '199', '200', '201', '202',
     &            '203', '204', '205', '206', '207', '208', '209'/

c SPECIFY SPECIES FOR CALCULATIONS
      data pname /'PALL','PTOC','SOA ','PEC ','PNH4','PNO3','PSO4',
     &         'CRST','SO2 ','NH3','SULF','NUM','NO2','NO ',
     &         'O3','HNO3','FORM','PAN','ETH','CO','TOL','XYL'/

      data pname2 /'PNH4','PNO3','PSO4','EC ','OC ','CRST','PCL ',
     &         'NA ','NUM','PALL','POC','SOA1','SOA2','SOA3','SOA4'/

      data cn_name /'CNTOT','CN3','CN10','CN50','CN100'/

c
c     Set chr variables
c
      do i=1, nhr
         if (i.le.9) then
            write(chr(i),'(a1,i1)')'0',i
         else
            write(chr(i),'(i2)')i
         endif
      enddo
c
c Loop for PM 2.5 and PM 10
c
      do icut=1,2 
        if (icut.eq.1) then
          msec = 35
          write(6,*)'PM 2.5 and gases'
        else
          msec = 41
          write(6,*)'PM 10'
        endif
c
c initialization
c
c iterate days
c
      idate=0
 99   idate=idate+1
      if(idate.gt.simdays) goto 1099

c SPECIFY INPUT FILE NAMES
cexample      fname1='/home/lnp/Output/NewEmissNoSO4Below205nm.092111/'
cexample      fname2='4rpos.baseE.'//cdate(idate)//
cexample     &        '.NewEmissNoSO4Below205nm.092111.avrg'

      fname1='/home/lnp/output/Base_NoNatGas_NucON_SA.07312012/'
      fname2='4rpos.baseE.'//cdate(idate)//
     &        '.Base_NoNatGas_NucON_SA.07312012.avrg'

      isblk=INDEX(fname1,' ')
      fname=fname1(1:isblk-1)//fname2
      isblk=INDEX(fname,' ')
      write(6,*)'INPUT FILE: ',fname(1:isblk-1)

c open input file
      open(unit=ifin,file=fname(1:isblk-1),form='UNFORMATTED',
     &     status='old')
c
c read header portion
c
      read(ifin,END=990)name,note,ione,nspec,ibdate,btime,iedate,etime
      write(6,*)(name(i)(1:1),i=1,10)
      write(6,*)(note(i)(1:1),i=1,60)
      write(6,*)ione,nspec,ibdate,btime,iedate,etime
      if(nspec.ne.nspc)stop'nspec is not equal to nspc'
      read(ifin,END=990)rdum,rdum,iutm,xorg,yorg,delx,dely,nx,ny,nz
     &                 ,idum,idum,rdum,rdum,rdum
      write(6,*)iutm,xorg,yorg,delx,dely
      write(6,*)nx,ny,nz
      read(ifin,END=990)izero,izero,nx,ny
      write(6,*)izero,nx,ny
      if(nx.ne.ncol)stop'nx is not equal to ncol'
      if(ny.ne.nrow)stop'ny is not equal to nrow'
      if(nz.ne.nlay)stop'nz is not equal to nlay'
      read(ifin,END=990)((mspec(i,j),i=1,10),j=1,nspec)
c
c read hourly portion
c
c
      do khr=1, nhr
      read(ifin,END=990)ibdate,btime,iedate,etime
      ihr = NINT(btime) +1
      write(6,*)'ihr = ',ihr
      if(ihr.gt.nhr)stop'incorrect hr index'
      if(ihr.ne.khr)stop'ihr ne khr'
      do isp=1,nspec
        do k=1,nz
        read(ifin,END=990)ione,(mspec(n,isp),n=1,10),
     &                    ((conh(isp,i,j,k),i=1,nx),j=1,ny)
        enddo
        write(spcname(isp),*)(mspec(n,isp)(1:1),n=1,10)
      enddo
c
c calculate desired quantities
c
c initialize conc variables
c
      do io=1,mspc
        do i=1,nx
          do j=1,ny
            do k=1,nz
              sconh(io,i,j,k)=0.0
            enddo
          enddo
        enddo
      enddo

      do io=1,cnsp
        do i=1,nx
          do j=1,ny
            do k=1,nz
              cn(io,i,j,k)=0.0
            enddo
          enddo
        enddo
      enddo

      do isp=1,nspec
        isund=INDEX(spcname(isp),'_')
        if(isund.gt.0) then ! AEROSOL
          if(spcname(isp)(isund+1:isund+2).eq.'1 ') then 
            !Start only from the 1st size section
            if((spcname(isp)(1:4).ne.'PH2O').and.
     &        (spcname(isp)(1:3).ne.'NUM')) then
              do i=1,nx
                do j=1,ny
                  do k=1,nz
                    do mm=1,msec
                      sconh(1,i,j,k) = sconh(1,i,j,k) 
     &                  + conh(isp-1+mm,i,j,k)
                    enddo
                  enddo ! k
                enddo ! j
              enddo ! i
            endif
            if((spcname(isp)(1:3).eq.'SOA').or.
     &        (spcname(isp)(1:3).eq.'POC')) then
              do i=1,nx
                do j=1,ny
                  do k=1,nz
                    do mm=1,msec
                      sconh(2,i,j,k) = sconh(2,i,j,k)
     &                  + conh(isp-1+mm,i,j,k)
                    enddo
                  enddo ! k
                enddo ! j
              enddo ! i
            endif
            if(spcname(isp)(1:3).eq.'SOA') then
              do i=1,nx
                do j=1,ny
                  do k=1,nz
                    do mm=1,msec
                      sconh(3,i,j,k) = sconh(3,i,j,k) 
     &                  + conh(isp-1+mm,i,j,k)
                    enddo
                  enddo ! k
                enddo ! j
              enddo ! i
            endif
            do kspc=4,8 ! Other aerosols JPD
              if(spcname(isp)(1:isund-1).eq.
     &          pname(kspc)(1:isund-1)) then
                do i=1,nx
                  do j=1,ny
                    do k=1,nz
                      do mm=1,msec
                        sconh(kspc,i,j,k) = sconh(kspc,i,j,k) + 
     &                    conh(isp-1+mm,i,j,k)
                      enddo
                    enddo ! k
                  enddo ! j
                enddo ! i
              endif
            enddo
            if(icut.eq.2) then ! Only for PM10
              if(spcname(isp)(1:3).eq.'NUM') then
                do i=1,nx
                  do j=1,ny
                    do k=1,nz
                      do mm=1,msec
                        sconh(12,i,j,k) = sconh(12,i,j,k) + 
     &                    conh(isp-1+mm,i,j,k)
                      enddo
                    enddo ! k
                  enddo ! j
                enddo ! i
              endif
            endif
          endif ! From 1st size bin
        endif ! AEROSOL
c
        if(icut.eq.1) then !Do only once when PM2.5 is calculated.
        if(isund.eq.0) then ! GAS
          do io=9,mspc
            if (io.eq.12) goto 300
            isblk=INDEX(pname(io),' ')
            if (isblk.eq.0) isblk = 5
            if(spcname(isp)(1:isblk).eq.pname(io)(1:isblk)) then 
              !To distinguish NO and NO2, isblk is used rather than isblk-1
              do i=1,nx
                do j=1,ny
                  do k=1,nz
                    sconh(io,i,j,k) = sconh(io,i,j,k) + 
     &                conh(isp,i,j,k)
                  enddo ! k
                enddo ! j
              enddo ! i
            endif
 300        continue
          enddo
        endif ! GAS
        endif ! icut = 1
c
c To get (CN)
c
        if(icut.eq.2) then !Do only once when PM10 is calculated.
          if(isund.gt.1) then ! AEROSOL
            if(spcname(isp)(isund+1:isund+2).eq.'1 ') then 
               !Start only from 1st size section
              if(spcname(isp)(1:3).eq.'NUM') then
                do i=1,nx
                  do j=1,ny
                    do k=1,nz
                      do mm=1,msec !cntot - all size sections (0.8nm-40um) 
                        cn(1,i,j,k) = cn(1,i,j,k) +
     &                   conh(isp-1+mm,i,j,k)
                      enddo
                      do mm=7,msec !cn3 - above 3 nm
                       cn(2,i,j,k) = cn(2,i,j,k) +
     &                   conh(isp-1+mm,i,j,k)
                      enddo
                      do mm=12,msec !cn10 - above 10 nm
                       cn(3,i,j,k) = cn(3,i,j,k) + 
     &                   conh(isp-1+mm,i,j,k)
                      enddo
                       do mm=19,msec !cn50 - above 50 nm
                        cn(4,i,j,k) = cn(4,i,j,k) +
     &                   conh(isp-1+mm,i,j,k)
                      enddo
                       do mm=22,msec !cn100 - above 100 nm
                        cn(5,i,j,k) = cn(5,i,j,k) +
     &                   conh(isp-1+mm,i,j,k)
                      enddo
                    enddo ! k
                  enddo ! j
                enddo ! i
              endif ! NUM
            endif ! 1st size section
          endif ! AEROSOL
        endif ! icut = 2

      enddo ! isp
c
c open output files for this hour
c
      if(icut.eq.1) then 
        do io=1,mspc 
          ! H2SO4, NH3, SO2, NO2, NO, O3, and HNO3, etc Gas   JGJ  5/22/2005
          if (io.eq.12) goto 10
          isblkx=INDEX(fdir,' ')
          isblk=INDEX(pname(io),' ')
          if(isblk.eq.0) isblk = 5
          if(io.le.8) then ! PM 2.5
            fnamex=fdir(1:isblkx-1)//pname(io)(1:isblk-1)//'25.'//
     &        cdate(idate)//'.'//chr(khr)
          endif
          isblk=INDEX(fnamex,' ')
          ifoutx = ifout + io
          open(unit=ifoutx,file=fnamex(1:isblk-1),status='unknown')
          do j=1,ny ! ROW 
            do i=1,nx ! COLUMN 
              do k=1,nz ! LAY 
                write(ifoutx,*)sconh(io,i,j,k)
cdbg                maxh(io) = AMAX1(maxh(io),sconh(io,i,j,k))
              enddo
            enddo
          enddo
          close(ifoutx)
 10       continue
        enddo
      endif ! icut = 1

      if(icut.eq.2) then 
        ! Only write CN when doing PM 10
        do io=1,cnsp
          isblkx=INDEX(fdir,' ')
          isblk=INDEX(cn_name(io),' ')
          fnamex=fdir(1:isblkx-1)//cn_name(io)(1:isblk-1)//'.'//
     &      cdate(idate)//'.'//chr(khr)
          isblk=INDEX(fnamex,' ')
          ifoutx = ifout + io
          open(unit=ifoutx,file=fnamex(1:isblk-1),status='unknown')
          do j=1,ny ! ROW 
            do i=1,nx ! COLUMN 
             do k=1,1 ! LAY 
                write(ifoutx,*)cn(io,i,j,k)
              enddo
            enddo
          enddo
          close(ifoutx)
        enddo
      endif ! icut = 2

      enddo ! khr
c
      goto 1000
c -------------------------------------------------------------------
 990  stop'END OF FILE: conc. average input file'

 1000 close(ifin)

      goto 99

 1099 write(6,*)'ALL DAYS DONE'

      enddo ! PM 2.5 and PM 10
c
      end

