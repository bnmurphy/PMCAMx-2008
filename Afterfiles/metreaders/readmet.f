      program readmet 
c     This program processes PMCAMx results, courtesy of JPD

      parameter(ifin = 10)
      parameter(ifout = 30)

      parameter(ncol = 97)
      parameter(nrow = 90)
      parameter(nlay = 14)
      parameter(nhr  = 24)
      parameter(nreg = 6)
      parameter(simdays=28)

      character*198 infolder
      character*4   year
      character*99  fname1
      character*99  fname2
      character*198 fname
      character*198 fnamex
      character*99  watercoverin
      character*3   cdate(simdays)
      character*2   chr(nhr)
      character*15  cldhdr

      real*4 btime,totwind,surfwind,surf,tot
      real*4 totvap(ncol,nrow),surfvap(ncol,nrow)
      real*4 temps(ncol,nrow)
      real*4 temp(ncol,nrow,nlay)
      real*4 temparray(ncol,nrow,nlay,simdays,nhr)
      real*4 dailymin(ncol,nrow,simdays)
      real*4 dailymax(ncol,nrow,simdays)
      real*4 minsum(ncol,nrow),maxsum(ncol,nrow),tsum(ncol,nrow,nlay)
      real*4 avgtmp(ncol,nrow,nlay),avgmin(ncol,nrow),avgmax(ncol,nrow)
      real*4 avgsurf(ncol,nrow)
      real*4 cldwtr(ncol,nrow,nlay)
      real*4 cldod(ncol,nrow,nlay)
      real*4 pcpwtr(ncol,nrow,nlay)
      real*4 locarray(ncol,nrow,simdays,nhr)
      real*4 dailypcp(ncol,nrow,simdays)
      real*4 hourlypcprate(ncol,nrow,simdays,nhr)
      real*4 monthlypcp(ncol,nrow)
      real*4 pcptot,mintot,maxtot,avgtot,surftot
      real   mixhgt(ncol,nrow)
      real   reglayertot(nreg)
      real*4 windspeed(ncol,nrow,nlay,nhr,simdays)
      real*4 unext(ncol,nrow,nlay),vnext(ncol,nrow,nlay)
      real*4 absu(ncol,nrow,nlay),absv(ncol,nrow,nlay) 
      real*4 humidity(ncol,nrow,nlay,simdays,nhr)
      real*4 vapor(ncol,nrow,nlay)
      real*4 regsurfsum(nreg),avgtempsum(nreg),regpcpsum(nreg)
      real*4 regtotwind(nreg),regsurfwind(nreg)
      real*4 vapsurfsum(nreg),avgvapsum(nreg)
      real   zero
 
      integer*4 ihr
      integer*4 fdate
      integer idate,ireg
      integer i,j,k,nxcl,nycl,nzcl
      integer waterfile,landcells
      integer landarray(ncol,nrow)
      integer regsize(nreg)


c      data cdate /'004','005','006','007','008','009','010','011',
c     &            '012','013','014','015','016','017','018','019',
c     &            '020','021','022','023','024','025','026','027',
c     &            '028','029'/
 
      data cdate /'185','186','187','188','189','190','191','192','193','194','195',
     &            '196','197','198','199','200','201','202','203','204','205','206',
     &            '207','208','209','210','211','212'/

      data chr /'01','02','03','04','05','06','07','08','09','10','11','12',
     &          '13','14','15','16','17','18','19','20','21','22','23','24'/

      infolder = '/home/mday/PMCAMx/met_hd/CAMx_Input/'
      year = '1991'

c Read in water cover file for averages
c      watercoverin='/usr/people/jpdawson/movie/4rpos.ladco.watercover'
c
c      waterfile = 63
c      landcells = 0
c      open(unit=waterfile,file=watercoverin,form='formatted',
c     &                            status='old')
c      do j=1,nrow
c       do i=1,ncol
c        read(waterfile,*)landarray(i,j)
c        if (landarray(i,j).eq.1) then
c         landcells = landcells + 1
c        endif
c       enddo
c      enddo
c      close(waterfile)

c      regsize(1)=1530 ! Old region sizes including more remote areas JPD
c      regsize(2)=1887
c      regsize(3)=1135
c      regsize(4)=1404
c      regsize(5)=1146
c      regsize(6)=1628

c      regsize(1)=1230 ! New regions sizes JPD
c      regsize(2)=1517
c      regsize(3)=835
c      regsize(4)=1102
c      regsize(5)=1101
c      regsize(6)=2945

c Begin temperature section

c Initialize variables
      do i=1,ncol
       do j=1,nrow
        do k=1,nlay
         tsum(i,j,k)=0.0
        enddo
        do idate=1,simdays
         dailymin(i,j,idate)=500.0
         dailymax(i,j,idate)=0.0
         do ihr=1,nhr
          do k=1,nlay
           temparray(i,j,k,idate,ihr)=0.0
          enddo
         enddo
        enddo
       enddo
      enddo

      write(6,*)'July ',year

       do idate=1,simdays
c
c open input file
c
      fname = infolder(1:INDEX(infolder,' ')-1)//year// '/' //'camx.tp.'//year//cdate(idate)//'.updt'
      isblk = INDEX(fname,' ')
      open(unit=ifin,file=fname(1:isblk-1),form='UNFORMATTED',status='old')

c
c iterate hours
c
      do ihr=1,nhr
c
c read surface temperatures
          read(ifin,END=990)btime,fdate,((temps(i,j),i=1,ncol),j=1,nrow)
c
c read hourly portion
c
       do k=1,nlay
        read(ifin,END=990)btime,fdate,((temp(i,j,k),i=1,ncol),j=1,nrow)
       enddo

       do i=1,ncol
        do j=1,nrow
         do k=1,nlay
          temparray(i,j,k,idate,ihr)=temp(i,j,k)
         enddo
        enddo
       enddo
c
       do i=1,ncol
        do j=1,nrow
         do k=1,nlay
          tsum(i,j,k)=tsum(i,j,k)+temparray(i,j,k,idate,ihr)
         enddo ! k
         if(temparray(i,j,1,idate,ihr).gt.dailymax(i,j,idate)) then
          dailymax(i,j,idate)=temparray(i,j,1,idate,ihr)
         endif
         if(temparray(i,j,1,idate,ihr).lt.dailymin(i,j,idate)) then
          dailymin(i,j,idate)=temparray(i,j,1,idate,ihr)
         endif
        enddo
       enddo 

      enddo ! ihr

      close(ifin)

      enddo ! idate

c -------------------------------------------------------------------
 990  do i=1,ncol
       do j=1,nrow
        minsum(i,j)=0.0
        maxsum(i,j)=0.0
        avgsurf(i,j)=0.0
        do k=1,nlay
         avgtmp(i,j,k)=0.0
        enddo
        avgmin(i,j)=0.0
        avgmax(i,j)=0.0
       enddo
      enddo

      do i=1,ncol
       do j=1,nrow
        do idate=1,simdays
         minsum(i,j)=minsum(i,j)+dailymin(i,j,idate)
         maxsum(i,j)=maxsum(i,j)+dailymax(i,j,idate)
        enddo
        avgmin(i,j)=minsum(i,j)/(simdays*1.0)
        avgmax(i,j)=maxsum(i,j)/(simdays*1.0)
        avgsurf(i,j)=tsum(i,j,1)/(simdays*nhr*1.0)
        do k=1,nlay
         avgtmp(i,j,k)=tsum(i,j,k)/(simdays*nhr*1.0)
        enddo
        avgmin(i,j)=avgmin(i,j)-273.15
        avgmax(i,j)=avgmax(i,j)-273.15
        avgsurf(i,j)=avgsurf(i,j)-273.15
        do k=1,nlay
         avgtmp(i,j,k)=avgtmp(i,j,k)-273.15
        enddo
       enddo
      enddo
 
      do ireg = 1,nreg
       regsurfsum(ireg) = 0.0
       avgtempsum(ireg) = 0.0
      enddo

      do i=1,ncol
       do j=1,nrow
        ireg=6

c        if(i.le.30.and.j.ge.40) ireg=1 ! Old regional definitions JPD
c        if(i.ge.31.and.i.le.67.and.j.ge.40) ireg=2
c        if(i.ge.68.and.i.le.79.and.j.ge.40.and.j.le.53) ireg=3
c        if(i.ge.68.and.i.le.86.and.j.ge.54.and.j.le.66) ireg=3
c        if(i.ge.68.and.j.ge.67) ireg=3
c        if(i.le.36.and.j.le.39) ireg=4
c        if(i.ge.37.and.i.le.73.and.j.ge.13.and.j.le.39) ireg=5
c        if(i.ge.63.and.i.le.73.and.j.ge.4.and.j.le.12) ireg=5
c        if(i.ge.74.and.i.le.79.and.j.ge.32.and.j.le.39) ireg=5

        if(i.le.30.and.j.ge.40.and.j.le.80) ireg=1 ! New regions JPD
        if(i.ge.31.and.i.le.67.and.j.ge.40.and.j.le.80) ireg=2
        if(i.ge.68.and.i.le.79.and.j.ge.40.and.j.le.53) ireg=3
        if(i.ge.68.and.i.le.86.and.j.ge.54.and.j.le.66) ireg=3
        if(i.ge.68.and.j.ge.67.and.j.le.80) ireg=3
        if(i.le.15.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.16.and.i.le.28.and.j.ge.3.and.j.le.39) ireg=4
        if(i.ge.29.and.i.le.36.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.37.and.i.le.47.and.j.ge.13.and.j.le.39) ireg=5
        if(i.ge.48.and.i.le.62.and.j.ge.16.and.j.le.39) ireg=5
        if(i.ge.63.and.i.le.73.and.j.ge.4.and.j.le.39) ireg=5
        if(i.ge.74.and.i.le.79.and.j.ge.32.and.j.le.39) ireg=5


        regsurfsum(ireg) = regsurfsum(ireg) + avgsurf(i,j)
        do k=1,nlay
         avgtempsum(ireg) = avgtempsum(ireg) + avgtmp(i,j,k)
        enddo
       enddo
      enddo

      do ireg = 1,nreg
       write(6,*)'Region',ireg
       write(6,*)'Avg surface T =', (1.0*regsurfsum(ireg))/
     &             (1.0*regsize(ireg)),'deg C'
       write(6,*)'Avg overall T =', (1.0*avgtempsum(ireg))/
     &             (1.0*regsize(ireg)*nlay),'deg C'
      enddo

      mintot=0.0
      maxtot=0.0
      surftot=0.0
      avgtot=0.0
 
      do i=1,ncol
       do j=1,nrow
        if (landarray(i,j).eq.1) then
         mintot = mintot + avgmin(i,j)
         maxtot = maxtot + avgmax(i,j)
         surftot = surftot + avgsurf(i,j)
         do k=1,nlay
          avgtot = avgtot + avgtmp(i,j,k)
         enddo
        endif
       enddo
      enddo

      write(6,*)' '
      write(6,*)'Average min T  =',mintot/(landcells*1.0),'deg C'       
      write(6,*)'Average surf T =',surftot/(landcells*1.0),'deg C'
      write(6,*)'Average max T  =',maxtot/(landcells*1.0),'deg C'
      write(6,*)'Average T      =',avgtot/(landcells*nlay*1.0),'deg C'
 
c End of temperature section

c Beginning of cloud/rain section

       do idate=1,simdays
        do ihr=1,nhr
         do i=1,ncol
          do j=1,nrow
           locarray(i,j,idate,ihr)=0.0
           dailypcp(i,j,idate)=0.0
          enddo
         enddo
        enddo
       enddo

       do idate=1,simdays

c open input file
      fname=infolder(1:INDEX(infolder,' ')-1)//year// '/' //'camx.cr.'//year//cdate(idate)//'.updt'
      isblk=INDEX(fname,' ')
      open(unit=ifin,file=fname(1:isblk-1),form='UNFORMATTED',
     &     status='old')

c read header
      read(ifin,END=991)cldhdr,nxcl,nycl,nzcl

c iterate hours
c 
      do ihr=1,nhr
c 
c read hourly portion
c
       read(ifin,END=991)btime,fdate
       do k=1,nlay
        read(ifin,END=991)((cldwtr(i,j,k),i=1,ncol),j=1,nrow)
        read(ifin,END=991)((pcpwtr(i,j,k),i=1,ncol),j=1,nrow)
        read(ifin,END=991)((cldod(i,j,k),i=1,ncol),j=1,nrow)
       enddo
c
       do i=1,ncol
        do j=1,nrow
          locarray(i,j,idate,ihr)=pcpwtr(i,j,1)
          locarray(i,j,idate,ihr)=pcpwtr(i,j,1)
          locarray(i,j,idate,ihr)=pcpwtr(i,j,1)
          locarray(i,j,idate,ihr)=pcpwtr(i,j,1)
        enddo
       enddo
c
       do i=1,ncol
        do j=1,nrow
         hourlypcprate(i,j,idate,ihr) =
     &               18.62*(locarray(i,j,idate,ihr)**1.27) ! convert to mm/hr
        enddo
       enddo

      enddo ! ihr
 
      close(ifin)
 
      do i=1,ncol
       do j=1,nrow
        do ihr=1,nhr
         dailypcp(i,j,idate)=dailypcp(i,j,idate)+
     &                          hourlypcprate(i,j,idate,ihr) ! mm/day
        enddo ! ihr

       enddo ! nrow
      enddo ! ncol

      enddo ! idate
 
c -------------------------------------------------------------------
 991  do i=1,ncol
       do j=1,nrow
       monthlypcp(i,j)=0.0
       enddo
      enddo

      do i=1,ncol
       do j=1,nrow
        do idate=1,simdays
         monthlypcp(i,j)=monthlypcp(i,j)+dailypcp(i,j,idate)
        enddo
       enddo
      enddo

      do ireg = 1,nreg
       regpcpsum(ireg) = 0.0
      enddo

      do i=1,ncol
       do j=1,nrow
        ireg=6
        if(i.le.30.and.j.ge.40.and.j.le.80) ireg=1 ! New regions JPD
        if(i.ge.31.and.i.le.67.and.j.ge.40.and.j.le.80) ireg=2
        if(i.ge.68.and.i.le.79.and.j.ge.40.and.j.le.53) ireg=3
        if(i.ge.68.and.i.le.86.and.j.ge.54.and.j.le.66) ireg=3
        if(i.ge.68.and.j.ge.67.and.j.le.80) ireg=3
        if(i.le.15.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.16.and.i.le.28.and.j.ge.3.and.j.le.39) ireg=4
        if(i.ge.29.and.i.le.36.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.37.and.i.le.47.and.j.ge.13.and.j.le.39) ireg=5
        if(i.ge.48.and.i.le.62.and.j.ge.16.and.j.le.39) ireg=5
        if(i.ge.63.and.i.le.73.and.j.ge.4.and.j.le.39) ireg=5
        if(i.ge.74.and.i.le.79.and.j.ge.32.and.j.le.39) ireg=5

        regpcpsum(ireg) = regpcpsum(ireg) + monthlypcp(i,j)
       enddo
      enddo

      do ireg = 1,nreg
       write(6,*)'Region',ireg
       write(6,*)'Avg precip =', (1.0*regpcpsum(ireg))/
     &             (1.0*regsize(ireg)),'mm'
      enddo

      pcptot=0.0

      do i=1,ncol
       do j=1,nrow
        if (landarray(i,j).eq.1) then
         pcptot = pcptot + monthlypcp(i,j)
        endif
       enddo
      enddo

      write(6,*)
      write(6,*)'Average monthly precip = ',pcptot/(landcells*1.0),
     &          ' mm'
      
c End cloud/precip section

c Begin mixing height section

        do i=1,ncol
         do j=1,nrow
          mixhgt(i,j)=0.0
         enddo
        enddo
c
c open input file
c
c      fname='/usr/people/jpdawson/CoupledModelInputs/JulyPbl/'//
c     &        'July.1992.pbl.avg'
c      isblk=INDEX(fname,' ')
c      open(unit=ifin,file=fname(1:isblk-1),form='FORMATTED',
c     &     status='old')
c
c      do j=1,nrow
c       do i=1,ncol
c        read(ifin,*)mixhgt(i,j)
c       enddo
c      enddo
c
c      close(ifin)

c -------------------------------------------------------------------
c 992  do ireg=1,nreg
c       reglayertot(ireg)=0.0
c      enddo 
c
c      do i=1,ncol
c       do j=1,nrow
c        ireg=6
c        if(i.le.30.and.j.ge.40.and.j.le.80) ireg=1 ! New regions JPD
c        if(i.ge.31.and.i.le.67.and.j.ge.40.and.j.le.80) ireg=2
c        if(i.ge.68.and.i.le.79.and.j.ge.40.and.j.le.53) ireg=3
c        if(i.ge.68.and.i.le.86.and.j.ge.54.and.j.le.66) ireg=3
c        if(i.ge.68.and.j.ge.67.and.j.le.80) ireg=3
c        if(i.le.15.and.j.ge.13.and.j.le.39) ireg=4
c        if(i.ge.16.and.i.le.28.and.j.ge.3.and.j.le.39) ireg=4
c        if(i.ge.29.and.i.le.36.and.j.ge.13.and.j.le.39) ireg=4
c        if(i.ge.37.and.i.le.47.and.j.ge.13.and.j.le.39) ireg=5
c        if(i.ge.48.and.i.le.62.and.j.ge.16.and.j.le.39) ireg=5
c        if(i.ge.63.and.i.le.73.and.j.ge.4.and.j.le.39) ireg=5
c        if(i.ge.74.and.i.le.79.and.j.ge.32.and.j.le.39) ireg=5
c
c        if(mixhgt(i,j).gt.2500) mixhgt(i,j)=2500.0
c        if(mixhgt(i,j).lt.1) mixhgt(i,j)=1.0
c
c        reglayertot(ireg)=reglayertot(ireg)+mixhgt(i,j)
c
c       enddo
c      enddo 
c
c      do ireg = 1,nreg
c       write(6,*)'Region',ireg
c       write(6,*)'Avg mix hgt  =', reglayertot(ireg)/
c     &             (1.0*regsize(ireg)),' m'
c      enddo

c End mixing height section

c Begin wind speed section

      do i=1,ncol
       do j=1,nrow
        do ihr=1,nhr
         do idate=1,simdays
          do k=1,nlay
           windspeed(i,j,k,ihr,idate)=0.0
          enddo
         enddo
        enddo
       enddo
      enddo

      do idate=1,simdays

      fname=infolder(1:INDEX(infolder,' ')-1)// 'camx.uv.'//year//cdate(idate)//'.updt'
      isblk=INDEX(fname,' ')
      open(unit=ifin,file=fname(1:isblk-1),form='UNFORMATTED',
     &     status='old')

      do ihr=1,nhr
      
       read(ifin,end=993)btime,fdate 
       do k = 1,nlay
         read(ifin,end=993) ((unext(i,j,k),i=1,ncol),j=1,nrow)
         read(ifin,end=993) ((vnext(i,j,k),i=1,ncol),j=1,nrow)

          do i=1,ncol
           do j=1,nrow
            absu(i,j,k)=abs(unext(i,j,k))
            absv(i,j,k)=abs(vnext(i,j,k))
            windspeed(i,j,k,ihr,idate)=sqrt((absu(i,j,k)*absu(i,j,k))+
     &                                    (absv(i,j,k)*absv(i,j,k)))      
           enddo
          enddo
       
       enddo ! nlay

      read(ifin,end=993)zero

      enddo ! ihr  

      close(ifin)

      enddo ! idate

c ------------------------------------------------------------
 993  totwind = 0.0
      surfwind = 0.0
      do ireg = 1,nreg
       regtotwind(ireg) = 0.0
       regsurfwind(ireg) = 0.0
      enddo
      
      do i=1,ncol
       do j=1,nrow

        ireg=6
        if(i.le.30.and.j.ge.40.and.j.le.80) ireg=1 ! New regions JPD
        if(i.ge.31.and.i.le.67.and.j.ge.40.and.j.le.80) ireg=2
        if(i.ge.68.and.i.le.79.and.j.ge.40.and.j.le.53) ireg=3
        if(i.ge.68.and.i.le.86.and.j.ge.54.and.j.le.66) ireg=3
        if(i.ge.68.and.j.ge.67.and.j.le.80) ireg=3
        if(i.le.15.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.16.and.i.le.28.and.j.ge.3.and.j.le.39) ireg=4
        if(i.ge.29.and.i.le.36.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.37.and.i.le.47.and.j.ge.13.and.j.le.39) ireg=5
        if(i.ge.48.and.i.le.62.and.j.ge.16.and.j.le.39) ireg=5
        if(i.ge.63.and.i.le.73.and.j.ge.4.and.j.le.39) ireg=5
        if(i.ge.74.and.i.le.79.and.j.ge.32.and.j.le.39) ireg=5

        do k=1,nlay
         do ihr=1,nhr
          do idate=1,simdays
           regtotwind(ireg)=regtotwind(ireg)+windspeed(i,j,k,ihr,idate)
           if (k.eq.1) then
            regsurfwind(ireg)=regsurfwind(ireg)+
     &               windspeed(i,j,k,ihr,idate)
           endif
           if (landarray(i,j).eq.1) then
            totwind = totwind + windspeed(i,j,k,ihr,idate)
            if (k.eq.1) then
             surfwind = surfwind + windspeed(i,j,k,ihr,idate)
            endif
           endif
          enddo
         enddo
        enddo
       enddo
      enddo

      do ireg = 1,nreg
       write(6,*)'Region',ireg
       write(6,*)'Reg surface wind spd =',
     &         (regsurfwind(ireg)/(regsize(ireg)*nhr*simdays*1.0)),'m/s'
       write(6,*)'Reg avg wind spd =',
     &     (regtotwind(ireg)/(regsize(ireg)*nhr*simdays*nlay*1.0)),'m/s'
      enddo

      write(6,*)' '
      write(6,*)'Avg surface wind speed =',
     &              (surfwind/(landcells*nhr*simdays*1.0)),'m/s'
      write(6,*)'Avg overall wind speed =',
     &              (totwind/(landcells*nhr*simdays*nlay*1.0)),'m/s'

c End wind speed section

c Begin absolute humidity section

      do i=1,ncol
       do j=1,nrow
        do ihr=1,nhr
         do idate=1,simdays
          do k=1,nlay
           humidity(i,j,k,ihr,idate)=0.0
          enddo
         enddo
        enddo
       enddo
      enddo

      do idate=1,simdays

      fname=infolder(1:INDEX(infolder,' ')-1)// year// '/' //'camx.qa.'//year//cdate(idate)//'.updt'
      isblk=INDEX(fname,' ')
      open(unit=ifin,file=fname(1:isblk-1),form='UNFORMATTED',
     &     status='old')

      do ihr=1,nhr

       do k=1,nlay
        read(ifin,END=994)btime,fdate,((vapor(i,j,k),i=1,ncol),
     &                                  j=1,nrow) 
       enddo ! k

       do k=1,nlay
        do i=1,ncol
         do j=1,nrow
           humidity(i,j,k,ihr,idate)=vapor(i,j,k)
         enddo
        enddo
       enddo

      enddo ! ihr

      close(ifin)
      enddo ! idate

c --------------------------------------------------------
 994  do i=1,ncol
       do j=1,nrow
        totvap(i,j) = 0.0
        surfvap(i,j) = 0.0
       enddo
      enddo

c      watercoverin='/usr/people/jpdawson/movie/4rpos.ladco.watercover'

c      waterfile = 63
      landcells = 90*97*2/3
c      open(unit=waterfile,file=watercoverin,form='formatted',
c     &                            status='old')
c      do j=1,nrow
c       do i=1,ncol
c        read(waterfile,*)landarray(i,j)
c        if (landarray(i,j).eq.1) then
c         landcells = landcells + 1
c        endif
c       enddo
c      enddo
c      close(waterfile)

      do idate=1,simdays
       do ihr=1,nhr
        do i=1,ncol
         do j=1,nrow
          surfvap(i,j) = surfvap(i,j) + humidity(i,j,1,ihr,idate)
          do k=1,nlay
           totvap(i,j) = totvap(i,j) + humidity(i,j,k,ihr,idate)
          enddo
         enddo
        enddo
       enddo
      enddo

      do ireg = 1,nreg
       vapsurfsum(ireg) = 0.0
       avgvapsum(ireg) = 0.0
      enddo

      do i=1,ncol
       do j=1,nrow
        ireg=6
        if(i.le.30.and.j.ge.40.and.j.le.80) ireg=1 ! New regions JPD
        if(i.ge.31.and.i.le.67.and.j.ge.40.and.j.le.80) ireg=2
        if(i.ge.68.and.i.le.79.and.j.ge.40.and.j.le.53) ireg=3
        if(i.ge.68.and.i.le.86.and.j.ge.54.and.j.le.66) ireg=3
        if(i.ge.68.and.j.ge.67.and.j.le.80) ireg=3
        if(i.le.15.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.16.and.i.le.28.and.j.ge.3.and.j.le.39) ireg=4
        if(i.ge.29.and.i.le.36.and.j.ge.13.and.j.le.39) ireg=4
        if(i.ge.37.and.i.le.47.and.j.ge.13.and.j.le.39) ireg=5
        if(i.ge.48.and.i.le.62.and.j.ge.16.and.j.le.39) ireg=5
        if(i.ge.63.and.i.le.73.and.j.ge.4.and.j.le.39) ireg=5
        if(i.ge.74.and.i.le.79.and.j.ge.32.and.j.le.39) ireg=5

        vapsurfsum(ireg) = vapsurfsum(ireg) + surfvap(i,j)
        avgvapsum(ireg) = avgvapsum(ireg) + totvap(i,j)
       enddo
      enddo

      do ireg = 1,nreg
       write(6,*)'Region',ireg
       write(6,*)'Avg surface humidity =', (1.0*vapsurfsum(ireg))/
     &             (1.0*nhr*simdays*regsize(ireg)),'ppm'
       write(6,*)'Avg overall humidity =', (1.0*avgvapsum(ireg))/
     &             (1.0*nhr*simdays*regsize(ireg)*nlay),'ppm'
      enddo


      do i=1,ncol
       do j=1,nrow
        surf = surf + surfvap(i,j)
        tot  = tot  + totvap(i,j)
       enddo
      enddo
 
      write(6,*)' '
      write(6,*)'Avg surface humidity =',
     &           ((1.0*surf)/(landcells*nhr*simdays*1.0)),'ppm'
      
      write(6,*)'Avg overall humidity =',
     &            ((1.0*tot)/(landcells*nhr*simdays*nlay*1.0)),'ppm'
      write(6,*)' '

c End water vapor section

      end

