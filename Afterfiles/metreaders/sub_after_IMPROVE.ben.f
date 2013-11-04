      subroutine after(dirname, filename, pname, totcon, ncomp)
c     This program processes PMCAMx results
c     The parameters are hardcoded for the specific output file
      
      parameter(nfile = 1) ! number of output file to read
                           ! if 1 read HYBR only
c     parameter(isgas = 0) ! if 1 - GAS is the species to be processed
      parameter(isavg = 1) ! 1 if 24-hr avg; 2 if 1-hr avg peak; 3 if hourly

      parameter(ifinh = 10)
      parameter(ifine = 20)
      parameter(ifout = 30)

      parameter(ncol = 97)
      parameter(nrow = 90)
      parameter(nlay = 1)
c      parameter(nlay = 14) ! PITT only
      parameter(nspc = 459) ! Number of avg species
c      parameter(nspc = 204) ! Number of avg species
      parameter(nhr  = 24)
      parameter(mspc = 1) ! Number of species to be processed
      parameter(nms  = 1) ! Number of Monitoring Sites



      character*99 infile 
      character*4 pname(mspc)
      character*3 cdate(1)
      character*99 dirname
      character*99 filename
      integer ncomp ! Number of comparisons
      integer isgas ! if 1 GAS is the species to be processed

      character*2 hnam(24)
      character*99 fname1
      character*99 fname2
      character*198 fname
      character*99 fname3
      character*99 fname4
      character*198 fnamex
      character*4  name(10),note(60),mspec(10,nspc)
      character*10 spcname(nspc),tspcname(nspc)


      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 conh(nspc,ncol,nrow,nhr),cone(nspc,ncol,nrow,nhr)
      real*4 sconh(mspc,ncol,nrow),scone(mspc,ncol,nrow)
      real*4 mbias(mspc),mnbias(mspc),merr(mspc),mnerr(mspc)
      real*4 maxh(mspc),maxe(mspc)
      real*4 tmpcon(mspc,nlay,nhr)
      real*4 pitt(nspc,nlay,nhr) ! PITT
      real*4 totcon(ncomp,ncol,nrow)
      real Temperature(17,24,8730)

      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate
      integer*4 idxspc,ihr,msec
      integer*4 nsum(mspc)
      integer*4 maxih,maxjh,maxie,maxje
      integer*4 msx,msy
      integer*2 bin

c      data cdate /'193','194','195','196','197','198','199'/
c'200',     &      '201','202','203','204','205','206','207','208','209'/
c      data pname /'O3  ','NO  ','NO2 ','SO2 ','NH3 ','HNO3','CO  '/ ! Species to be processed
c      data pname /'NO  ','NO2 ','CG1 ','CG2 ','CG3 ','CG4  '/ ! Species to be processed
c      data pname /'MEOH','ETOH','CRES','XYL ','FORM','ALD2 '/ ! Species to be processed
c       data pname /'HNO3'/
c      data pname /'PALL','PTOC','PEC ','PNH4','PNO3','PSO4','PH2O',
c     &     'NH3','HNO3','O3','CO'/
c      data pname /'PALL','PTOC','SOA1','SOA2','SOA3','SOA4'/
      data msec  / 6 /
c PITT      data msx / 65 / ! x-coordinate of monitoring sites
c PITT      data msy / 51 / ! y-coordinate of monitoring sites
c      data msx / 65 / ! x-coordinate of monitoring sites
c      data msy / 51 / ! y-coordinate of monitoring sites
       data hnam /'01','02','03','04','05','06','07','08','09','10','11',
     & '12','13','14','15','16','17','18','19','20','21','22','23','24'/
c
c iterate days
c
      idate=0
 99   idate=idate+1
      if(idate.gt.1) goto 1099

       
c       do ih=1,24
c        infile='/home/mkshriva/work/Temperature/temp.'
c     &          //cdate(idate)//'.'//hnam(ih)
c        open(unit=525,FILE=infile,STATUS='old',FORM='FORMATTED')
c        do i=1,8730
c         read(525,*)Temperature(idate,ih,i)
c        enddo
c       enddo
      


c
c open output file if necessary
c
c      if(ifout.ne.6) then
c       do i=1,mspc
c         isblk=INDEX(pname(i),' ')
c         if(isblk.eq.0) isblk = 5
c         if(isavg.eq.1) then
c          fnamex=pname(i)(1:isblk-1)//'.avg.'//cdate(idate)
c         elseif(isavg.eq.2) then
c          fnamex=pname(i)(1:isblk-1)//'.peak.'//cdate(idate)
c         elseif(isavg.eq.3) then
c          fnamex=pname(i)(1:isblk-1)//'.orig.hourly.'//cdate(idate)
c         endif
c         isblk=INDEX(fnamex,' ')
c         ifoutx = ifout + i
c         open(unit=ifoutx,file=fnamex(1:isblk-1),status='new')
c       enddo
c      endif
c
c open input file - HYBR
c
c      fname1='/home/bkoo/pmcamxv3.11/outputs/2001_'//cdate//'/'
c      fname2='4rpos.basecase.'//cdate//'.2001.avrg'
c      fname1='/usr/people/tmg/kfahey/NEUS-settling/NEUS-settling/outputs/2001/'
c      fname1='/usr/people/tmg/albert/2001/'
      fname1=dirname
      fname2=filename
      isblk=INDEX(fname1,' ')
      fname=fname1(1:isblk-1)//fname2
      isblk=INDEX(fname,' ')
c      write(6,*)'HYBR: ',fname(1:isblk-1)
      open(unit=ifinh,file=fname(1:isblk-1),form='UNFORMATTED',
     &     status='old')
c
c read header portion - HYBR
c
      read(ifinh,END=990)name,note,ione,nspec,ibdate,btime,iedate,etime
      write(6,*)(name(i)(1:1),i=1,10)
      write(6,*)(note(i)(1:1),i=1,60)
      write(6,*)ione,nspec,ibdate,btime,iedate,etime
      if(nspec.ne.nspc)stop'nspec is not equal to nspc'
      read(ifinh,END=990)rdum,rdum,iutm,xorg,yorg,delx,dely,nx,ny,nz
     &                 ,idum,idum,rdum,rdum,rdum
      write(6,*)iutm,xorg,yorg,delx,dely
      write(6,*)nx,ny,nz
      read(ifinh,END=990)izero,izero,nx,ny
      write(6,*)izero,nx,ny
      if(nx.ne.ncol)stop'nx is not equal to ncol'
      if(ny.ne.nrow)stop'ny is not equal to nrow'
      if(nz.ne.nlay)stop'nz is not equal to nlay'

      read(ifinh,END=990)((mspec(i,j),i=1,10),j=1,nspec)

c    START of hourly loop
 100  continue
      read(ifinh,END=990)ibdate,btime,iedate,etime
      ihr = NINT(btime)
      write(6,*)'ihr = ',ihr
      if(ihr+1.gt.nhr)stop'incorrect hr index'
      do isp=1,nspec
       do ilay=1,nz
        read(ifinh,END=990)ione,(mspec(n,isp),n=1,10),
     &                    ((conh(isp,i,j,ihr+1),i=1,nx),j=1,ny)
        write(spcname(isp),*)(mspec(n,isp)(1:1),n=1,10)
        pitt(isp,ilay,ihr+1) = conh(isp,msx,msy,ihr+1) ! PITT
       enddo
      enddo

      if(ihr+1.eq.nhr) goto 150 ! done reading hourly data

      goto 100
c
c calculate desired quantities
c
 150  continue

      do m=1,ncomp
         do i=1,nx
         do j=1,ny
            totcon(m,i,j)=0.0
c           do khr=1,nhr
               sconh(m,i,j)=0.0
c           enddo
         enddo
         enddo
         do khr=1,nhr
            do ilay=1,nlay
               tmpcon(m,ilay,khr)=0.0
            enddo
         enddo

      do isp=1,nspec
       if(isgas.eq.1) then ! GASES
        do i=1,nx
        do j=1,ny
         do khr=1,nhr
           if(spcname(isp)(1:4).eq.pname(m)) then
            if(isavg.eq.1) then
             totcon(m,i,j)=totcon(m,i,j)+conh(isp,i,j,khr)/24.0 ! sum
            elseif(isavg.eq.2) then
             sconh(m,i,j)=max(sconh(m,i,j),conh(isp,i,j,khr)) ! 1-hr peak
            elseif(isavg.eq.3) then
             if(i.eq.msx.and.j.eq.msy) then
              ifoutx = ifout + m
c              write(ifoutx,'(A3,I11,E14.5)')cdate(idate),khr,conh(isp,i,j,khr)
c              write(ifoutx,777)cdate(idate),khr,(pitt(isp,ilay,khr)
c     &             ,ilay=1,nlay) ! PITT only
             endif
            endif
           endif
         enddo   !khr
        enddo    !j
        enddo    !i
       else                ! AEROSOLS
       isund=INDEX(spcname(isp),'_')
       if (isund.gt.0) 
     &	  read (spcname(isp)(isund+1:isund+2), *) bin
       if(isund.gt.0 .and. bin.eq.1) then
          do i=1,nx
          do j=1,ny
            do khr=1,nhr
               if(pname(m).eq.'PTOC'.and.
     &		  spcname(isp)(1:isund-1).eq.'POC') then
          	   
		   if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                   do mm=1,msec ! section 1 through section msec
                      totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                   enddo
               endif
               if(pname(m).eq.'PNH4' .and.
     &	            spcname(isp)(1:isund-1).eq.'PNH4') then
c     &         	     spcname(isp)(1:isund-1).eq.'AA52'.or.
c     &         	     spcname(isp)(1:isund-2).eq.'AO1'.or.
c     &         	     spcname(isp)(1:isund-2).eq.'AO2')
c     &             ) then
c          	  write(6,*), pname(m), spcname(isp)
		   if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
           	  do mm=1,msec ! section 1 through section msec
             	     totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
           	  enddo
               endif
               if(pname(m).eq.'PNO3' .and.
     &	            spcname(isp)(1:isund-1).eq.'PNO3' ) then
c     &               spcname(isp)(1:isund-2).eq.'AR2')
c     &             ) then
c          	  write(6,*), pname(m), spcname(isp)
		   if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                  do mm=1,msec ! section 1 through section msec
                    totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                  enddo
               endif
               if(pname(m).eq.'PSO4'.and.
     &	           spcname(isp)(1:isund-1).eq.'PSO4') then
c     &              spcname(isp)(1:isund-2).eq.'ATR'.or.
c     &              spcname(isp)(1:isund-2).eq.'ANP')
c     &            ) then
c          	  write(6,*), pname(m), spcname(isp)
		   if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                  do mm=1,msec ! section 1 through section msec
                   totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                  enddo
               endif
               if(pname(m).eq.'PEC'.and.
     &	           spcname(isp)(1:isund-1).eq.'PEC') then
c     &              spcname(isp)(1:isund-2).eq.'ATR'.or.
c     &              spcname(isp)(1:isund-2).eq.'ANP')
c     &            ) then
c          	  write(6,*), pname(m), spcname(isp)
		   if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                  do mm=1,msec ! section 1 through section msec
                   totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                  enddo
               endif
               if (pname(m).eq.'PM25'.and.
     &	          (spcname(isp)(1:isund-2).eq.'APO'.or.
     &             spcname(isp)(1:isund-2).eq.'AOO'.or.
     &		   spcname(isp)(1:isund-1).eq.'POC'.or.
     &		   spcname(isp)(1:isund-2).eq.'ABS'.or.
     &		   spcname(isp)(1:isund-2).eq.'AAS'.or.
     &		   spcname(isp)(1:isund-1).eq.'PEC'.or.
     &		   spcname(isp)(1:isund-1).eq.'CRST'.or.
     &		   spcname(isp)(1:isund-1).eq.'PCL'.or.
     &		   spcname(isp)(1:isund-1).eq.'NA'.or.
     &		   spcname(isp)(1:isund-1).eq.'PH2O'.or.
     &		   spcname(isp)(1:isund-1).eq.'PNH4'.or.
     &             spcname(isp)(1:isund-2).eq.'PNO3'.or.
     &             spcname(isp)(1:isund-2).eq.'PSO4')
     &            ) then
c          	  write(6,*), pname(m), spcname(isp)
		   if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                  do mm=1,msec ! section 1 through section msec
                   totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)*1.5/24.0
                  enddo
               endif
            enddo ! khr
          enddo ! j
          enddo ! i
       endif ! isund > 0 ?
       endif ! isgas = 1 ?
      enddo ! isp
      enddo ! m
 990  continue
c     write(6,*)'END OF FILE: conc. average input file'

      close(ifinh)
      if(nfile.ne.1) close(ifine)
      if(ifout.ne.6) then
       do i=1,mspc
         ifoutx = ifout + i
c         close(ifoutx)
       enddo
      endif

      goto 99

 1099 continue
c     write(6,*)'ALL DAYS DONE'

      end


