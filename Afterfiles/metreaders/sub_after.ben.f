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
      parameter(nspc = 459) ! Number of avg species
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

      data msec  / 6 /

c
c iterate days
c
      idate=0
 99   idate=idate+1
      if(idate.gt.1) goto 1099

      
c
c open input file - HYBR
c
      fname1=dirname
      fname2=filename
      isblk=INDEX(fname1,' ')
      fname=fname1(1:isblk-1)//fname2
      isblk=INDEX(fname,' ')
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

      if(ihr+1.eq.nhr) goto 150   !done reading hourly data

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
             endif
            endif
           endif
         enddo   !khr
        enddo    !j
        enddo    !i
       else                ! AEROSOLS
       isund=INDEX(spcname(isp),'_')
       if (isund.gt.0) 
     &	  read (spcname(isp)(isund+1:isund+3), *) bin
       if(isund.gt.0 .and. bin.eq.1) then
          do i=1,nx
          do j=1,ny
            do khr=1,nhr
               if(pname(m).eq.'PTOC'.and.
     &		  (spcname(isp)(1:isund-1).eq.'POC'.or.
     &		   spcname(isp)(1:isund-2).eq.'APO'.or.
     &             spcname(isp)(1:isund-3).eq.'APO'.or.
     &		   spcname(isp)(1:isund-2).eq.'AOO'.or.
     &		   spcname(isp)(1:isund-3).eq.'AOO'.or.
     &		   spcname(isp)(1:isund-2).eq.'ABS'.or.
     &		   spcname(isp)(1:isund-2).eq.'AAS')
     &		  ) then
		      if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                      do mm=1,msec ! section 1 through section msec
                         totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                      enddo
               endif
               if(pname(m).eq.'PNH4' .and.
     &	         spcname(isp)(1:isund-1).eq.'PNH4') then
		      if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
           	      do mm=1,msec ! section 1 through section msec
             	         totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
           	      enddo
               endif
               if(pname(m).eq.'PNO3' .and.
     &	         spcname(isp)(1:isund-1).eq.'PNO3' ) then
		      if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                      do mm=1,msec ! section 1 through section msec
                         totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                      enddo
               endif
               if(pname(m).eq.'PSO4'.and.
     &	         spcname(isp)(1:isund-1).eq.'PSO4') then
		      if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                      do mm=1,msec ! section 1 through section msec
                         totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                      enddo
               endif
               if(pname(m).eq.'PEC'.and.
     &	         spcname(isp)(1:isund-1).eq.'PEC') then
		      if (i.eq.1.and.j.eq.1.and.khr.eq.1) write(6,*), pname(m), spcname(isp)
                      do mm=1,msec ! section 1 through section msec
                         totcon(m,i,j)=totcon(m,i,j)+conh(isp-1+mm,i,j,khr)/24.0
                      enddo
               endif
               if (pname(m).eq.'PM25'.and.
     &	          (spcname(isp)(1:isund-2).eq.'APO' .or.
     &		   spcname(isp)(1:isund-3).eq.'APO' .or.
     &             spcname(isp)(1:isund-2).eq.'AOO' .or.
     &		   spcname(isp)(1:isund-3).eq.'AOO' .or.
     &		   spcname(isp)(1:isund-1).eq.'POC' .or.
     &		   spcname(isp)(1:isund-2).eq.'ABS' .or.
     &		   spcname(isp)(1:isund-2).eq.'AAS' .or.
     &		   spcname(isp)(1:isund-1).eq.'PEC' .or.
     &		   spcname(isp)(1:isund-1).eq.'CRST'.or.
     &		   spcname(isp)(1:isund-1).eq.'PCL' .or.
     &		   spcname(isp)(1:isund-1).eq.'NA'  .or.
     &		   spcname(isp)(1:isund-1).eq.'PNH4'.or.
     &             spcname(isp)(1:isund-1).eq.'PNO3'.or.
     &             spcname(isp)(1:isund-1).eq.'PSO4')
     &            ) then
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


