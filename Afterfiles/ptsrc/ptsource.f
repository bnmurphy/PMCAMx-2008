      PROGRAM ptsource

c to determine where the point sources are 

      parameter(ifinp = 10)
      parameter(ifina = 20)
      parameter(ifoutp = 30)
      parameter(ifouta = 40)
      parameter(simdays = 1) !num of simulation days
      parameter(max = 6000)
      parameter(kmax = 10)
      parameter(nispmax = 800)
      parameter(bmax=50)

      character*198 pemis, aemis, fname, p_out, a_out
      character*198 outdir, indir, scenario
      character*3   cdate(simdays)   !Julian Date
      character*4   name(10),note(60)

      integer n,t,isp,k,b,nisp
      real*4 nspec,ibdate,btime,ietime,iedate,etime,iutm
      real*4 rdum,xorg,yorg,delx,dely,nx,ny,nz
      real*4 nptsrc,idum, xloc(max),yloc(max),hstk(max),dstk(max)
      real*4 tstk(max),vstk(max),xstk(max),ystk(max),ihr,flowrat(max),effph(max)
      real*4 map(max,5),indx,DAY(365)
      real*4 ione,mspec(max,nispmax),ptemis1(max),gamma,p0,press(bmax,kmax),dstkabs
      real*4 hght1d(kmax),height(kmax),temp1d(kmax),temp(bmax,kmax)
      real*4 w2,xwind(bmax),ywind(bmax)
      real*4 wind1d(kmax),dtheta,dtdz1d(kmax)
      real*4 pispid(nispmax),pptemis(bmax,nispmax,kmax)
      real*4 gispid(nispmax),gptemis(bmax,nispmax,kmax),x(bmax),y(bmax)

      intrinsic amax1,amin1
      external plumeris

      data cdate /'193'/
c ,'194','195','196','197','198','199','200',
c     &            '201','202','203','204','205','206','207','208',
c     &            '209'/

      indir = '/home/BaseInput/July/SAPRC/prep99/semivol/'
      outdir = '/home/mday/processed_output/'
      scenario = 'sites/'    !output folder

      do idate = 1,simdays

c Open output files
      p_out = outdir(1:INDEX(outdir,' ')-1)//
     &             scenario(1:INDEX(scenario,' ')-1)//'pointemis.'//cdate(idate)
          open(unit=ifoutp,file=p_out,status='UNKNOWN')

      a_out = outdir(1:INDEX(outdir,' ')-1)//
     &             scenario(1:INDEX(scenario,' ')-1)//'areaemis.'//cdate(idate)
          open(unit=ifouta,file=a_out,status='UNKNOWN')

c Open point emissions file

      pemis = indir(1:INDEX(indir,' ')-1)//'point.camx.9902.2001'//cdate(idate)//'.bin'
      write(6,*)'INPUT FILE: ',pemis
      open(unit=ifinp, file=pemis, form='UNFORMATTED', status='old')

        read(ifinp,End=990)name,note,ione,nspec,ibdate,btime,iedate,etime
c       write(6,*)(name(i)(1:1),i=1,10),(note(i)(1:1),i=1,60)
cc      write(6,*)(note(i)(1:1),i=1,60)
cc      write(6,*)nspec,ibdate,btime,iedate,etime
        read(ifinp,end=990)rdum,rdum
c,iutm,xorg,yorg,delx,dely,nx,ny,nz
c     &                  ,idum,idum,rdum,rdum,rdum
	read(ifinp,END=990)ione,ione,nx,ny
cc	write(6,*)rdum, iutm
cc      write(6,*)iutm,xorg,yorg,delx,dely
cc      write(6,*)nx,ny,nz
c      	read(ifin, END=990)ione,nstk
c	read(ifin,END=990)(xstk(n),ystk(n),hstk(n),dstk(n),tstk(n),vstk(n),
c     &			n=1,nstk)

c      print *,'\n','Stack # ',nstk,'\n\n'

c   READ STACK INFORMATION, TEST LOCATION, AND STORE FOR LATER IF PASSED
         read(ifinp) idum, nptsrc
         read(ifinp) (xloc(n),yloc(n),hstk(n),dstk(n),tstk(n),vstk(n),
     &             n=1,nptsrc)

C      Map pt sources to grid
       do n=1,nptsrc
         xstk(n)=xloc(n)/1000.+900.
         ystk(n)=yloc(n)/1000.+1620.
         map(n,1)=1
         map(n,2)=1+INT(xstk(n)/36)
         map(n,3)=1+INT(ystk(n)/36)
	 write(ifinp,*) map
       enddo

        print *,'NUMBER OF POINT SOURCES= ',nptsrc

C   The rest of this has not been debugged.........(also missing header parts)
C   READ EMISSIONS AND CLOSE FILE
c              do ihr = 1,24
c          indx = (DAY(t)-DAY(1))*24*6 + (ihr-1)*6 + 1
c          read(ifinp)ibdate,btime,iedate,etime        !TIMES
c
c          read(ifinp) idum,npts                        !NUMBER OF POINT SOURCES
c          read(ifinp) (idum,idum,idum,flowrat(n),effph(n),n=1,npts) !PT SOURCE SPECS

C DEBUG
c        print *,'ibdate= ',ibdate,'  btime= ',btime,'  iedate= ',iedate,
c     &		'  etime= ',etime
c
c          do isp = 1,nspec
c                       read(ifinp)ione,(mspec(n,isp),n=1,10),
c     &            (ptemis1(n),n=1,npts)
c
c                do n = 1,npts
c                  do b = indx,indx+5
c                  if (map(n,2).eq.x(b) .and. map(n,3).eq.y(b) ) then

C                CALCULATE PLUME RISE

c                if (effph(n) .lt. 0.) then
c                      zstk = abs(effph(n))
c                      goto 14
c                  endif
c                  do k = 1,10
c                      hght1d(k) = height(k)
c                      temp1d(k) = temp(b,k)
c                      w2 = xwind(b)*xwind(b) + ywind(b)*ywind(b)
c                      wind1d(k) = amax1(sqrt(w2),0.1)
c                      if (k.lt.10) then
c                          dz = height(k+1)/2.
c                          if (k.gt.1) dz = (height(k+1) - height(k-1))/2.  
c                            dtheta = (temp(b,k+1)*(p0/press(b,k+1))**gamma -
c     &                  temp(b,k)*(p0/press(b,k))**gamma)
c                            dtdz1d(k) = dtheta/dz
c                        else
c                            dtdz1d(k) = dtdz1d(k-1)
c                      endif
c                  enddo
c                  dstkabs = abs(dstk(n))
c                  call plumeris(nlay,hght1d,tempk1d,dtdz1d,wind1d,hstk(n),
c     &                  dstkabs,tstk(n),vstk(n),zstk)
c  14              continue
c                  do k = 1,10
c                      if (height(k).gt.zstk) goto 15
c                  enddo
c                  k = 10
c  15              continue
c
cC                ASSIGN EMISSIONS TO CORRECT CELL
c                      do nisp = 1,pisp
c                        if (isp.eq.pispid(nisp))
c     &                          pptemis(b,nisp,k) = pptemis(b,nisp,k) + ptemis1(n)
c                      enddo
c                      do nisp = 1,gisp
c                        if (isp.eq.gispid(nisp))
c     &                          gptemis(b,nisp,k) = gptemis(b,nisp,k) + ptemis1(n)
c                      enddo
c                  endif
c                  enddo
cc                enddo
c          end do  !isp
c              end do  !ihr

      close(IFIN)
      enddo

c Open area emissions file

c      aemis = indir(1:INDEX(indir,' ')-1)//'area.camx.9902.2001'//cdate(idate)//'.bin'
c      write(6,*)'INPUT FILE: ',aemis
c      open(unit=ifina, file=aemis, form='UNFORMATTED', status='old')

c      end do

 990   end
