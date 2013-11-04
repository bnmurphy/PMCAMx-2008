      SUBROUTINE  after_saprc(outdir,scenario,nlay,clay,ilay,ndate,
     &                  cdate,mspc,cspc,nspc,mconc)
c######################################################################
c   BNM
c     This program processes PMCAMx results
c     The parameters are hardcoded for the specific output file
c     Ben Murphy rewrote this post-processor 3-11-08
c       outputs Gas, PM 2.5, and NOT(PM 10) averaged over each day
c
c     Modifications
c    nmol to ppb conversion for CPO,COO,CBS,CAS implemented 7-1-08
c
c######################################################################

c
c Define constants and arrays
c

      parameter(ifin  = 10)
      parameter(ifout = 30)

      parameter(ncol = 97)
      parameter(nrow = 90)
      parameter(nhr  = 24)    ! Hours in each day
c      parameter(mspc = 95)    ! Number of output species
      parameter(gasfac = 1000)

      parameter(R    = 83.145)  !cm^3*bar/mol/K
      parameter(T    = 300)     !Kelvin
      parameter(P    = 1.013)   !bar

      integer MW
      integer ndate, mspc, nspc
      integer nlay, ilay
      character*2 clay(nlay)

      character*198 fname1
      character*198 outdir, indir, indir1
      character*190 fname
      character*198 scenario   !the name of the output folder
      character*198 scenario2  !name of input file

      character*4   name(10),note(60),mspec(10,nspc)
      character*2   chr(nhr)         !Hour
      character*3   cdate(ndate)   !Julian Date
      character*10  spcname(nspc)    !Input species name string
      character*7   cspc(mspc)      !Output species name vector

      real*4 maxh(mspc)
      real*4 btime,etime,rdum,xorg,yorg,delx,dely
      real*4 tbtime
      real*4 sconh(mspc,ncol,nrow)            !Avg conc to be output
      real*4 mconc(ndate,31,ncol,nrow,mspc)  !Conc read from input file
      real*4 date1
      
      integer*4 ione,nspec,ibdate,iedate,idum,izero,nx,ny,nz
      integer*4 itbdate,col,row
      integer*4 ihr,msec,icut,line
      integer*4 iday,ix,iy,iz

      data chr /'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &          '21','22','23','24'/

c
c Open Output Files
c
      print *,'\n','Opening Output Files...','\n\n'
      do io = 1,mspc
      isblk = INDEX(cspc(io),' ')
          if(isblk.eq.0) isblk = 5
          fname = outdir(1:INDEX(outdir,' ')-1)//
     &           scenario(1:INDEX(scenario,' ')-1)//'/'//
     &           cspc(io)(1:isblk-1)//'.daily.'//clay(ilay)
        
      isblk = INDEX(fname,' ')
      open(unit=(ifout+io),file=fname(1:isblk-1),status='UNKNOWN')
          write(ifout+io, 46), 'Date', 'X','Y','CONC'

      
      end do
c
c Iterate Days (input files)
c
      do idate=1,ndate
       print *,'After_SAPRC subroutine says: idate = ',idate,'  ndate = ',ndate
       print *
       print *,'!! Processing Daily Average Files for Maps (Task 1)!!'
       print *

c
c   Add up daily averages
c
      do io = 1,mspc
        do i = 1,ncol
        do j = 1,nrow
            sconh(io,i,j) = 0.0
            do k = 1,nhr
            sconh(io,i,j) = sconh(io,i,j) + mconc(idate,k,i,j,io)
            end do
        end do
        end do
      end do

c
c  Write to Output Files
c
      print *,'Exporting Data to directory: ',outdir(1:INDEX(outdir,' ')-1)
      do io = 1,mspc
        print *,'sconh: ',sconh(io,35,46)

        do i=1,ncol   ! COLUMN FIRST
      do j=1,nrow     ! ROW LATER
            write( (ifout+io), 47),cdate(idate),i,j,sconh(io,i,j)/nhr
          end do
        end do

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
        do io = 1,mspc
        close ((ifout+io))    !Close Output Files
        end do

  46   format(A4,2x,A1,3x,A1,3x,A4)
  47   format(A3,3x,i2,2x,i2,2x,e12.5)

       write(6,*)'ALL DAYS DONE'

 990   end
