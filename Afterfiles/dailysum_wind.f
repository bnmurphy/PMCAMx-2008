      program mgnout2camx
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c  mgnout2camx reads from the MEGANv2.04 Output files.
c  It then outputs data for the CAMx domain.
c
c  MDay 5/2011
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c      implicit none

c    netCDF ID for file and data variable
c      include 'netcdf.inc'   !from $MGNLIB/netcdf-3.6.3/fortran/netcdf.inc ?

      integer ifin, ifout, ifout2
      parameter(ifin   = 10)
      parameter(ifout  = 20)
      parameteR(ifout2 = 90)

      integer ny, nx, nlay, nvar, ndate
      parameter(ny = 97) !columns, or 147
      parameter(nx = 90) !rows,    or 111
      parameter(nlay = 1)
      parameter(nvar = 1)!filenames you want to run
      parameter(ndate = 31.)
      parameter(nyr = 1.)

      integer ihr, idate, day, ix, iy, icol, irow, itime, io, it, incr, io_edit, iyr
      integer m,d,i,j
      integer ncid, retval
      real var_in(ny,nx,ndate) !, var(nvar,nt,nx,ny,ndate), YYYYDDD(2*nt*59*ndate), t_in(nt,nx,ny)
      real uvar(ny,nx,ndate), vvar(ny,nx,ndate)
      real sconh(ny,nx,ndate,nvar), ctot(ny,nx,nvar)

      character*199 fname,fnameout,fnamedaily,isblk, fnametot
      character*198 indir,outdir,scenario 
      character*3   cdate(ndate)
      character*4   year,cyr(nyr)
      character*6   dum
      character*8   add
      character*10  varname(nvar)
      character*2   varnum(nvar), io2

      data cdate /'182','183','184','185','186','187','188','189','190','191','192',
     &           '193','194','195','196','197','198','199','200','201','202','203',
     &           '204','205','206','207','208','209','210','211','212'/ !,'213'/

      data cyr /'1991'/
c      data cyr /'1991','1992','1993','1994','1995','1996','1997','1998','1999','2000'/
c      data cyr /'2051','2052','2053','2054','2055','2056','2057','2058','2059','2060'/
c      data cyr /'2051','2052','2054','2055','2056','2057','2058','2059','2060'/

      data varname /'uvwind'/
c      data varname /'TEMPG','TEMP2','GLW','GSW','RGRND','CFRAC','SNOCOV','VEG','LAI'/
c      data varnum /'11','12','16','17','18','21','25','26','27'/

c
c Define input/output locations
c
      indir = '/home/mday/PMCAMx/met_hd/CAMx_Input/readable/'
      outdir = '/home/mday/PMCAMx/met_hd/CAMx_Input/readable/'

c
c Open Output Files (daily files)
c
        print *,'Opening Output Files...','\n\n'
      do io = 1,nvar
        fnamedaily = outdir(1:INDEX(outdir,' ')-1)//
     &               varname(io)(1:INDEX(varname(io),' ')-1) // '_' // cyr(1)//'-'//cyr(nyr) // '.daily'
        open(unit=(ifout+io), file=fnamedaily(1:INDEX(fnamedaily,' ')-1),status='UNKNOWN')
        write((ifout+io),91),'DAY','ROW','COL',varname(io)
        print *, 'Output file: ', fnamedaily
      end do

c
c   Initialize Output Species Matrix
c
        print *,'Initializing Output Concentration Arrays...'
        do m=1,nvar
         do d=1,ndate
          do j=1,ny
           do i=1,nx
             sconh(j,i,d,m) = 0.0
           end do
          end do
         end do
        end do

c Start year loop
      do iyr = 1,nyr
      year = cyr(iyr)  !'2060'
       print *, 'Starting ',year
      add = '.Jul' // year

c Open Input File
      do io = 1,nvar
        fname = indir(1:INDEX(indir,' ')-1) //
     &          varname(io)(1:INDEX(varname(io),' ')-1) // add // '.daily.01'
        open(unit=(ifin+io), file=fname(1:INDEX(fname,' ')-1),status='OLD')
        read(ifin+io,*)dum,dum,dum,dum,dum  !'DAY','ROW','COL',uwind,vwind
c        print *,io, dum
        print *,'Input file: ',fname

c Open decadal output file
        fnametot = outdir(1:INDEX(outdir,' ')-1) // varname(io)(1:INDEX(varname(io),' ')-1)//
     &            '_' // year // '.tot'
        open(unit=(ifout2+io),file=fnametot(1:INDEX(fnametot,' ')-1), status='UNKNOWN')
         print *, 'Output file: ',fnametot
      end do

c   Initialize tot output file
        print *,'Initializing Yearly Output Concentration Arrays...'
        do m=1,nvar
          do j=1,ny
           do i=1,nx
             ctot(j,i,m) = 0.0
           end do
          end do
        end do

c
c   Write to output file
c
      print *,'Reading concentrations...'
      do io = 1,nvar
       do idate = 1,ndate
          do iy = 1,ny   !columns
	    do ix = 1,nx !rows
              read((ifin+io),*),day,i,j,uvar(iy,ix,idate),vvar(iy,ix,idate)
              if (ix.eq.1 .and. iy.eq.1 .and. idate.eq.1 .and. io.eq.1) then
              print *,'Check file start date: ', day, ' versus read date: ', cdate(idate)
c              elseif (ix.ge.(nx-1) .and. iy.ge.(ny-1)) then
c              print *,cdate(idate),day,iy,i,ix,j,var_in(iy,ix,idate)
              end if
              var_in(iy,ix,idate) = sqrt(uvar(iy,ix,idate)**2.+vvar(iy,ix,idate)**2.) !compute magnitude
              sconh(iy,ix,idate,io) = sconh(iy,ix,idate,io) + var_in(iy,ix,idate)
              ctot(iy,ix,io) = ctot(iy,ix,io) + var_in(iy,ix,idate)
              if (ix.eq.1 .and. iy.eq.1 .and. io.eq.1) then
              print *,iy,ix,uvar(iy,ix,idate),vvar(iy,ix,idate),ctot(iy,ix,io)
              endif
	    enddo !ix
          enddo   !iy
       enddo      !ndate
      enddo       !io

       print *, 'Writing output files...','\n'
      do io = 1,nvar
        do iy = 1,ny
          do ix = 1,nx
            write((ifout2+io),*),ix,iy,ctot(iy,ix,io)/ndate
          enddo
        enddo
      enddo  !io

c
c   Close the file, freeing all resources
c
      do io = 1,nvar
        close ((ifin+io))
      end do

c  
c  Done with Year -> Start Next One
c
       print *,'Done with ',cyr(iyr),'\n\n'
      enddo !Year -> iyr iterates

c Final concentration average!
       print *, 'Writing output files...','\n'
      do io = 1,nvar
       do idate = 1,ndate
        do iy = 1,ny
          do ix = 1,nx
            write((ifout+io),*),cdate(idate),ix,iy,sconh(iy,ix,idate,io)/nyr
          enddo
        enddo        
       end do
      enddo  !io
        
c       
c  Close Output Files
c       
       print *,'Closing Output Files...','\n'
       do io = 1,nvar
          close ((io+ifout))      !Close Output Files
       end do


 90   format(A3,10x, A4,9x, A3,9x, A3,2x, A10)
c 91   format(A3,9x, A3,9x, A10)
 91   format(A3,10x, A3,9x, A3,9x, A10)
c 91   format(I3,2x, I2,2x, I3,2x, I3,2x, E10.4)

      end 
