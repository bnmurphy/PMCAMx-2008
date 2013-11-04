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

      integer ifin, ifout, ifout2, ifout3
      parameter(ifin = 10)
      parameter(ifout = 20)
      parameter(ifout2 = 60)
      parameter(ifout3 = 95)

      integer ny, nx, nt, nlay, nvar, ndate
      parameter(ny = 97) !columns, or 147
      parameter(nx = 88) !rows,    or 111
      parameter(nt = 24) !hours
      parameter(nlay = 1)
      parameter(nvar = 9) !METCRO2D files have 59 but we want only 9
      parameter(ndate = 32)
      parameter(nyr = 1)

      integer ihr, idate, day, ix, iy, icol, irow, itime, io, it, incr, io_edit, iyr
      integer m,d,i,j
      integer ncid, retval
      real var_in(ny,nx,ndate), var(nvar,nt,nx,ny,ndate), YYYYDDD(2*nt*59*ndate), t_in(nt,nx,ny)
      real sconh(ny,nx,ndate,nvar), ctot(nvar,nt)

      character*199 fname,fnameout,fnamedaily,isblk, fnametot
      character*198 indir,outdir,scenario 
      character*3   cdate(ndate)
      character*4   year,cyr(nyr)
      character*6   dum
      character*10  varname(nvar)
      character*2   varnum(nvar), io2

      data cdate /'182','183','184','185','186','187','188','189','190','191','192',
     &           '193','194','195','196','197','198','199','200',
     &           '201','202','203','204','205','206','207','208',
     &           '209','210','211','212','213'/

c      data cyr /'1991','1992','1993','1994','1995','1996','1997','1998','1999','2000'/
c      data cyr /'2051','2052','2053','2054','2055','2056','2057','2058','2059','2060'/
c      data cyr /'2051','2052','2054','2055','2056','2057','2058','2059','2060'/
      data cyr /'2053'/

      data varname /'TEMPG','TEMP2','GLW','GSW','RGRND','CFRAC','SNOCOV','VEG','LAI'/
c      data varnum /'11','12','16','17','18','21','25','26','27'/

c
c Define input/output locations
c
      indir = '/home/mday/mday/met_hd/MCIP_Output/readable'
      outdir = '/home/mday/mday/met_hd/MCIP_Output/readable'

c
c Open Output Files (daily files)
c
        print *,'Opening Output Files...','\n\n'
      do io = 1,nvar
        fnamedaily = outdir(1:INDEX(outdir,' ')-1)// '/' //
     &               varname(io)(1:INDEX(varname(io),' ')-1) // '_' // cyr(1)//'-'//cyr(nyr) // '.daily'
        open(unit=(ifout+io), file=fnamedaily(1:INDEX(fnamedaily,' ')-1),
     &                                status='UNKNOWN')
        write((ifout+io),91),'DAY','ROW','COL',varname(io)
        print *, 'Output file: ', fnamedaily
      end do

c
c   Initialize Output Species Matrix
c
        print *,'Initializing Output Concentration Arrays...'
        do m=1,nvar
         do d=1,ndate
          do i=1,nx
           do j=1,ny
             sconh(j,i,d,m) = 0.0
           end do
          end do
         end do
        end do

c Start year loop
      do iyr = 1,nyr
      year = cyr(iyr)  !'2060'
       print *, 'Starting ',year

c Open Input File
      do io = 1,nvar
        fname = indir(1:INDEX(indir,' ')-1)// '/' // varname(io)(1:INDEX(varname(io),' ')-1) 
     &          // '_' // year // '.daily'
        open(unit=(ifin+io), file=fname(1:INDEX(fname,' ')-1),status='UNKNOWN')
        read(ifin+io,*)dum,dum,dum,dum  !'DAY','ROW','COL',varname(io)
c        print *,io, dum
c        print *,'Input file: ',fname
      end do

c
c   Write to output file
c
      print *,'Reading concentrations...'
      do io = 1,nvar
       do idate = 1,ndate
          do ix = 1,nx   !rows
	    do iy = 1,ny !columns
              read((ifin+io),*),day,i,j,var_in(iy,ix,idate)
              if (ix.eq.1 .and. iy.eq.1 .and. idate.eq.1 .and. io.eq.1) then
              print *,'Check file start date: ', day, ' versus read date: ', cdate(idate)
              end if
              sconh(iy,ix,idate,io) = sconh(iy,ix,idate,io) + var_in(iy,ix,idate)
	    enddo !iy
          enddo   !ix
       enddo      !ndate
      enddo       !io

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
        do ix = 1,nx
          do iy = 1,ny
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
