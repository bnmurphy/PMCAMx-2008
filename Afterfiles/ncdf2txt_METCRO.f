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
      include 'netcdf.inc'   !from $MGNLIB/netcdf-3.6.3/fortran/netcdf.inc ?

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

      integer ihr, idate, ix, iy, icol, irow, itime, io, it, incr, io_edit
      integer ncid, retval
      real var_in(ny,nx,nt,ndate), var(nvar,nt,nx,ny,ndate), YYYYDDD(2*nt*59*ndate), t_in(nt,nx,ny)
      real sconh(ny,nx,nvar,ndate), ctot(nvar,nt)

      character*199 fname,fnameout,fnamedaily,isblk, fnametot
      character*198 indir,outdir,scenario 
      character*3   cdate(ndate)
      character*4   year
      character*10  varname(nvar)
      character*2   varnum(nvar), io2

      data cdate /'182','183','184','185','186','187','188','189','190','191','192',
     &           '193','194','195','196','197','198','199','200',
     &           '201','202','203','204','205','206','207','208',
     &           '209','210','211','212','213'/

      data varname /'TEMPG','TEMP2','GLW','GSW','RGRND','CFRAC','SNOCOV','VEG','LAI'/
      data varnum /'11','12','16','17','18','21','25','26','27'/

c
c Define input/output locations
c
      indir = '/home/mday/PMCAMx/met_hd/MCIP_Output' !'/home2/mday/MCIP_out'
      outdir = '/home/mday/PMCAMx/met_hd/MCIP_Output/readable'
      scenario = '53b'
      year = '2053'

c
c Open Output Files (full dataset and daily files)
c
        print *,'Opening Output Files...','\n\n'
      do io = 1,nvar
        fnameout = outdir(1:INDEX(outdir,' ')-1)// '/' //
     &             varname(io)(1:INDEX(varname(io),' ')-1) // '_' // year // '.txt'
        open(unit=(ifout+io), file=fnameout(1:INDEX(fnameout,' ')-1),
     &                              status='UNKNOWN')
        write((ifout+io),90),'DAY','TIME','ROW','COL',varname(io)
c        print *, 'Total data out file: ',fnameout

        fnamedaily = outdir(1:INDEX(outdir,' ')-1)// '/' //
     &               varname(io)(1:INDEX(varname(io),' ')-1) // '_' // year // '.daily'
        open(unit=(ifout2+io), file=fnamedaily(1:INDEX(fnamedaily,' ')-1),
     &                                status='UNKNOWN')
        write((ifout2+io),91),'DAY','ROW','COL',varname(io)
c        print *, 'Daily data out file: ', fnamedaily

        fnametot =  outdir(1:INDEX(outdir,' ')-1)// '/' // 
     &               varname(io)(1:INDEX(varname(io),' ')-1) // '_' // year // '.tot'
        open(unit=(ifout3+io),file=fnametot(1:INDEX(fnametot,' ')-1),
     &                                status='UNKNOWN')
c        print *, 'Grid sum data out file: ',fnametot
      end do

c      do idate = 1,ndate
c      print *,'Day ',cdate(idate)

c
c Open Input File
c
       fname = indir(1:INDEX(indir,' ')-1)// '/' // 'METCRO2D_JUL2_' // scenario(1:INDEX(scenario,' ')-1)
       print *,'Input file: ',fname
c       open(unit=ifin, file=fname, status='old')

c
c   Open the netCDF file with read-only access
c
       retval = nf_open(fname, NF_NOWRITE, ncid)
       retval = nf_get_var_real(ncid,1,YYYYDDD)
        print *, YYYYDDD(1), YYYYDDD(2)

c           
c   Initialize Output Species Matrix
c              
        print *,'Initializing Output Concentration Arrays...'
        do m=1,nvar  
         do idate=1,ndate
           do it=1,nt
             ctot(m,it) = 0.0
           enddo
           do i=1,nx
              do j=1,ny 
                 sconh(j,i,m,idate) = 0.0
              end do    
           end do    
         end do   
        end do     

c
c   Write to output file
c
      print *,'Reading concentrations and writing output files...'
      do io = 1,nvar
        io2 = varnum(io)
        read(io2,'(i2)') io_edit
        print *,'Compound: ',varname(io), ' number: ',io_edit
         retval = nf_get_var_real(ncid,(io_edit+1),var_in)
       do idate = 1,ndate
        do it = 1,nt
          do ix = 1,nx   !rows
	    do iy = 1,ny !columns
              write((ifout+io),*),cdate(idate),it,ix,iy,var_in(iy,ix,it,idate)
              sconh(iy,ix,io,idate) = sconh(iy,ix,io,idate) + var_in(iy,ix,it,idate)
c              ctot(io,it) = ctot(io,it) + sconh(iy,ix,io,idate)
	    enddo !iy
          enddo   !ix
        enddo     !it
       enddo      !ndate
c      enddo       !io
c      print *,ctot

c      do io = 1,nvar
       do idate = 1,ndate
        do ix = 1,nx
          do iy = 1,ny
            write((ifout2+io),*),cdate(idate),ix,iy,sconh(iy,ix,io,idate)/nt
          enddo
        enddo

c        do it = 1,nt
c          write((ifout3+io),*) cdate(idate),it,ctot(io,it)
c          print *, ctot(io,it)
c        enddo
       enddo !ndate
      enddo  !io

c
c   Close the file, freeing all resources
c
      retval = nf_close(ncid)

c  
c  Done with Day -> Start Next One
c
c       print *,'Done with Day ',cdate(idate),'\n\n'
c      enddo !Day -> idate iterates
        
c       
c  Close Output Files
c       
       print *,'Closing Output Files...','\n'
       do io = 1,nvar
          close ((io+ifout))      !Close Output Files
          close ((io+ifout2))
          close ((io+ifout3))
       end do


 90   format(A3,10x, A4,9x, A3,9x, A3,2x, A10)
 91   format(A3,10x, A3,9x, A3,9x, A10)
c 91   format(I3,2x, I2,2x, I3,2x, I3,2x, E10.4)

      end 

