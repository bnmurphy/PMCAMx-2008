      subroutine finwind(iunit,ncol,nrow,nlay,io,jo,nmesh,nmshv,ncolf,
     &       nrowf,nlayf,windu,windv,winduf,windvf,igrd,date,time,iout)
c
c-----CAMx v4.02 030709
c
c     FINWIND estimates wind speed for a fine grid from its parent
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        11/06/01  Input dates are now Julian
c
c     Input arguments:
c        iunit             child grid wind file unit number
c        ncol              number of columns on parent grid
c        nrow              number of rows on parent grid
c        nlay              number of layers on parent grid
c        io,jo             parent starting grid cell indices
c        nmesh             fine grid meshing factor
c        nmshv             fine grid vertical meshing factor
c        ncolf             number of columns on fine grid
c        nrowf             number of rows on fine grid
c        nlayf             number of layers on fine grid
c        windu             u-component wind speed on parent grid (m/s)
c        windv             v-component wind speed on parent grid (m/s)
c        igrd              current grid number
c        date              current date (YYJJJ)
c        time              current time
c        iout              unit number of .out file
c
c     Output arguments:
c        windu             u-component wind speed on fine grid (m/s)
c        windv             v-component wind speed on fine grid (m/s)
c
c     Subroutine called:
C        none
c
c     Called by:
c        STARTUP
c        INTRPDAT
c
      dimension windu(ncol,nrow,nlay),windv(ncol,nrow,nlay),
     &          winduf(ncolf,nrowf,nlayf),windvf(ncolf,nrowf,nlayf),
     &          nmshv(nlay)
      integer   date
c
c-----Entry point
c
c-----Assign wind speed to the fine grid: mass is conserved
c
      kg1 = 1
      do 50 kp = 1,nlay
        do 40 kg = kg1,kg1+nmshv(kp)-1
c
c-----In the case that fine grid winds are supplied via input file,
c     set only the fine boundary winds from parent winds
c
          if (iunit.ne.0) then
            do j = 2,nrowf-1
              ii = io - 1
              jj = jo + (j-2)/nmesh
              winduf(1,j,kg) = windu(ii,jj,kp)
              ii = io - 1 + (ncolf-2)/nmesh
              winduf(ncolf-1,j,kg) = windu(ii,jj,kp)
            enddo
            do i = 2,ncolf-1
              ii = io + (i-2)/nmesh
              jj = jo - 1
              windvf(i,1,kg) = windv(ii,jj,kp)
              jj = jo - 1 + (nrowf-2)/nmesh
              windvf(i,nrowf-1,kg) = windv(ii,jj,kp)
            enddo
            goto 40
          endif
c
c-----Use the wind where it available from the parent grid
c
          if( kp .EQ. 1 ) write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating winds from parent grid',
     &                             time, date,' grid',igrd
          do j = 2,nrowf-1
            do i = 1,ncolf-1,nmesh
              ii = io - 1 + (i-1)/nmesh
              jj = jo + (j-2)/nmesh
              winduf(i,j,kg) = windu(ii,jj,kp)
            enddo
          enddo
c
          do j = 1,nrowf-1,nmesh
            do i = 2,ncolf-1
              ii = io + (i-2)/nmesh
              jj = jo - 1 + (j-1)/nmesh
              windvf(i,j,kg) = windv(ii,jj,kp)
            enddo
          enddo
c
c-----Interpolate in between
c
          do j = 2,nrowf-1
            do i = 1,ncolf-nmesh,nmesh
              du = (winduf(i+nmesh,j,kg)-winduf(i,j,kg))/nmesh
              do i1 = i+1,i+nmesh-1
                winduf(i1,j,kg) = winduf(i1-1,j,kg) + du
              enddo
            enddo
          enddo
c
          do j = 1,nrowf-nmesh,nmesh
            do i = 2,ncolf-1
              dv = (windvf(i,j+nmesh,kg)-windvf(i,j,kg))/nmesh
              do j1 = j+1,j+nmesh-1
                windvf(i,j1,kg) = windvf(i,j1-1,kg) + dv
              enddo
            enddo
          enddo
c
  40    continue
        kg1 = kg
  50  continue
c
      return
      end
