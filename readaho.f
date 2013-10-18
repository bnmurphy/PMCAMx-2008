      subroutine readaho(ncol,nrow,time,date,ly2k,name,ahotim,ahodate,
     &                   indx)
c
c-----CAMx v4.02 030709
c
c     READAHO reads the ozone or haze records of the albedo/haze/ozone
c     file to the current time/date
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        11/06/01  CAMx now assumes that all input file dates are in
c                  Julian format (YYJJJ) if the simulation year is 2000
c                  or greater.
c
c     Input arguments:
c        ncol                number of columns
c        nrow                number of rows
c        time                simulation time (HHMM)
c        date                simulation date (YYJJJ)
c        ly2k                Year 2000 flag (T is >=2000)
c        name                name of index to search (HAZE or OZONE COL)
c
c     Output arguments:
c        ahotim              next update time for AHO map (HHMM)
c        ahodate             next update date for AHO map (YYJJJ)
c        indx                coarse grid haze or albedo codes
c
c     Routines called:
c        JULDATE
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'filunit.com'
c     include 'grid.com'
c
      integer ahodate,date
      integer indx(ncol,nrow)
      character*10 title,namoz,namhz,name
      logical ly2k
c
      data namoz/'OZONE COL '/
      data namhz/'HAZE      '/
c
c-----Entry point
c
c-----Read coarse grid ozone column and haze indices
c
 100  read(iaho,'(a10,i10,f10.0,i10,f10.0)',end=900) title,id1,t1,id2,t2 
      if (title.ne.namoz .and. title.ne.namhz) then 
        write(iout,'(//,a)') 'ERROR in READAHO:'
        write(iout,*) 'Expecting HAZE or OZONE COL' 
        write(iout,*) 'Read from IAHO: ',title
        call camxerr()
      endif 
      if (.not.ly2k .and. id1.gt.100000) call juldate(id1)
      if (.not.ly2k .and. id2.gt.100000) call juldate(id2)
c
c-----Check title against specified name
c
      if (title.eq.name) then
        do j = nrow,1,-1
          read(iaho,'(999i1)') (indx(i,j),i=1,ncol)
        enddo
      else
        do j = nrow,1,-1
          read(iaho,'(999i1)') idum
        enddo
        goto 100
      endif
      write(iout,'(a40,2(f7.0,i8.5),1x,A9)')
     &       'Read albedo/haze/ozone file at ',t1,id1,t2,id2,name
c
c-----Check time/date
c
      if ((id1.lt.date .or. (id1.eq.date .and. t1.le.time)) .and.
     &    (id2.gt.date .or. (id2.eq.date .and. t2.gt.time))) then
        ahodate = id2
        ahotim = t2
        if (ahotim.ge.2400.) then
          ahotim = ahotim - 2400.
          ahodate = ahodate + 1
          if( MOD(ahodate,1000) .GT. 365 ) then
            if( MOD(INT(ahodate/1000),4) .EQ. 0 ) then
               if( MOD(ahodate,1000) .EQ. 367 )
     &                     ahodate = (INT(ahodate/1000)+1)*1000 + 1
            else
               ahodate = (INT(ahodate/1000)+1)*1000 + 1
            endif
         endif
        endif
        goto 999
      else
        goto 100
      endif
c
c-----Reached end of AHO file; rewind and read header/albedo records again
c
 900  continue
      rewind(iaho)
      do n = 1,4
        read(iaho,*)
      enddo
      read(iaho,'(10x,3i10)') idd,nx,ny
      do n = 1,ny
        read(iaho,*)
      enddo 
 901  read(iaho,'(10x,3i10)') idd,nx,ny
      if (idd.ne.0) then
        do n = 1,ny
          read(iaho,*)
        enddo 
        goto 901
      endif
      goto 100
c
 999  return
      end
