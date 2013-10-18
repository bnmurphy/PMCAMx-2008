      subroutine getdepth(ncols,nrows,nlays,ly2k,ibgdhp,idtnow,btimhp,
     &                    timnow,iunit,height,depth)
c
c-----CAMx v4.02 030709
c
c     GETDEPTH reads the height/pressure file until to the current hour,
c     calculates layer depth, and then save the data for OSAT calculations
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        11/06/01  Added Y2K flag
c 
c     Input arguments: 
c        ncols               number of columns
c        nrows               number of rows
c        nlays               number of layers
c        l2yk                Year 2000 flag (T is >=2000)
c        idtnow              current date (YYJJJ)
c        timnow              current time (HHMM)
c        iunit                height/pressure file unit number
c             
c     Output arguments: 
c        ibgdhp              date from height/pressure file (YYJJJ)
c        btimhp              time from height/pressure file (HHMM) 
c        height              array of gridded layer heights (m)
c        depth               array of gridded layer depths (m)
c             
c     Routines Called: 
c        JULDATE
c             
c     Called by: 
c        CLCIWT
c        RDSUMBC
c
      include "camx.prm"
      include "filunit.com"
c
      real height(ncols,nrows,nlays),depth(ncols,nrows,nlays)
      logical ly2k
c
c-----Read Height pressure file to current time/date
c
  333 continue
      if( (ibgdhp .EQ. idtnow .AND. btimhp .LT. timnow) .OR.
     &     ibgdhp .LT. idtnow ) then
         do 20 k=1,nlays
           read(iunit,ERR=7001) btimhp, ibgdhp, 
     &                         ((height(i,j,k),i=1,ncols),j=1,nrows)
           read(iunit,ERR=7001) idum,idum,((dum,i=1,ncols),j=1,nrows)
   20    continue
         if (.not.ly2k .and. ibgdhp .GT. 100000) call juldate(ibgdhp)
         btimhp = btimhp/100.0
         if( idtnow .EQ. ibgdhp .AND. timnow .GT. btimhp ) goto 333
         if( idtnow .LT. ibgdhp ) goto 333
c
c-----Calculate the depths from the heights
c
         do 30 k=1,nlays
            do 40 j=1,nrows
               do 50 i=1,ncols
                  depth(i,j,k) = height(i,j,k)
                  if( k .GT. 1 ) 
     &              depth(i,j,k) = height(i,j,k) - height(i,j,k-1)
   50           continue
   40       continue
   30    continue
      endif
c
      return
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in GETDEPTH:'
      write(iout,'(/,1X,A)') 'Reading height/pressure file'
      call camxerr()
c
      end
