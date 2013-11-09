      subroutine readbnd(bndtim,bnddate)
c 
c-----CAMx v4.02 030709
c 
c     READBND reads and cycles through the BOUNDARY file to the
c     current time/date, and loads coarse grid boundary concentrations
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none
c 
c     Input arguments: 
c        bndtim                 model simulation time (HHMM)
c        bnddate                model simulation date (YYJJJ)
c             
c     Output arguments: 
c        bndtim                 next boundary update time (HHMM)
c        bnddate                next boundary update date (YYJJJ)
c             
c     Routines Called: 
c        none
c             
c     Called by: 
c        CAMx
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'bndary.com'
      include 'grid.com'
      include 'chmstry.com'
c
      character*4 bcspec(10)
      integer bnddate
      real bctmp(MX1D,MXLAYA,4,MXSPEC)
c
c-----Entry point
c
      nz = nlay(1)
c
c-----Read through coarse grid concentration records until current time/date
c
 100  read(ibc,end=900) idat1,tim1,idat2,tim2
      if (INT(tim2) .EQ. 24) then
        tim2 = 0.
        idat2 = idat2 + 1
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      do l = 1,nbcspc
        do n = 1,4
          nc = nrow(1)
          if (n.gt.2) nc = ncol(1)
         read(ibc) idum,(bcspec(j),j=1,10),iedge,
     &              ((bctmp(i,k,n,l),k=1,nz),i=1,nc)
c BNM             read(ibc) (bcspec(j),j=1,10),iedge,
c     &              ((bctmp(i,k,n,l),k=1,nz),i=1,nc)
	enddo
      enddo
c

      write(iout,'(a40,2(f7.0,i8.5))')
     &  'Read boundary condition file at ',tim1,idat1,tim2,idat2
      if ((idat1.lt.date .or. (idat1.eq.date .and. tim1.le.time)) .and.
     &    (idat2.gt.date .or. (idat2.eq.date .and. tim2.gt.time))) then
c
c-----Load boundary concentrations; convert gasses from ppm to umol/m3,
c     PM stays at ug/m3
c
        do 60 l = 1,nspec
          lbc = lbcmap(l)
          do 50 k = 1,nz
c
            do 40 j = 1,nrow(1)
              if (ibeg(j).eq.-999) goto 40
              i = ibeg(j) - 1
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(l - 1)
              if (l.le.ngas) then
                convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              else
                convfac = 1.
              endif
              conc(n4d) = convfac*bdnl(l)
              if (lbc.gt.0 .and. bctmp(j,k,1,lbc).gt.bdnl(l)) 
     &          conc(n4d) = convfac*bctmp(j,k,1,lbc)
c-----Set POC and PEC BC to zero-----------------------------------
c              if (l.ge.548.and.l.le.553) then
c                conc(n4d) = 0.18 
c              endif
c              if (l.ge.118.and.l.le.123) then
c                conc(n4d) = 0.0
c              endif
c------------------------------------------------------------------              
c
              i = iend(j) + 1
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1)
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(l - 1)
              if (l.le.ngas) then
                convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              else
                convfac = 1.
              endif
              conc(n4d) = convfac*bdnl(l)
              if (lbc.gt.0 .and. bctmp(j,k,2,lbc).gt.bdnl(l)) 
     &          conc(n4d) = convfac*bctmp(j,k,2,lbc)
c-----Set POC and PEC BC to zero-----------------------------------
c              if (l.ge.548.and.l.le.553) then
c                conc(n4d) = 0.18
c              endif
c              if (l.ge.118.and.l.le.123) then
c                conc(n4d) = 0.0
c              endif 
c------------------------------------------------------------------

 40         continue
c
            do 45 i = 1,ncol(1) 
              if (jbeg(i).eq.-999) goto 45
              j = jbeg(i) - 1 
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(l - 1)
              if (l.le.ngas) then
                convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              else
                convfac = 1.
              endif
              conc(n4d) = convfac*bdnl(l)
              if (lbc.gt.0 .and. bctmp(i,k,3,lbc).gt.bdnl(l))
     &          conc(n4d) = convfac*bctmp(i,k,3,lbc)
c-----Set POC and PEC BC to zero-----------------------------------
c              if (l.ge.548.and.l.le.553) then
c                conc(n4d) = 0.18
c              endif
c              if (l.ge.118.and.l.le.123) then
c                conc(n4d) = 0.0
c              endif
c------------------------------------------------------------------

c
              j = jend(i) + 1 
              n3d = i + ncol(1)*(j - 1) + ncol(1)*nrow(1)*(k - 1) 
              n4d = n3d + ncol(1)*nrow(1)*nlay(1)*(l - 1)
              if (l.le.ngas) then
                convfac = densfac*273./tempk(n3d)*press(n3d)/1013.
              else
                convfac = 1.
              endif
              conc(n4d) = convfac*bdnl(l)
              if (lbc.gt.0 .and. bctmp(i,k,4,lbc).gt.bdnl(l))
     &            conc(n4d) = convfac*bctmp(i,k,4,lbc) 
c-----Set POC and PEC BC to zero-----------------------------------
c              if (l.ge.548.and.l.le.553) then
c                conc(n4d) = 0.18
c              endif
c              if (l.ge.118.and.l.le.123) then
c                conc(n4d) = 0.0
c              endif
c------------------------------------------------------------------

c                    if (i.eq.5.and.k.eq.1) then
c                      write(6,*)lbc,l,conc(n4d)
c                    endif
 45         continue 
 50       continue
 60     continue
c        pause
      else
        goto 100
      endif
c
c-----Set next boundary update time
c
      bndtim = tim2
      bnddate = idat2
      if (bndtim.ge.2400.) then
        bndtim = bndtim - 2400.
        bnddate = bnddate + 1
        if( MOD(bnddate,1000) .GT. 365 ) then
           if( MOD(INT(bnddate/1000),4) .EQ. 0 ) then
              if( MOD(bnddate,1000) .EQ. 367 )
     &                     bnddate = (INT(bnddate/1000)+1)*1000 + 1
           else
              bnddate = (INT(bnddate/1000)+1)*1000 + 1
           endif
        endif
      endif
      goto 999
c
c-----End of BC file reached
c
 900  write(iout,'(//,a)') 'ERROR in READBND:'
      write(iout,*)'Premature End of BC file reached.'
      write(iout,*)'Make sure boundary file contains simulation ',
     &                                               'time period.'
      call camxerr()
c
 999  return
      end
