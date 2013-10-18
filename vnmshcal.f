      subroutine vnmshcal(igrd,ncol,nrow,nlay,i1,j1,ncolf,nrowf,nlayf,
     &                    height,heightf,nmshv)
c
c-----CAMx v4.02 030709
c
c     VMESHCAL determines the vertical mesh number using coarse and
c     fine height data
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c 
c     Modifications:  
c        1/28/99   Changed check on fine/coarse layer interface matching
c                  from absolute (0.1 m) to relative (5%)
c  
c     Input arguments:  
c        igrd                grid number of fine grid
c        ncol                number of columns in coarse grid
c        nrow                number of rows in coarse grid
c        nlay                number of layers on coarse grid
c        i1                  starting i index for the fine grid 
c        j1                  starting j index for the fine grid 
c        ncolf               number of columns in fine grid
c        nrowf               number of rows in fine grid
c        nlayf               number of layers on fine grid
c        height              layer interface height on coarse grid(m)
c        heightf             layer interface height on fine grid(m)
c             
c     Output arguments:  
c        nmshv               vertical fine mesh number by coarse layer
c
c     Routines called:
c        none
c
c     Called by:
c        STARTUP 
c
      include "camx.prm"
      include "filunit.com"
c
      dimension height(ncol,nrow,nlay),heightf(ncolf,nrowf,nlayf),
     &          nmshv(nlay)
c
c------Entry point
c
c
c------If ZP file was not supplied, set meshign factor to 1 ----
c
      if( ihtp(igrd) .LE. 0 ) then
         do kp=1,nlay
              nmshv(kp) = 1
         enddo
         goto 9999
      endif
c
c------Otherwise figure it out---
c
      kg1 = 1
      do 20 kp = 1,nlay
        ntmp = 0
        do kg = kg1,nlayf
          ntmp = ntmp + 1
          dht = abs(heightf(2,2,kg) - height(i1,j1,kp))/height(i1,j1,kp)
          if (dht.lt.0.05) then
            nmshv(kp) = ntmp
            kg1 = kg + 1
            goto 20
          endif
        enddo
  20  continue
c
      do 30 kp = 1,nlay
        if (nmshv(kp).lt.1) then
          write(iout,'(//,a)') 'ERROR in VNMSHCAL:'
          write(iout,*) 'Vertical meshing factor < 1!'
          write(iout,*) 'In layer ',kp
          call camxerr()
        endif
 30   continue
c
 9999 continue
      return
      end
