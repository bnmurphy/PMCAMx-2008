      subroutine pigwalk(dt) 
c
c-----CAMx v4.02 030709
c
c     PIGWALK transports PiG puffs for the duration of one CG timestep
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        dt                  coarse grid time step (sec)     
c
c     Output arguments:
c        none
c
c     Subroutines called:
c        WALK1PIG
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "grid.com"
      include "filunit.com"
      include "camxfld.com"
      include "pigsty.com"
c
c-----Entry point
c
c-----Loop over all PiGs; find active puffs to move
c
      do 20 n = 1, npig
        if (ingrd(n).eq.0) goto 20
c
c-----New born PiGs walk for the duration of their age; set LNEWT = FALSE
c     Older puffs walk for the duration of the current timestep
c
        if (lnewt(n)) then
          dtcell = agepig(n)
        else
          dtcell = dt
          agepig(n) = agepig(n) + dt
        endif
c
c-----Walk the PiG
c
        kount = 0
        kpuff = 0
        igrdo = ingrd(n)
  10    kount = kount + 1
        igrd = ingrd(n)
        call walk1pig(kount,kpuff,dtcell,igrd,igrdo,n,ncol(igrd),
     &                nrow(igrd),nlay(igrd),windu(iptr3d(igrd)),
     &                windv(iptr3d(igrd)),height(iptr3d(igrd)))
        if (kount.gt.50 .and. dtcell.gt.0.) then
          write(iout,'(//,a)') 'ERROR in PIGWALK:'
          write(iout,*) 'Number of steps > 50'
          write(iout,*) 'Puff#,grid#,timestep left:'
          write(iout,*) n,igrd,dtcell
          call camxerr()
        endif
        if (dtcell.gt.0.) goto 10
  20  continue
c
      return
      end
