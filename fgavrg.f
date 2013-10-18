      subroutine fgavrg(igrd)
c
c-----CAMx v4.02 030709
c
c     FGAVRG passes arrays from common blocks to AVERAGE
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c  
c     Modifications:
c        none
c
c     Input arguments: 
c        igrd                current grid index
c
c     Output arguments: 
c        none
c      
c     Routines called: 
c        AVERAGE
c      
c     Called by: 
c        NESTING
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'grid.com'
      include 'chmstry.com'
      include 'flags.com'
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
c
c========================= Source Apportion End ========================
c
c
c======================== Process Analysis Begin ====================================
c
      include 'procan.com'
c
c========================= Process Analysis End =====================================
c
c-----Entry point
c
c-----Update cumulative time for this fine grid 
c
      call average(.FALSE.,igrd,deltat(igrd)/2.0,
     &             ncol(igrd),nrow(igrd),nlay(igrd),nlay(igrd),
     &             navspc,nspec,lavmap,tempk(iptr3d(igrd)),
     &             press(iptr3d(igrd)),
     &             conc(iptr4d(igrd)),avcnc(iptr4d(igrd)),
     &             ipacl_3d(iptr3d(igrd)) )
c
c---- update PiG contribution
c
      if( ipigflg .EQ. IRONPIG ) then
         call avepig(igrd,deltat(igrd)/2.0,ncol(igrd),nrow(igrd),
     &             nlay(igrd),nlay(igrd),deltax(1,igrd),deltay(igrd),
     &             mapscl(iptr2d(igrd)), height(iptr3d(igrd)),navspc,
     &             nspec,lavmap,tempk(iptr3d(igrd)),press(iptr3d(igrd)),
     &             conc(iptr4d(igrd)),avcnc(iptr4d(igrd)))
      endif
c
c======================== Source Apportion Begin =======================
c
c   --- call routine to update the running averages ---
c
      if( ltrace .OR. lddm ) then
         call average(.TRUE.,igrd,deltat(igrd)/2.0,
     &                ncol(igrd),nrow(igrd),nlay(igrd),1,
     &                ntotsp,ntotsp,lsamap,tempk(iptr3d(igrd)),
     &                press(iptr3d(igrd)),ptconc(ipsa3d(igrd)),
     &                ptavrg(ipsa2d(igrd)),
     &               ipacl_3d(iptr3d(igrd)) )
      endif
c
c========================= Source Apportion End ========================
c
      return
      end
