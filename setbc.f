      subroutine setbc(ip)
c
c-----CAMx v4.02 030709
c
c     SETBC sets up boundary conditions for children grids
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        01/22/02     Now only sets OSAT to lower bound if not in RTRAC
c        01/30/02     Added code for RTRAC probing tool
c
c     Input arguments:
c        ip                  parent grid index 
c
c     Output arguments:
c        none
c
c     Routines Called:
c        BC1GRD
c        BCMODFY
c
c     Called By:
c        NESTING 
c
      include "camx.prm"
      include "grid.com"
      include "camxfld.com"
      include "chmstry.com"
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
c
c========================= Source Apportion End ========================
c
c-----Entry point
c
c-----When children grids are not connected
c
      do ic = 1,nchdrn(ip)
        ig = idchdrn(ic,ip)
        call bc1grd(nspec,ncol(ip),nrow(ip),nlay(ip),ncol(ig),nrow(ig),
     &              nlay(ig),i1(ig),j1(ig),nmesh(ig),nmshv(1,ig),
     &              conc(iptr4d(ip)),conc(iptr4d(ig)) )
c
c======================== Source Apportion Begin =======================
c
        if( ltrace .OR. lddm ) then
            call bc1grd(ntotsp,ncol(ip),nrow(ip),nlay(ip),ncol(ig),
     &              nrow(ig),nlay(ig),i1(ig),j1(ig),nmesh(ig),
     &              nmshv(1,ig),ptconc(ipsa3d(ip)),ptconc(ipsa3d(ig)) )
            if( ltrace .AND. tectyp .NE. RTRAC ) then
               do 10 i=1,MXSA3D
                  ptconc(i) = MAX( ptconc(i),BNDLPT )
   10          continue
            endif
        endif
c
c========================= Source Apportion End ========================
c
      enddo
c
c-----Modify BC's if some children are attached to each other
c
      do ic1 = 1,nchdrn(ip)
        ig1 = idchdrn(ic1,ip)
        do ic2 = 1,nchdrn(ip)
          ig2 = idchdrn(ic2,ip)
          call bcmodfy(nspec,ncol(ig1),nrow(ig1),nlay(ig1),ncol(ig2),
     &                 nrow(ig2),nlay(ig1),i1(ig1),j1(ig1),i2(ig1),
     &                 j2(ig1),i1(ig2),j1(ig2),i2(ig2),j2(ig2),
     &                 nmesh(ig1),nmesh(ig2),conc(iptr4d(ig1)),
     &                 conc(iptr4d(ig2)) )
c
c======================== Source Apportion Begin =======================
c
           if( ltrace .OR. lddm ) then
               call bcmodfy(ntotsp,ncol(ig1),nrow(ig1),nlay(ig1),
     &                ncol(ig2),nrow(ig2),nlay(ig1),i1(ig1),j1(ig1),
     &                i2(ig1),j2(ig1),i1(ig2),j1(ig2),i2(ig2),j2(ig2),
     &                nmesh(ig1),nmesh(ig2),ptconc(ipsa3d(ig1)),
     &                ptconc(ipsa3d(ig2)) )
               if( ltrace .AND. tectyp .NE. RTRAC ) then
                  do 20 i=1,MXSA3D
                     ptconc(i) = MAX( ptconc(i),BNDLPT )
   20             continue
               endif
           endif
c
c========================= Source Apportion End ========================
c
        enddo
      enddo
c
      return
      end
