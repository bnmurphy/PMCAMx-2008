      subroutine aggr00(igrd,ip,icode)
c
c-----CAMx v4.02 030709
c
c     AGGR00 passes arrays from the common blocks to AGGREG
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        1/15/03     Added aggregation of deposition fields
c
c     Input arguments:
c        igrd                grid index
c        ip                  parent grid index
c        icode               flag to select when to sum mass on parent grid
c                              0 = sum both before and after aggregation
c                              1 = sum before aggregation
c                              2 = sum after aggregation
c                             >2 = do not sum
c
c     Output arguments:
c        none
c
c     Subroutine called:
c        MASSUM
c        AGGREG
c
c     Called by:
c        NESTING
c
      include "camx.prm"
      include "camx.com"
      include "camxfld.com"
      include "grid.com"
      include "chmstry.com"
      include "filunit.com"
c
c-----Entry point
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
c
c========================= Source Apportion End ========================
c
c
c========================= Process Analysis Begin ==============================
c
      include 'procan.com'
c
c========================= Process Analysis End ==============================
c
      write(*,'(a20,$)') 'aggrg ......'
      write(iout,'(a20,$)') 'aggrg ......'
c
      if (icode.le.1) then
        call massum(ip,nspec,ncol(ip),nrow(ip),nlay(ip),deltax(1,ip),
     &             deltay(ip),depth(iptr3d(ip)),conc(iptr4d(ip)),
     &             xmstmp(1,ip))
      endif
c
      call aggreg(ncol(igrd),nrow(igrd),nlay(igrd),ncol(ip),nrow(ip),
     &            nlay(ip),nspec,i1(igrd),j1(igrd),i2(igrd),j2(igrd),
     &            nmesh(igrd),nmshv(1,igrd),deltax(1,igrd),deltay(igrd),
     &            depth(iptr3d(igrd)),conc(iptr4d(igrd)),
     &            conc(iptr4d(ip)) )
      call aggdep(ncol(igrd),nrow(igrd),ncol(ip),nrow(ip),navspc*3,
     &            i1(igrd),j1(igrd),i2(igrd),j2(igrd),nmesh(igrd),
     &            deltax(1,igrd),deltay(igrd),depfld(iptrdp(igrd)),
     &            depfld(iptrdp(ip)) )
c
c======================== Source Apportion Begin =======================
c
      if( ltrace .OR. lddm ) then
         call aggreg(ncol(igrd),nrow(igrd),nlay(igrd),ncol(ip),nrow(ip),
     &            nlay(ip),ntotsp,i1(igrd),j1(igrd),i2(igrd),j2(igrd),
     &            nmesh(igrd),nmshv(1,igrd),deltax(1,igrd),deltay(igrd),
     &            depth(iptr3d(igrd)),ptconc(ipsa3d(igrd)),
     &            ptconc(ipsa3d(ip)) )
      endif
c
c========================= Source Apportion End ========================
c
c
c========================= Process Analysis Begin ==============================
c
      if( lirr ) then
         call aggreg(ncol(igrd),nrow(igrd),nlay(igrd),ncol(ip),nrow(ip),
     &            nlay(ip),ntotsp,i1(igrd),j1(igrd),i2(igrd),j2(igrd),
     &            nmesh(igrd),nmshv(1,igrd),deltax(1,igrd),deltay(igrd),
     &            depth(iptr3d(igrd)),ptconc(ipsa3d(igrd)),
     &            ptconc(ipsa3d(ip)) )
      endif
c
c========================= Process Analysis End ==============================
c
      if (icode.eq.0 .or. icode.eq.2) then
        call massum(ip,nspec,ncol(ip),nrow(ip),nlay(ip),deltax(1,ip),
     &              deltay(ip),depth(iptr3d(ip)),conc(iptr4d(ip)),
     &              xmass(1,ip))
        do l = 1,nspec 
          xmsfin(l,ip) = xmsfin(l,ip) + xmass(l,ip) - xmstmp(l,ip) 
        enddo
      endif
c
      tcpu = dtime(tarray2)
      write(*,'(a,f10.3)') '   CPU = ', tarray2(1)
      write(iout,'(a,f10.3)') '   CPU = ', tarray2(1)
c
      return
      end
