      subroutine massum(igrid,nospec,ncol,nrow,nlay,dx,dy,depth,conc,
     &                  xmass,subxmass)
c
c-----CAMx v4.02 030709
c
c     MASSUM sums up mass on a given grid, not including boundary cells
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        igrid               grid index
c        nospec              number of species
c        ncol                number of columns
c        nrow                number of rows
c        nlay                number of layers
c        dx                  cell size in x-direction (m)
c        dy                  cell size in y-direction (m)
c        depth               cell depth (m)
c        conc                concentration field (umol/m3)
c
c     Output arguments:
c        xmass               grid mass (umol)
c
c     Routines called:
c        none
c
c     Called by:
c        CAMx
c        AGGR00
c        CHEMRXN
c
      include "camx.prm"
      include "bndary.com"
c
      real*8 xmass,dtmp
      real*8 subxmass
      integer subd(ncol,nrow)
      dimension conc(ncol,nrow,nlay,nospec),xmass(nospec*nlay), subxmass(nospec*nlay)
      dimension dx(nrow),depth(ncol,nrow,nlay)
c
c-----Entry point
c
      call subdomain(subd)

      do 50 is = 1,nospec
        do 30 k = 1,nlay
          xmass(k+(is-1)*nlay) = 0.     !<-BNM altered for multi-layer mass budget
          subxmass(k+(is-1)*nlay) = 0.  !<-BNM altered for subdomain mass budget (9-24-09)
          do 20 j = 2,nrow-1
            i1 = 2
            i2 = ncol-1
            if (igrid.eq.1) then
              if (ibeg(j).eq.-999) goto 20
              i1 = ibeg(j)
              i2 = iend(j)
            endif
            do i = i1,i2
              dtmp = conc(i,j,k,is)*dx(j)*dy*depth(i,j,k)
              xmass(k+(is-1)*nlay) = xmass(k+(is-1)*nlay) 
     &                  + dtmp            !<-BNM altered for multi-layer mass budget
              subxmass(k+(is-1)*nlay) = subxmass(k+(is-1)*nlay) 
     &                  + dtmp*subd(i,j)  !BNM 9-23-09
            enddo
  20      continue
  30    continue
c
  50  continue
c
      return
      end
