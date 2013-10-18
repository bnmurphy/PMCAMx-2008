      subroutine interpv(ncol,nrow,nlay,ncolf,nrowf,nlayf,nmesh,nmshv,
     &                   i1,j1,htcrs,htfin,fval)
c
c-----CAMx v4.02 030709
c
c     INTERPV vertically interpolates coarse layer fields to a fine 
c     layer structure
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        none
c
c     Input arguments:
c        ncol              number of columns in coarse grid
c        nrow              number of rows in coarse grid
c        nlay              number of layers in coarse grid
c        ncolf             number of columns in fine grid
c        nrowf             number of rows in fine grid
c        nlayf             number of layers in fine grid
c        nmesh             meshing factor for fine grid
c        nmshv             vertical mesh number
c        i1,j1             starting indices for fine grid
c        htcrs             coarse grid layer interface height field (m)
c        htfin             fine grid layer interface height field (m)
c        fval              layer-centered field on coarse structure
c
c     Output arguments:
c        fval              layer-centered field on fine structure
c
c     Subroutine called:
c        none
c
c     Called by:
c        STARTUP
c        INTRPDAT
c
      include 'camx.prm'
      dimension htcrs(ncol,nrow,nlay),htfin(ncolf,nrowf,nlayf),
     &          fval(ncolf,nrowf,nlayf),nmshv(nlay)
      real dtdz(MXLAYA+1),dum(MXLAYA),zmid(MXLAYA)
c
c-----Entry Point
c
      do j = 1,nrowf
        do i = 1,ncolf
          ii = i1 + (i-2)/nmesh
          jj = j1 + (j-2)/nmesh
          do kp = 1,nlay
            zmid(1) = htcrs(ii,jj,kp)/2.
            dum(kp) = fval(i,j,kp)
            if (kp.gt.1) 
     &        zmid(kp) = (htcrs(ii,jj,kp) + htcrs(ii,jj,kp-1))/2.
          enddo
          dtdz(1) = (dum(2) - dum(1))/(zmid(2) - zmid(1))
          dtdz(nlay+1) = (dum(nlay) - dum(nlay-1))/
     &                   (zmid(nlay) - zmid(nlay-1))
          do kp = 2,nlay
            dtdz(kp) = (dum(kp) - dum(kp-1))/(zmid(kp) - zmid(kp-1))
          enddo
c
          kg = 1
          do 10 kp = 1,nlay
            if (nmshv(kp).eq.1) then
              fval(i,j,kg) = dum(kp)
              kg = kg + 1
              goto 10
            else
              do kfin = kg,kg+nmshv(kp)-1
                zmidf = htfin(i,j,kfin)/2.
                if (kfin.gt.1) 
     &            zmidf = (htfin(i,j,kfin) + htfin(i,j,kfin-1))/2.
                kgrad = kp
                if (zmidf.ge.zmid(kp)) kgrad = kp + 1
                fval(i,j,kfin) = dum(kp) + 
     &                           dtdz(kgrad)*(zmidf - zmid(kp))
              enddo
              kg = kg + nmshv(kp)
            endif
 10       continue
        enddo
      enddo
c
      return
      end
