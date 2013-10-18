      subroutine expndlay(ncolf,nrowf,nlayf,nlay,nmshv,fval)
c
c-----CAMx v4.02 030709
c
c     EXPNDLAY expands 3-D values from coarse grid layer number 'nlay'
c     to fine grid layer number 'nlayf'.  Assigns constant values to
c     the fine layers that make up a given coarse layer
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        none
c
c     Input arguments:
c        ncolf             number of columns in fine grid
c        nrowf             number of rows in fine grid
c        nlayf             number of layers in fine grid
c        nlay              number of layers in coarse grid
c        nmshv             vertical mesh number
c        fval              array of fine grid values
c
c     Output arguments:
c        fval              array of fine grid values
c
c     Subroutine called:
c        none
c
c     Called by:
c        INTRPDAT
c        INTRPCNC
c
      dimension fval(ncolf,nrowf,nlayf),nmshv(nlay)
c
c-----Entry point
c
      kg1 = 0
      do kp = 1,nlay
        kg1 = kg1 + nmshv(kp)
      enddo
c
      do kp = nlay,1,-1
        do kg = kg1,kg1-nmshv(kp)+1,-1
          do j=1,nrowf
            do i=1,ncolf
              fval(i,j,kg) = fval(i,j,kp)
            enddo
          enddo
        enddo
        kg1 = kg
      enddo
c
      return
      end
