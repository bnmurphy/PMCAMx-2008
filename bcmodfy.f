      subroutine bcmodfy(nspec,ncol1,nrow1,nlay1,ncol2,nrow2,nlay2,
     &                   i1g1,j1g1,i2g1,j2g1,i1g2,j1g2,i2g2,j2g2,
     &                   nmsh1,nmsh2,conc1,conc2)
c
c-----CAMx v4.02 030709
c
c     BCMODFY set the portion of boundary conditions where two fine 
c     grids are attached to each other.  CAMx requires adjacent nests
c     to have the same layer structure (nlay1 == nlay2) and this is
c     is enforced here.
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c        None
c
c     Input arguments:
c        nspec             number of species
c        ncol1             number of columns of grid 1
c        nrow1             number of rows of grid 1
c        nlay1             number of layers of grid 1
c        ncol2             number of columns of grid 2
c        nrow2             number of rows of grid 2
c        nlay2             number of layers of grid 2
c        i*g1,j*g1         starting/ending row/columns for grid 1
c        i*g2,j*g2         starting/ending row/columns for grid 2
c        nnmsh1            meshing factor for grid 1
c        nnmsh2            mehs number of grid 2
c        conc1             concentrations on grid 1 (umol/m3)
c        conc2             concentrations on grid 2 (umol/m3)
c
c     Output arguments:
c        conc1             concentrations on grid 1 (umol/m3)
c        conc2             concentrations on grid 2 (umol/m3)
c
c
c     Routines called:
c        SETBC1D
c
c     Called by:
c        SETBC
c
      include "camx.prm"
      include "filunit.com"
c
      dimension cnew1(MX1D),cold1(MX1D),cnew2(MX1D),cold2(MX1D)
      dimension conc1(ncol1,nrow1,nlay1,nspec),
     &          conc2(ncol2,nrow2,nlay2,nspec)
c
c-----Entry point
c
      do 50 ispc = 1,nspec
        nlay = nlay1
c
c-----Case 1: grid 2 is to the east of grid 1
c
        jbeg0 = max0(j1g1,j1g2)
        jend0 = min0(j2g1,j2g2)
        if (i2g1+1.eq.i1g2 .and. jend0.ge.jbeg0) then
c
          if (nlay1.ne.nlay2) then
            write(iout,'(//,a)') 'ERROR in BCMODFY:'
            write(iout,*) 'Illegal nest configuration'
            write(iout,*) 'nlay1 != nlay2: ',nlay1,nlay2
            call camxerr()
          endif
c
          do 10 k=1,nlay
c
c-----Load old concentrations
c
            jbeg1 = (jbeg0 - j1g1)*nmsh1
            jend1 = (jend0 - j1g1)*nmsh1
            n1 = 0
            do j=jbeg1,jend1
              n1 = n1 +1
              cold1(n1) = conc1(ncol2-1,j,k,ispc)
            enddo
c
            jbeg2 = (jbeg0 - j1g2)*nmsh2
            jend2 = (jend0 - j1g2)*nmsh2
            n2 = 0
            do j=jbeg2,jend2
              n2 = n2 + 1
              cold2(n2) = conc2(2,j+1,k,ispc)
            enddo
c
c-----Compute new concentrations
c
            call setbc1d(n1,n2,cnew1,cold2)
            call setbc1d(n2,n1,cnew2,cold1)
c
c-----Load new concentrations for the boundary cells
c
            n1 = 0
            do j=jbeg1,jend1
              n1 = n1 +1
              conc1(ncol2,j+1,k,ispc) = cold1(n1)
            enddo
c
            n2 = 0
            do j=jbeg2,jend2
              n2 = n2 + 1
              conc2(1,j+1,k,ispc) = cold2(n2)
            enddo
  10      continue
        endif
c
c-----Case 2: grid 2 is to the north of grid 1
c
        ibeg0 = max0(i1g1,i1g2)
        iend0 = min0(i2g1,i2g2)
        if (j2g1+1.eq.j1g2 .and. iend0.ge.ibeg0) then
c
          if (nlay1.ne.nlay2) then
            write(iout,'(//,a)') 'ERROR in BCMODFY:'
            write(iout,*) 'Illegal nest configuration'
            write(iout,*) 'nlay1 != nlay2: ',nlay1,nlay2
            call camxerr()
          endif
c
          do 20 k=1,nlay
c
c-----load old concentrations
c
            ibeg1 = (ibeg0 - i1g1)*nmsh1
            iend1 = (iend0 - i1g1)*nmsh1
            n1 = 0
            do i=ibeg1,iend1
              n1 = n1 + 1
              cold1(n1) = conc1(i,nrow1-1,k,ispc)
            enddo
c
            ibeg2 = (ibeg0 - i1g2)*nmsh2
            iend2 = (iend0 - i1g2)*nmsh2
            n2 = 0
            do i=ibeg2,iend2
              n2 = n2 + 1
              cold2(n2) = conc2(i+1,2,k,ispc)
            enddo
c
c-----compute new concentrations
c
            call setbc1d(n1,n2,cnew1,cold2)
            call setbc1d(n2,n1,cnew2,cold1)
c
c-----load new concentrations for the boundary cells
c
            n1 = 0
            do i=ibeg1,iend1
              n1 = n1 + 1
              conc1(i+1,nrow1,k,ispc) = cold1(n1)
            enddo
c
            n2 = 0
            do i=ibeg2,iend2
              n2 = n2 + 1
              conc2(i+1,1,k,ispc) = cold2(n1)
            enddo
c
  20      continue
        endif
c
  50  continue
c
      return
      end
