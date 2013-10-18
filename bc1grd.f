      subroutine bc1grd(nspec,ncol,nrow,nlay,ncolf,nrowf,nlayf,i1,j1,
     &                  nmesh,nmshv,conc,concf)
c
c-----CAMx v4.02 030709
c
c     BC1GRD sets up boundary conditions for one fine grid using 
c     concentrations of its parent grid
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        none
c
c     Input arguments:
c        nspec             Number of species
c        ncol              Number of columns on the parent grid
c        nrow              Number of rows on the parent grid
c        nlay              Number of layers on the parent grid
c        ncolf             Number of columns on the fine grid
c        nrowf             Number of rows on the fine grid
c        nlayf             Number of layers on the fine grid
c        i1                Starting i-index of fine grid
c        j1                Starting j-index of fine grid
c        nmesh             meshing factor relative to parent grid
c        nmshv             vert factor relative to parent grid
c        conc              concentration on the coarse grid (umol/m3)
c
c     Output arguments:
c        concf             concentration on the fine grid (umol/m3)
c
c     Routines called:
c        none
c
c     Called by:
c        SETBC 
c
      dimension conc(ncol,nrow,nlay,nspec),
     &          concf(ncolf,nrowf,nlayf,nspec),nmshv(nlay)
c
c-----Entry point
c
      do 30 ispc = 1,nspec
        kg1 = 1
        do 20 kp=1,nlay
          do 10 kg=kg1,kg1+nmshv(kp)-1
c
c-----Southern boundary
c
            j=1
            jj=j1 - 1
            do i=2,ncolf-1
              ii=i1+(i-1)/nmesh
              concf(i,j,kg,ispc) = conc(ii,jj,kp,ispc)
            enddo
c
c-----Northern boundary
c
            j=nrowf
            jj= j1 + (nrowf-1)/nmesh
            do i=2,ncolf-1
              ii=i1+(i-1)/nmesh
              concf(i,j,kg,ispc) = conc(ii,jj,kp,ispc)
            enddo
c
c-----Western boundary
c
            i=1
            ii=i1 - 1
            do j=2,nrowf-1
              jj=j1+(j-1)/nmesh
              concf(i,j,kg,ispc) = conc(ii,jj,kp,ispc)
            enddo
c
c-----Eastern boundary
c
            i=ncolf
            ii=i1 + (ncolf - 2)/nmesh
            do j=2,nrowf-1
              jj=j1+(j-1)/nmesh
              concf(i,j,kg,ispc) = conc(ii,jj,kp,ispc)
            enddo
  10      continue
          kg1 = kg
  20    continue
c
  30  continue
c
      return
      end
