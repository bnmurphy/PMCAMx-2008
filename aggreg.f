      subroutine aggreg(ncolf,nrowf,nlayf,ncol,nrow,nlay,nspec,i1,j1,
     &                  i2,j2,nmesh,nmshv,dxfin,dyfin,depfin,
     &                  concf,conc)
c
c-----CAMx v4.02 030709
c  
c     AGGREG aggregates child (fine) grid concentrations to its parent (coarse)
c     grid using an arithmitic average
c                            
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c            
c     Modifications:  
c        1/29/99   Fixed a bug in layer indexing of fine layer depth
c                  (kp vs. kg)  
c
c     Input arguments:
c        ncolf               number of columns on fine grid
c        nrowf               number of rows on fine grid
c        nlayf               number of layers on fine grid
c        ncol                number of columns on parent grid
c        nrow                number of rows on parent grid
c        nlay                number of layers on parent grid
c        nspec               number of species
c        i1,j1               starting indices of fine grid on parent grid
c        i2,j2               ending indices of fine grid on parent grid
c        nmesh               fine grid meshing factor relative to parent
c        nmshv               fine grid vertical mesh factor relative to parent
c        dxfin               cell size in x-direction on fine grid
c        dyfin               cell size in y-direction on fine grid
c        depfin              layer depths on fine grid
c        concf               concentration array on fine grid
c
c     Output arguments:
c        conc                concentration array on parent grid
c
c     Routines called:
c        none
c
c     Called by:
c        AGGR00
c
      dimension dxfin(nrowf),depfin(ncolf,nrowf,nlayf),
     &          conc(ncol,nrow,nlay,nspec),
     &          concf(ncolf,nrowf,nlayf,nspec),nmshv(nlay)
c
c-----Entry point
c
      do 60 ispc = 1,nspec
        do 50 j = j1,j2
          do 40 i = i1,i2
            kg1 = 1
            do 30 kp = 1,nlay
              tvol = 0.
              ctmp = 0.
              do  kg = kg1,kg1+nmshv(kp)-1
                do jtmp = 1,nmesh
                  j0 = (j-j1)*nmesh+jtmp + 1
                  do itmp = 1,nmesh
                    i0 = (i-i1)*nmesh+itmp + 1
                    vol = dxfin(j0)*dyfin*depfin(i0,j0,kg)
                    tvol = tvol + vol
                    ctmp = ctmp + concf(i0,j0,kg,ispc)*vol
                  enddo
                enddo
              enddo
              conc(i,j,kp,ispc) = ctmp/tvol
              kg1 = kg1+nmshv(kp)
  30        continue
  40      continue
  50    continue
  60  continue
c
      return
      end
