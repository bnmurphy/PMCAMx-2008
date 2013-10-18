      subroutine aggdep(ncolf,nrowf,ncol,nrow,ndpspc,i1,j1,
     &                  i2,j2,nmesh,dxfin,dyfin,depfldf,depfld)
c
c-----CAMx v4.02 030709
c  
c     AGGDEP aggregates child (fine) grid deposition fields to its parent
c     (coarse) grid using an arithmitic average
c                            
c     Copyright 2003 
c     ENVIRON International Corporation
c            
c     Modifications:  
c        none
c
c     Input arguments:
c        ncolf               number of columns on fine grid
c        nrowf               number of rows on fine grid
c        ncol                number of columns on parent grid
c        nrow                number of rows on parent grid
c        ndpspc              number of species (3*navspc)
c        i1,j1               starting indices of fine grid on parent grid
c        i2,j2               ending indices of fine grid on parent grid
c        nmesh               fine grid meshing factor relative to parent
c        dxfin               cell size in x-direction on fine grid
c        dyfin               cell size in y-direction on fine grid
c        depfldf             deposition array on fine grid
c
c     Output arguments:
c        depfld              deposition array on parent grid
c
c     Routines called:
c        none
c
c     Called by:
c        AGGR00
c
      dimension dxfin(nrowf),depfld(ncol,nrow,ndpspc),
     &          depfldf(ncolf,nrowf,ndpspc)
c
c-----Entry point
c
      do 60 ispc = 1,ndpspc
        do 50 j = j1,j2
          do 40 i = i1,i2
            tarea = 0.
            dtmp = 0.
            do jtmp = 1,nmesh
              j0 = (j-j1)*nmesh+jtmp + 1
              do itmp = 1,nmesh
                i0 = (i-i1)*nmesh+itmp + 1
                area = dxfin(j0)*dyfin
                tarea = tarea + area
                dtmp = dtmp + depfldf(i0,j0,ispc)*area
              enddo
            enddo
            depfld(i,j,ispc) = dtmp/tarea
  40      continue
  50    continue
  60  continue
c
      return
      end
