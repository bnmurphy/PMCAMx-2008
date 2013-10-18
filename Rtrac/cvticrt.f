c**** CVTICRT.F
c
      subroutine cvticrt(igrd,nox,noy,noz,nspsa,saconc,tpgrd,prgrd)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine converts the intial conditions in the gridded 
c   array from PPM to ug/m.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c      Argument description:
c       Outputs:
c           saconc   R  tracer concentrations
c       Inputs:
c           igrd     I  grid number
c           nox      I  number of X cells in the grid
c           noy      I  number of Y cells in the grid
c           noz      I  number of layers in the grid
c           nspsa    I  number of tracer species
c           tpgrd    I  3-D temperature field
c           prgrd    I  3-D pressure field
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/21/02   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'bndary.com'
      include 'camx.com'
      include 'tracer.com'
      include 'rtracchm.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer*4 igrd
      integer   nox
      integer   noy
      integer   noz
      integer   nspsa
      real      saconc(nox,noy,noz,nspsa)
      real      tpgrd(nox,noy,noz)
      real      prgrd(nox,noy,noz)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer ispc, i, j, ilay
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- loop over the species ---
c
      do ispc=1,nspsa
c
c   --- loop over layers ---
c
         do ilay=1,noz
c
c  --- loop over columns, skip the boundary cells ---
c
            do 10 j=2,noy-1
               icel1 = 2
               icel2 = nox - 1
               if( igrd .EQ. 1 ) then
                  if( ibeg(j) .EQ. -999 ) goto 10
                  icel1 = ibeg(j)
                  icel2 = iend(j)
               endif
c
c  --- loop over rows in this column ----
c
               do i=icel1,icel2
c
c  --- get the conversion factor to umol/m3, only if not lower bound ---
c
                    if( saconc(i,j,ilay,ispc) .GT. rtlbnd(ispc) ) then
                        if( ispc .LE. nrtgas ) then
                            convfac = densfac*273./ tpgrd(i,j,ilay)*
     &                                       prgrd(i,j,ilay)/1013.
                        else
                            convfac = 1.0
                        endif
                        saconc(i,j,ilay,ispc) = 
     &                              saconc(i,j,ilay,ispc) * convfac
                    endif
c
c  --- next cell ---
c
              enddo
 10         continue
c
c  --- if this is a nest, zero out the concentrations on the
c      boundary ----
c
            if( igrd .GT. 1 ) then
                do j=1,noy
                   saconc(1,j,ilay,ispc) =  0.
                   saconc(nox,j,ilay,ispc) =  0.
                enddo
                do i=1,nox
                   saconc(i,1,ilay,ispc) =  0.
                   saconc(i,noy,ilay,ispc) =  0.
                enddo
            endif
c
c  --- next layer ---
c
         enddo
c
c  --- next species ---
c
      enddo
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
