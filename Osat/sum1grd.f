      subroutine sum1grd(igroup,igrd,idx,ibeg,iend,emssum,emsgrd,
     &                   emsbas,emsoth,emslft,emstot)
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.02 030709
c
c     SUM1GRD sums up the area emission of one species for a given group
c     in a given grid
c
c       07/19/02  --gwilson-- Added seperate source area map for each grids.
c
c     Input argument:
c        igroup            group ID
c        igrd              grid ID
c        idx               specie ID
c        ibeg              begining row of the domain
c        iend              ending row of the domain
c        emsgrd            the species emission in the grid
c
c     Output arguments:
c        emssum            emission summed over grid
c        emsbas            base emission
c        emsoth            "otherwise" emission
c        emslft            leftover emission
c        emstot            total emission
c        
      include "camx.prm"
      include "grid.com"
      include "tracer.com"
c
c
      dimension ibeg(MXROWA),iend(MXROWA)
      real   emssum(MXSPEC,MXTRSP)
      real   emsgrd(MXCOLA,MXROWA)
      real   emsbas(MXSPEC,MXTRSP), emsoth(MXSPEC,MXTRSP)
      real   emslft(MXCOLA,MXROWA), emstot(MXCOLA,MXROWA)
c
c   --- sum up the emissions excluding the boundary for this grid ---
c
      do 40 j=2,nrow(igrd)-1
        istrt = 2
        ifin = ncol(igrd)-1
        if(igrd .EQ. 1) then
          if(ibeg(j) .EQ. -999) goto 40
          istrt = ibeg(j)
          ifin = iend(j)
        endif
        do 50 i=istrt,ifin
c
c   --- skip cell if this grid has a child in this cell ---
c
           ijcl = i + (j-1)*ncol(igrd)
           if( idfin(iptr2d(igrd)-1+ijcl) .NE. 0 ) goto 50
c
c  --- get the region for this cell from mapping array,
c      the grid cell should be the coarse grid  ----
c
           imap = igrmap(igrd,i,j)
           if( imap .LT. 0 .OR. imap .GT. nregin ) goto 50
c
c  --- calculate the index into the tracer species for this gruoup/region ---
c
           if( ngroup .GT. 0 ) then
c
c   --- if group is base emissions, add to "leftover" group ----
c
              if( igroup .EQ. 0 ) then
                if( leftovr ) then
                   inox = iemnox - 1 + imap + ngroup*nregin
                   emsbas(idx,inox) = emsbas(idx,inox) + emsgrd(i,j)
                   ivoc = iemvoc - 1 + imap + ngroup*nregin
                   emsbas(idx,ivoc) = emsbas(idx,ivoc) + emsgrd(i,j)
                endif
                emstot(i,j) = emstot(i,j) + emsgrd(i,j)
c
c   --- otherwise, add to this group/region and subtract from "leftover" ---
c
              else
                inox = iemnox-1 + imap + (igroup-1)*nregin
                emssum(idx,inox) = emssum(idx,inox) + emsgrd(i,j)
                ivoc = iemvoc-1 + imap + (igroup-1)*nregin
                emssum(idx,ivoc) = emssum(idx,ivoc) + emsgrd(i,j)
                if( leftovr ) then
                   inox = iemnox - 1 + imap + ngroup*nregin
                   emsoth(idx,inox) = emsoth(idx,inox) + emsgrd(i,j)
                   ivoc = iemvoc - 1 + imap + ngroup*nregin
                   emsoth(idx,ivoc) = emsoth(idx,ivoc) + emsgrd(i,j)
                endif
                emslft(i,j) = emslft(i,j) + emsgrd(i,j)
              endif
c
c   --- only using regular model emissions ---
c
           else
              inox = iemnox - 1 + imap
              emssum(idx,inox) = emssum(idx,inox) + 
     &                                        emsgrd(i,j)
              ivoc = iemvoc - 1 + imap
              emssum(idx,ivoc) = emssum(idx,ivoc) + 
     &                                        emsgrd(i,j)
           endif
  50    continue
  40  continue
c
      return
      end
