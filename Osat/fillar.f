      subroutine fillar(igrd,igroup,mxtemf,mxcol,mxrow,nox,noy,
     &                  ibeg,iend,lvocsp,lnoxsp,crbnum,emsgrd,emsvoc,
     &                  emsnox)
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.02 030709
c 
c     FILLAR adds area emissions to nox and vos emission chatagories
c
c
c     Input arguments:
c        igrd             grid ID
c        igroup           emission group id
c        mxtemf           maximum emission group number
c        mxcol            maximum column number
c        mxrow            maximum row number
c        nox              number of columns
c        noy              number of rows
c        ibeg             begining row number
c        iend             ending row number
c        lvocsp           flag for voc species
c        lnoxsp           flag for nox species
c        crbnum           carbon number if the species
c        emsgrd           area emision for the species
c
c     Output arguments:
c        emsvoc           voc emission
c        emsnox           nox emission
c
      logical   lvocsp,lnoxsp
      dimension ibeg(mxrow),iend(mxrow)
      dimension emsgrd(mxcol,mxrow)
      real      emsvoc(0:mxtemf,mxcol,mxrow)
      real      emsnox(0:mxtemf,mxcol,mxrow)
c
c   --- add to totals for VOC or NOx ---
c
      do 30 j=2,noy-1
        ibegcl = 2
        iendcl = nox-1
        if( igrd .EQ. 1 ) then
           if( ibeg(j) .EQ. -999 ) goto 30
           ibegcl = ibeg(j)
           iendcl = iend(j)
        endif
        do 40 i=ibegcl,iendcl
          if( lvocsp ) then
             emsvoc(igroup,i,j) = emsvoc(igroup,i,j) + 
     &                                 emsgrd(i,j) * crbnum
          else if( lnoxsp ) then
             emsnox(igroup,i,j) = emsnox(igroup,i,j)+ emsgrd(i,j)
          endif
  40    continue
  30  continue
c
      return
      end
