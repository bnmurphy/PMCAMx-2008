      subroutine fillpt(igroup,nptsrc,mxtemf,mxptsrc,lpig,lpigsa,
     &                  lvocsp,lnoxsp,crbnum,emsvoc,emsnox,emspts)
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.02 030709
c 
c     FILLPT adds point emissions to nox and vos emission catagories
c
c     Input arguments:
c        igroup           emission group id
c        mxtemf           maximum emission group number
c        mxptsrc          maximum point source number
c        lpig             flag for pig module
c        lpigsa           flag for pig modele in source apportionment
c        lvocsp           flag for voc species
c        lnoxsp           flag for nox species
c        crbnum           carbon number if the species
c        emsgrd           area emision for the species
c
c     Output arguments:
c        emsvoc           voc emission
c        emsnox           nox emission
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Variable declarartion:
c-----------------------------------------------------------------------
c
      logical   lpig, lpigsa(nptsrc),lvocsp,lnoxsp
      real      emsvoc(0:mxtemf,mxptsrc), emsnox(0:mxtemf,mxptsrc)
      real      emspts(mxptsrc)
c
c   --- add to toals for VOC or NOx ---
c
      do 50 i=1,nptsrc
c
c  --- skip the adding of emissions if it is used by PiG ---
c
         if( lpig .AND. lpigsa(i) .AND. lnoxsp ) goto 50
         if( lvocsp ) then
             emsvoc(igroup,i) = emsvoc(igroup,i) + 
     &                                 emspts(i) * crbnum
         else if( lnoxsp ) then
             emsnox(igroup,i) = emsnox(igroup,i) + emspts(i)
         endif
  50  continue
c
      return
      end
