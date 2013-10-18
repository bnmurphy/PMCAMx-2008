c**** STABSA
c
      subroutine stabsa()
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine sets the flags that determine if each species in the
c   species list should be added to the VOC or NOx tracers.  The carbon
c   number for the VOC species is also put into the appropriate place
c   in the array.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/27/96   --gwilson--    Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'chmstry.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   NSACBIV    I    the number of VOC species names to be checked (CBIV)
c   NSASPRC    I    the number of VOC species names to be checked (SAPRC)
c   MXSANAM    I    maximum sixe of the arrays
c
      integer   NSACBIV
      integer   NSASPRC
      integer   MXSANAM
c
      parameter( NSACBIV = 10 )
      parameter( NSASPRC = 14 )
      parameter( MXSANAM = NSACBIV + NSASPRC )

c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 spcnam(MXSANAM), spccbiv(NSACBIV), spcsprc(NSASPRC)
      character*10 nonam, no2nam, o3name
      integer      ispc, nnames, i
      real         carbon(MXSANAM), react(MXSANAM)
      real         carbcbiv(NSACBIV), reactcbiv(NSACBIV)
      real         carbsprc(NSASPRC), reactsprc(NSASPRC)
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data spccbiv /'PAR       ','OLE       ','ETH       ',
     &              'TOL       ','XYL       ','FORM      ',
     &              'ALD2      ','ISOP      ','MEOH      ',
     &              'ETOH      '/
      data carbcbiv /1.,            2.,           2.,
     &               7.,            8.,           1.,
     &               2.,            5.,           1.,
     &               2./
      data reactcbiv  /1203.,        42000.,      11920.,
     &                 9150.,        36200.,      15000.,
     &                 24000.,       142000.,     1363.,
     &                 4791./
c
      data spcsprc /'ACET      ','ALK1      ','ALK2      ',
     &              'ARO1      ','ARO2      ','CCHO      ',
     &              'ETH       ','FORM      ','ISOP      ',
     &              'MEK       ','OLE1      ','OLE2      ',
     &              'OLE3      ','RCHO      '/
      data carbsprc /3.00,        4.55,         8.82,
     &               7.10,        8.82,         2.00,
     &               2.00,        1.00,         5.00,
     &               4.00,        3.34,         5.40,
     &              10.00,        3.00/
      data reactsprc  /   324.,   4690.,      15350.,
     &                   8487.,  41126.,      23289.,
     &                  12594.,  14428.,     147570.,
     &                   1708.,  40000.,      93737.,
     &                 134970.,  29267./
      data nonam  /'NO        '/
      data no2nam /'NO2       '/
      data o3name /'O3        '/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- load the local arrays based on the mechanism ---
c
      if( idmech .EQ. 5 ) then
         nnames = NSASPRC
         do i=1,nnames
            spcnam(i) = spcsprc(i)
            carbon(i) = carbsprc(i)
            react(i) = reactsprc(i)
         enddo
      else
         nnames = NSACBIV
         do i=1,nnames
            spcnam(i) = spccbiv(i)
            carbon(i) = carbcbiv(i)
            react(i) = reactcbiv(i)
         enddo
      endif
c
c   --- loop over the species in the chemparam file and 
c       intialize the flags ---
c
      do 10 ispc=1,nspec
         lvocsp(ispc) = .FALSE.
         lnoxsp(ispc) = .FALSE.
         lo3sp(ispc) = .FALSE.
c
c   --- check for O3 species ----
c
         if( spname(ispc) .EQ. o3name ) then
              lo3sp(ispc) = .TRUE.
              goto 10
         endif
c
c   --- check for NOx species ----
c
         if( spname(ispc) .EQ. nonam .OR. spname(ispc) .EQ. no2nam) then
              lnoxsp(ispc) = .TRUE.
              goto 10
         endif
c
c   ---- loop over species that are known to be VOC species ---
c
         do 20 i=1,nnames
            if( spname(ispc) .EQ. spcnam(i) ) then
               lvocsp(ispc) = .TRUE.
               crbnum(ispc) = carbon(i)
               reacrt(ispc) = react(i)
               goto 10
            endif
   20    continue
c
c  --- next species in the list ---
c
   10 continue
c      
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
