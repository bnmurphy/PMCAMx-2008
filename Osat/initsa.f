C**** INITSA
c
      subroutine initsa(version,nxx,nxy,idate,begtim,jdate,endtim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine calls all of the routines that set up the run for the
c   passive tracer algorithm.  The first thing done on routine OPENPT
c   is to set the flag for using passive tracer alogithm.  If this
c   flag is false coming out of the subroutine then all other routines
c   are skipped.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Argument description:
c      Inputs:
c         version C  version number of the model
c         nxx     I  number of cells in the X-direction
c         nxy     I  number of cells in the Y-direction
c         jdate   I  beginning date of the simulation (YYJJJ)
c         begtim  R  beginning time of the simulation
c         idate   I  ending date of the simulation (YYJJJ)
c         endtim  R  ending time of the simulation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     05/27/96   --gwilson--    Original development
c     11/06/01   --cemery--     Input dates are now Julian
c     01/30/02   --gwilson--    Added code for RTRAC probing tool
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      character*20 version
      integer      nxx(*)
      integer      nxy(*)
      integer      idate
      integer      jdate
      real         begtim
      real         endtim
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  ---- set the global variables to the passed arguments ---
c
      do igrd=1,ngrid
        nxcell(igrd) = nxx(igrd)
        nycell(igrd) = nxy(igrd)
      enddo 
c
c  --- if doing RTRAC, call routine to read chemistry definitions file ---
c 
      if( tectyp .EQ. RTRAC ) then
         call rdchmrt(version)
         ntrtim = 0 
         nsaspc = ntotsp
         iptnox = 1
         iptvoc = ntotsp + 1
         ipto3n = 1
         ipto3v = 1
         ipttim = ntotsp + 1
c
c  ---- call routine to read the source mapping file ---
c
      else
          if( nregin .GT. 0 ) call resmap()
c
c   ---- call routine to set the species flags ----
c
          call stabsa()
c
c   ---- call routine to set up the tracer species lists ----
c
          if( ltrace ) then
            call specsa(idate,begtim,jdate,endtim)
          else if( lddm ) then
            call specddm( )
          endif
      endif
c
c  ---- return to the calling routine ---
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
