c*** PEKRCP
c
      subroutine pekrcp(conc,nox,noy,noz,nspc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c   Description:
c     This routine looks in the concentraction array and
c     searches for the regular model peak cell.  A global variable
c     is set to this cell value. The tracer concentrations are also
c     accumulated in the entire grid.  This is necessary because
c     the hourly peak may change from one time step to another.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c   Argument descriptions:
c     Inputs:
c         conc   R  concentrations from the regular model
c         nox    I  number of X cells in grid
c         noy    I  number of Y cells in grid
c         noz    I  number of layers
c         nspc   I  number of regular model species
c         dtime  R  length of this time step
c
c-----------------------------------------------------------------------
c   LOG:
c-----------------------------------------------------------------------
c
c    12/11/2000  --gwilson--  Corrected bug in accessing species list.
c                             Now uses the map for average species. 
c
c-----------------------------------------------------------------------
c   Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'bndary.com'
      include 'chmstry.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declarations:
c-----------------------------------------------------------------------
c
      integer   nox
      integer   noy
      integer   noz
      integer   nspc
      real      conc(nox,noy,noz,nspc)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer   idxo3, icl, jcl, i
      real      o3max
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- find the index of the ozone species ---
c
      idxo3 = 0
      do 10 i=1,nspc
        if( lo3sp(lavmap(i)) ) idxo3 = i
   10 continue
c
c  --- Initialize to dummy values and just bail if not
c      doing ozone as average species ---
c
      ipekcl(1) = 0
      ipekcl(2) = 0
      o3max = -999999.
      if( idxo3 .LE. 0 ) goto 9999
c
c  --- loop over all of the cells, looking for the peak ---
c
      do 20 jcl=2,noy-1
         if( ibeg(jcl) .EQ. -999 ) goto 20
         do 30 icl=ibeg(jcl),iend(jcl)
            if( conc(icl,jcl,1,idxo3) .GT. o3max ) then
               o3max = conc(icl,jcl,1,idxo3)
               ipekcl(1) = icl
               ipekcl(2) = jcl
            endif
   30    continue
   20 continue
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
c
