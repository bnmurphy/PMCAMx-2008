       subroutine wrtirr(iendat,endtim)
c
c-----CAMx v4.02 030709
c
c     This routine writes to the output file for the Integrated Reaction
c     Rates (IRR) data for the Process Analysis algorithm.  Each record 
c     contains the all of the data for once cell and species for the 
c     specified time period.
c     This routine also calls the subroutines that will write the 
c     gridded Chemical Process Analysis data.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        endtim     ending time for this period
c        iendat     ending date for this period
c
c     Subroutines Called:
c       WRTCGCPA
c       WRTFGCPA
c
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'chmstry.com'
      include 'filunit.com'
      include 'grid.com'
      include 'tracer.com'
      include 'procan.com'
c
c-----Argument declarations
c
      integer iendat
      real    endtim 
c
c-----Local variables
c
c
c-----Entry point
c
c  --- loop over the number of sub-domain cells ---
c
      do icel=1,npa_cels
c
c  --- write out the record for this cell ---
c
           write(irr_unit,ERR=7000) iendat, endtim, 
     &              ipadom(icel), ipanst(icel),
     &                  ipax(icel), ipay(icel), ipaz(icel),
     &                                  (cirr(icel,i),i=1,nirrrxn)
c
c  --- next subdomain cell ---
c
      enddo
c
c  --- call routine to write the gridded CPA array for the coarse grid ---
c
      if(lsfcfl(IDXCRS) ) call wrtcgcpa(iendat,endtim,ncol(1),
     &                              nrow(1),nlay(1),ntotsp,ptconc(1))
c
c  --- call routine to write the gridded CPA array for the fine grids ---
c
      if(lsfcfl(IDXFIN) .AND. ngrid .GT. 1) call wrtfgcpa(iendat,endtim)
c
      goto 9999
c
c----Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in WRTIRR:'
      write(iout,'(1X,2A)',ERR=9999) 'Writing data to the output ',
     &                         'Integrated Reaction Rates (.irr) file.'
      write(iout,'(10X,A,I8.5,5X,A,F8.1)') 
     &      'Date: ',iendat,'Time: ',endtim
      call camxerr()
c
 9999 continue
      return
      end
