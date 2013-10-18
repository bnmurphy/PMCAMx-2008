       subroutine wrtipr(iendat,endtim)
c
c-----CAMx v4.02 030709
c
c     This routine writes to the output file for the Integrated Process
c     Rate (IPR) data for the Process Analysis algorithm.  Each record 
c     contains the all of the data for once cell and species for the 
c     specified time period.
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
c
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "chmstry.com"
      include "filunit.com"
      include "procan.com"
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
c
c  --- loop over the species that should br written to output file
c      (as defines by regular model average species list)
c
      do l=1,navspc
        ispc = lavmap(l)
c
c  --- loop over the number of sub-domain cells ---
c
        do icel=1,npa_cels
c
c  --- calculate the average volume of cell,
c      and re-initialize the counter ---
c
           if( npastep(icel,ispc) .GT. 0 ) then
              cipr(IPR_VOL,icel,ispc) = cipr(IPR_VOL,icel,ispc) / 
     &                                          npastep(icel,ispc)
           else
              cipr(IPR_VOL,icel,ispc) = 0.0
           endif
c
c  --- write out the record for this cell ---
c
           write(ipr_unit,ERR=7000) iendat, endtim, 
     &                    spname(ispc), ipadom(icel), ipanst(icel),
     &                           ipax(icel), ipay(icel), ipaz(icel),
     &                                 (cipr(i,icel,ispc),i=1,NPAPRC)
c
c  --- next subdomain cell ---
c
        enddo
c
c  --- next species
c
      enddo
c
c  --- re-initialize counters to zero ---
c
      do icel=1,npa_cels
         do ispc=1,MXSPEC
           npastep(icel,ispc) = 0
         enddo
      enddo
      goto 9999
c
c----Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in WRTIPR:'
      write(iout,'(1X,2A)',ERR=9999) 'Writing data to the output ',
     &                         'Integrated Process Rate (.ipr) file.'
      write(iout,'(10X,A,I8.5,5X,A,F8.1)') 
     &      'Date: ',iendat,'Time: ',endtim
      call camxerr()
c
 9999 continue
      return
      end
