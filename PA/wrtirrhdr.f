       subroutine wrtirrhdr(ibgdat,btim,iendat,etim)
c
c-----CAMx v4.03 031205
c
c     This routine writes the header to the output file for the Integrated
c     Reaction Rates (IRR) data for the Process Analysis algorithm.  The 
c     header data is designed to provide information about what data is 
c     contained in the rest of the file.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        11/06/01  Input dates are now Julian
c        12/03/03  Fixed bug in writing orgin of coarse grid
c
c     Input arguments:
c        btim       starting time for this period
c        ibgdat     starting date for this period (YYJJJ)
c        etim       ending time for this period
c        iendat     ending date for this period (YYJJJ)
c
c     Subroutines Called:
c        none
c
c     Called by:
c        STARTUP
c
      include "camx.prm"
      include "camx.com"
      include "chmstry.com"
      include "filunit.com"
      include "grid.com"
      include "procan.com"
c
c-----Argument declarations
c
      integer ibgdat,iendat
      real    btim,etim 
c
c-----Local variables
c
c-----Entry point
c
c --- write the date and time of the simulation ---
c
      write(irr_unit,ERR=7000) runmsg
      write(irr_unit,ERR=7000) ibgdat, btim, iendat, etim 
c
c --- write the information about the grids ---
c
c
c --- write the information about the grids ---
c
      write(irr_unit,ERR=7000) ngrid
      write(irr_unit,ERR=7000) xorg, yorg, ncol(1), nrow(1), 
     &                                              delx, dely, iuzon
      do igrd=2,ngrid
        xsize = delx / FLOAT( meshold(igrd) )
        ysize = dely / FLOAT( meshold(igrd) )
        orgx = (inst1(igrd)-1)*delx + xorg + xsize
        orgy = (jnst1(igrd)-1)*dely + yorg + ysize
        write(irr_unit,ERR=7000) orgx, orgy, ncol(igrd), nrow(igrd),
     &                                          xsize, ysize, iuzon
      enddo
c
c --- write out the sub-domain information ---
c
      write(irr_unit,ERR=7000) npadom
      do idom=1,npadom
         write(irr_unit,ERR=7000) ipagrd(idom), i_sw(idom), i_ne(idom),
     &                  j_sw(idom), j_ne(idom), b_lay(idom), t_lay(idom)
      enddo
c
c ---- write the number of reactions ----
c
      write(irr_unit,ERR=7000) nirrrxn
c      
      goto 9999
c
c----Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in WRTIRRHDR:'
      write(iout,'(1X,2A)',ERR=9999) 'Writing header to the output ',
     &                         'Integrated Reaction Rates (.irr) file.'
      call camxerr()
c
 9999 continue
      return
      end
