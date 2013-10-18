       subroutine wrtiprhdr(ibgdat,btime,iendat,etime)
c
c-----CAMx v4.02 030709
c
c     This routine writes the header to the output file for the Integrated
c     Process Rates (IPR) data for the Process Analysis algorithm.  The header 
c     data is designed to provide information about what data is contained
c     in the rest of the file.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        11/06/01  Input dates are now Julian
c
c     Input arguments:
c        btime      ending time for this period
c        ibgdat     ending date for this period (YYJJJ)
c        etime      ending time for this period
c        iendat     ending date for this period (YYJJJ)
c
c     Subroutines Called:
c
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
      integer ibgdat
      real    btime 
      integer iendat
      real    etime 
c
c-----Local variables
c
      character*25 prcnam(NPAPRC)
c
c-----Entry point
c
c  --- load process names into local array ---
c
       prcnam(IPR_INIT)    = 'Initial concentration    '
       prcnam(IPR_CHEM)    = 'Chemistry                '
       prcnam(IPR_AEMIS)   = 'Area emissions           '
       prcnam(IPR_PTEMIS)  = 'Point source emissions   '
       prcnam(IPR_PIGEMIS) = 'Plume-in-Grid change     '
       prcnam(IPR_WADV)    = 'West boundary advection  '
       prcnam(IPR_EADV)    = 'East boundary advection  '
       prcnam(IPR_SADV)    = 'South boundary advection '
       prcnam(IPR_NADV)    = 'North boundary advection '
       prcnam(IPR_BADV)    = 'Bottom boundary advection'
       prcnam(IPR_TADV)    = 'Top boundary advection   '
       prcnam(IPR_DADV)    = 'Dilution in the vertical '
       prcnam(IPR_WDIF)    = 'West boundary diffusion  '
       prcnam(IPR_EDIF)    = 'East boundary diffusion  '
       prcnam(IPR_SDIF)    = 'South boundary diffusion '
       prcnam(IPR_NDIF)    = 'North boundary diffusion '
       prcnam(IPR_BDIF)    = 'Bottom boundary diffusion'
       prcnam(IPR_TDIF)    = 'Top boundary diffusion   '
       prcnam(IPR_DDEP)    = 'Dry deposition           '
       prcnam(IPR_WDEP)    = 'Wet deposition           '
       prcnam(IPR_FAERO)   = 'Aerosol chemistry        '
       prcnam(IPR_FINAL)   = 'Final concentration      '
       prcnam(IPR_CONV)    = 'Units conversion         '
       prcnam(IPR_VOL)     = 'Average cell volume      '
c
c --- write the date and time of the simulation ---
c
      write(ipr_unit,ERR=7000) runmsg
      write(ipr_unit,ERR=7000) ibgdat, btime, iendat, etime 
c
c --- write the information about the grids ---
c
      write(ipr_unit,ERR=7000) ngrid
      write(ipr_unit,ERR=7000) xorg, yorg, ncol(1), nrow(1), delx, dely
      do igrd=2,ngrid
        xsize = delx / FLOAT( meshold(igrd) )
        ysize = dely / FLOAT( meshold(igrd) )
        orgx = (inst1(igrd)-1)*delx + xorg + xsize
        orgy = (jnst1(igrd)-1)*dely + yorg + ysize
        write(ipr_unit,ERR=7000) orgx, orgy, ncol(igrd), nrow(igrd), 
     &                                                    xsize, ysize
      enddo
c
c --- write the species output to file ---
c
      write(ipr_unit,ERR=7000) navspc
      do ispc=1,navspc
         write(ipr_unit,ERR=7000) spname(lavmap(ispc))
      enddo
c
c --- write out the sub-domain information ---
c
      write(ipr_unit,ERR=7000) npadom
      do idom=1,npadom
         write(ipr_unit,ERR=7000) ipagrd(idom), i_sw(idom), i_ne(idom),
     &                  j_sw(idom), j_ne(idom), b_lay(idom), t_lay(idom)
      enddo
c
c  ---- write the descriptions of the processes ---
c
      write(ipr_unit,ERR=7000) NPAPRC
      do i=1,NPAPRC
         write(ipr_unit,ERR=7000) prcnam(i)
      enddo
c
      goto 9999
c
c----Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in WRTIPRHDR:'
      write(iout,'(1X,2A)',ERR=9999) 'Writing header to the output ',
     &                         'Integrated Process Rate (.ipr) file.'
      call camxerr()
c
 9999 continue
      return
      end
