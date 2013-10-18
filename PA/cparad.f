      subroutine cparad( rad, nr, pa, npa, nn, dt )
c
c-----CAMx v4.02 030709
c
c     CPARAD saves radical concentrations for CPA output
c     and sets the CPA names and position numbers on the 
c     first call from PASETUP
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        none
c
c     Input arguments:
c        rad                 radical concentrations array (ppm)
c        nr                  dimension of radical array
c        pa                  CPA parameter array (ppb units)
c        npa                 dimension of pa array
c        nn                  counter of CPA parameters filled
c        dt                  time step (hours)
c
c     Output arguments:
c        pa                  CPA parameter array (ppb units)
c        nn                  counter of CPA parameters filled
c
c     Routines called:
c        none
c
c     Called by:
c        CHEMDRIV
c        PASETUP
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'tracer.com'
      include 'procan.com'
      include 'chmstry.com'
c
      integer npa, nr, nn
      real    dt, dtfact, ppbfact
      real    pa(npa), rad(nr)
c
c-----Entry point
c
      ppbfact = 1000.
c
c --- first calculate factor for time averaging concentrations ---
c
      dtfact = dt/dtout/60.
c
      nn = nn + 1
      ptname(nn)  = 'OH'
      pa(nn) = dtfact*rad(kOH)*ppbfact
c
      nn = nn + 1
      ptname(nn)  = 'HO2'
      pa(nn) = dtfact*rad(kHO2)*ppbfact
c
      nn = nn + 1
      ptname(nn)  = 'NO3'
      pa(nn) = dtfact*rad(kNO3)*ppbfact
c
      nn = nn + 1
      ptname(nn)  = 'N2O5'
      pa(nn) = dtfact*rad(kN2O5)*ppbfact
c
c
      if( nn .GT. MXCPA ) then
         write(iout,'(//,a)') 'ERROR in CPARAD:'
         write(iout,*) 'Number of outputs requested exceeds limit.'
         write(iout,*) 'Increase parameter MXCPA to', nn
         call camxerr()
      end if
c
      return
      end
