c-----CAMx v4.02 030709
c
c     CAMx.COM contains general model variables
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        11/05/01  Added Y2K flag
c        4/10/03   Moved dtmax from control file to parameter in timestep.f
c
c-----------------------------------------------------------------------
c     Variables for simulation date/time and update intervals:
c
c     date    --  Julian date (YYJJJ) of the current time step
c     time    --  time (HHMM) of the current step
c     ly2k    --  Year 2000 flag (T is >=2000)
c     begdate --  Julian date (YYJJJ) of the beginning of the simulation
c     begtim  --  time (HHMM) of the beginning of the simulation
c     dtinp   --  time interval for updating meteorological inputs (min)
c     dtems   --  time interval for updating emission inputs (min)
c     dtout   --  time interval to output concentrations (min)
c     dtmax   --  maximum allowable time step (min)
c     deltat  --  time step size for each grid (seconds)
c     xyordr  --  binary (0 or 1) variable to define horizontal advection
c                 order (0=X-Y, 1=Y-X) 
c     MAXDT   --  parameter for maximum allowable time step (in minutes)
c-----------------------------------------------------------------------
c
      integer   date
      real      time
      logical   ly2k
      integer   begdate
      real      begtim
      real      dtinp
      real      dtout
      real      dtmax
      real      deltat(MXGRID)
      integer   xyordr
      real      MAXDT
c      
      common /timstf/ date, time, ly2k, begdate, begtim, 
     &                dtinp, dtems, dtout, dtmax, deltat, xyordr
      parameter( MAXDT = 15. )
c
c-----------------------------------------------------------------------
c     Variables for simulation identification:
c
c     runmsg  -- string conatining the simulation run ID
c-----------------------------------------------------------------------
c
      character*80 runmsg
c
      common /runid/ runmsg
c
c-----------------------------------------------------------------------
c     Variables for conversion constants:
c
c     densfac -- conversion factor to convert from PPM to umol/m3 (mol/m3)
c-----------------------------------------------------------------------
c
      real   densfac
c
      common /convrt/ densfac
c
c-----------------------------------------------------------------------
c     Variables for keeping track of CPU time:
c
c     tarray  -- array for cumulative CPU of entire simulation (seconds)
c     tarray2 -- array for elapsed CPU of the current process (seconds)
c-----------------------------------------------------------------------
c
      real   tarray(2) 
      real   tarray2(2) 
c
      common /runcpu/ tarray, tarray2
