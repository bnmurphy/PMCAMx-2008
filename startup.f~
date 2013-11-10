      subroutine startup(version,inptim,inpdate,emstim,emsdate,haztim,
     &                   hazdate,ozntim,ozndate,bndtim,bnddate,wrttim,
     &                   wrtdate,endtim,enddate)
c
c-----CAMx v4.03 031205
c
c     STARTUP is the main initialization and setup routine for CAMx.
c     It performs the following tasks:
c        - initializes certain vector/array and scalar variables 
c        - reads and checks the user input file
c        - sets model simulation clock
c        - opens all I/O files
c        - reads all time-invariant files and look-up tables
c        - reads/writes headers from/to UAM-formatted I/O files
c        - calculates grid parameters
c        - initializes the PiG submodel
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        1/29/99   Added diagnostic printout of layer nesting, and error
c                  messages if layer nesting calculation is invalid
c        2/12/99   Removed assignment of negative longitude for xorg
c        4/26/99   Added Piecewise Parabolic Method for horizontal advection
c        5/24/99   Fixed bug in filling idfin array
c        10/20/00  Added CAMx version as first record on control file
c        11/06/01  Added Y2K flag; conversion of simulation date to Julian 
c                   now done immediately after reading from control file
c        1/21/02   Added RTRAC tech type to OSAT
c        1/25/02   Revised input I/O frequencies and max time step to minutes,
c                  and improved checks on values
c        7/5/02    Added code to handle new IRON-PiG option
c        1/10/03   Added prep of deposition output files 
c        03/21/03  Removed the OSAT technology type OPPAT
c
c     Input arguments:
c        version             model version character string
c
c     Output arguments:
c        inptim              next time to read environmental fields (HHMM)
c        inpdate             next date to read environmental fields (YYJJJ)
c        emstim              next time to read emissions (HHMM)
c        emsdate             next date to read emissions (YYJJJ)
c        haztim              next update time for haze map (HHMM)  
c        hazdate             next update date for haze map (YYJJJ)  
c        ozntim              next update time for ozone map (HHMM)  
c        ozndate             next update date for ozone map (YYJJJ) 
c        bndtim              next update time for boundary conditions (HHMM)  
c        bnddate             next update date for boundary conditions (YYJJJ) 
c        wrttim              next time to output concentrations (HHMM)
c        wrtdate             next date to output concentrations (YYJJJ)
c        endtim              model end time (HHMM)
c        enddate             model end date (YYJJJ)
c
c     Routines Called:
c        CALDATE,  OPENFILS, CHMPREP, METPREP,  BNDPREP,  GRDPREP, NSTPREP,
c        IASSGN2D, INTERP2D, READPI0, READZPWT, VNMSHCAL, SRFPREP, LUASSGN,
c        PNTPREP,  AREAPREP, CNCPREP, PIGPREP,  JULDATE,  INTERPV, FINWIND
c        RASSGN3D, DEPPREP
c
c     Called by:
c        CAMx
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'chmstry.com'
      include 'ahomap.com'
      include 'grid.com'
      include 'pigsty.com'
      include 'flags.com'
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
c
c========================= Source Apportion End ========================
c
c
c======================== Process Analysis Begin =======================
c
      include 'procan.com'
c
c========================= Process Analysis End ========================
c
c     integer begdate,inpdate,emsdate,wrtdate,enddate,hazdate,ozndate,
      integer inpdate,emsdate,wrtdate,enddate,hazdate,ozndate,
     &        bnddate
      real inptim,emstim,wrttim,endtim,haztim,ozntim,bndtim
      character*200 inprec,filroot,ctlfil
      character*80  action
      character*20 version
      character*10 grdflg,advflg,tecsav,slvstr,camxv,camxvin,pigstr
      integer   nemiss, nlines
      integer   cbdate,cedate
      logical   lexist,lbot,lfound, lpig
c
      data inp /3/
      data ctlfil /'CAMx.in'/
      data camxv /'VERSION4.0'/

      integer domlen_vec2d      ! Pavan
c
c-----Entry point
c
      llatlon = .FALSE.
      lutm = .FALSE.
      lpolar = .FALSE.
      lambrt = .FALSE.
      lbot = .FALSE.
      irec = 0
c
c-----Open user control file
c
      filroot = ' '
      inquire(file=ctlfil,exist=lexist)
      if( .NOT. lexist ) then
          write(*,'(//,a)') 'ERROR in STARTUP:'
          write(*,'(2A)') 'Cannot find control file: ',
     &                     ctlfil(:istrln(ctlfil))
          stop
      endif
      filroot = ctlfil
      open(unit=inp,file=ctlfil,STATUS='UNKNOWN',ERR=7005)
c
c-----Check CAMx version number on control file
c
      irec = irec + 1
      read(inp,'(A)',ERR=7006,END=7004) inprec
      camxvin = inprec(21:30)
      call jstlft(camxvin)
      call toupper(camxvin)
      if (camxvin.ne.camxv) then
        write(*,'(//,a)') ' ERROR in STARTUP:'
        write(*,'(//,a)') ' CAMx version in control file INVALID'
        write(*,'(a,a)')  ' Expecting: ',camxv
        write(*,'(a,a)')  '     Found: ',camxvin
        write(*,'(1X,3A)') 'Make sure you are using a control ',
     &                               'file designed for ',version
        stop
      endif
c
c-----Read run message and I/O file root name
c
      irec = irec + 1
      action = 'Reading the run message.'
      read(inp,'(A)',ERR=7006,END=7004) inprec
      runmsg = inprec(21:80)
      irec = irec + 1
      action = 'Reading the output filename root.'
      read(inp,'(A)',ERR=7006,END=7004) inprec
      filroot = inprec(21:)
      call jstlft( filroot )
      ii = istrln( filroot )
c
c-----Open ASCII output, diagnostic, and mass reporting files
c
      filroot(ii+1:ii+4) = '.out'
      iout = 7
      open(unit=iout,file=filroot(1:ii+4),status='UNKNOWN',ERR=7005)
c
      filroot(ii+1:ii+5) = '.diag'
      idiag = 8
      open(unit=idiag,file=filroot(1:ii+5),status='UNKNOWN',ERR=7000)
c
      filroot(ii+1:ii+5) = '.mass'
      imass = 9 
      open(unit=imass,file=filroot(1:ii+5),status='UNKNOWN',ERR=7000)
      iifroot = ii
c
c-----Write model version to output and diagnostic files
c
      write(iout,8000) version(:istrln(version))
      write(idiag,8000) version(:istrln(version))
 8000 format(//,30x,20('*'),/,30x,a,/,30x,20('*'),//) 
      write(iout,8001) runmsg(:istrln(runmsg))
      write(idiag,8001) runmsg(:istrln(runmsg))
 8001 format(/,a,/)
c
c-----Read run control
c
      irec = irec + 1
      action = 'Reading the start date/time.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),*,ERR=7002,END=7002) ibyr,ibmo,ibdy,begtim
      call cvtdate(ibyr,ibmo,ibdy,begtim,cbdate,begdate,ly2k)
c
      irec = irec + 1
      action = 'Reading the end date/time.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),*,ERR=7002,END=7002) ieyr,iemo,iedy,endtim
      call cvtdate(ieyr,iemo,iedy,endtim,cedate,enddate,ly2k)
c
      irec = irec + 1
      action = 'Reading the max timestep parameters.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),*,ERR=7002,END=7002) dtmax,dtinp,dtems,dtout
c
      irec = irec + 1
      action = 'Reading the number of cells for coarse grid.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),*,ERR=7002,END=7002) ncol(1),nrow(1),nlay(1)
      inst1(1) = 2
      inst2(1) = ncol(1)-1
      jnst1(1) = 2
      jnst2(1) = nrow(1)-1
c
      irec = irec + 1
      action = 'Reading the grid projection flag.'
      read(inp,'(A)',ERR=7001,END=7004) inprec
      grdflg = inprec(21:30)
      call jstlft( grdflg )
      call toupper( grdflg )
      if( grdflg .EQ. 'LATLON    ' ) then
         llatlon = .TRUE.
      elseif( grdflg .EQ. 'UTM       ' ) then
         lutm = .TRUE.
      elseif( grdflg .EQ. 'POLAR     ' ) then
        lpolar = .TRUE.
      elseif( grdflg .EQ. 'LAMBERT   ' ) then
        lambrt = .TRUE.
      else
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,'(3A)') 'Incorrect coordinate ID read from ',
     &                     'control file: ',inprec(21:30)
        write(iout,'(A)') 'CAMx supports: LATLON, UTM, POLAR, LAMBERT'
        call camxerr()
      endif
c
      irec = irec + 1
      action = 'Reading the coarse grid origin definition.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      if( llatlon ) then
         read(inprec(21:),*,ERR=7002,END=7002) xorg,yorg,delx,dely
      elseif( lutm ) then
         read(inprec(21:),*,ERR=7002,END=7002) xorg,yorg,delx,dely,iuzon
         if( iuzon .eq. 0 ) then
           write(iout,'(//,A)')'ERROR in STARTUP:'
           write(iout,'(A)')   '  The UTM zone can not be set to zero'
           write(iout,'(2A)')  '  Use +60 for the northern or -60 for',
     &                         ' the southern hemisphere'
           call camxerr()
         endif
C*****Pavan: added  -> psp_stdlon, psp_truelat1, psp_lon1, psp_lat1
C*****       remove -> polelon, polelat
      elseif( lpolar ) then
         read(inprec(21:),*,ERR=7002,END=7002) xorg, yorg, delx, dely,
     &                   psp_stdlon, psp_truelat1, psp_lon1, psp_lat1
C*****
      else
         read(inprec(21:),*,ERR=7002,END=7002) xorg,yorg,delx,dely,
     &                                       xlonc,ylatc,tlat1,tlat2
      endif
c
      irec = irec + 1
      action = 'Reading the time zone definition.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),*,ERR=7002,END=7002) itzon
c
      irec = irec + 1
      action = 'Reading the PiG parameters.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),*,ERR=7002,END=7002) xlmax,agemax
c
      irec = irec + 1
      action = 'Reading the number of output species.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(I10)',ERR=7002,END=7002) navspc
      if( navspc .GT. 0 ) then
         if( MOD(navspc,6) .EQ. 0 ) then
            nlines = navspc/6
         else
            nlines = INT(navspc/6) + 1
         endif
         action = 'Reading the output species names.'
         do i=1,nlines
           ibeg = (i - 1) * 6 + 1
           irec = irec + 1
           read(inp,'(A)',ERR=7001,END=7004) inprec 
           read(inprec(21:),'(6A10)',ERR=7002,END=7002) 
     &                                 (spavg(l),l=ibeg,ibeg+5)
         enddo
      endif
c
      irec = irec + 1
      action = 'Reading the number of nested grids.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(I10)',ERR=7002,END=7002) nnest
c
c-----Read fine grid parameters
c
      ngrid = nnest + 1
      if( ngrid .GT. 1 ) then
        if( ngrid .GT. MXGRID ) then
          write(iout,'(//,a)') 'ERROR in STARTUP:'
          write(iout,'(a)') 'Number of grids exceeds internal dimension'
          write(iout,'(a,i2)') 'Number of grids to be run:  ',ngrid
          write(iout,'(a,i2)') 'Maximum dimension (MXGRID): ',MXGRID
          write(iout,'(a,a)')'Increase MXGRID and expand the grid',
     &                       ' dimension lists in CAMx.PRM'
          write(iout,'(/,a,a)')'Besides MXGRID, CAMx is coded to treat',
     &                         ' up to 5 generations of nested grids'
          write(iout,'(a,a)')'You will need to modify NESTING.F to',
     &                       ' increase this capacity'
          call camxerr()
        endif
        do n = 2,ngrid
          irec = irec + 1
          write(action,'(A,I4)') 'Reading the definition for nest: ',n-1
          read(inp,'(A)',ERR=7001,END=7004) inprec 
          read(inprec(21:),*,ERR=7002,END=7002) inst1(n),inst2(n),
     &                            jnst1(n),jnst2(n),nlay(n),meshold(n)
        enddo
      endif
c
c-----Flags
c
      irec = irec + 1
      action = 'Reading the flag for advection solvers.'
      read(inp,'(A)',ERR=7001,END=7004) inprec
      advflg = inprec(21:30)
      call jstlft( advflg )
      call toupper( advflg )
      if( advflg .eq. 'SMOLAR    ' ) then
        write(iout,'(//,a)') 'ERROR in STARTUP:' 
        write(iout,'(3A)') 'Incorrect horizontal advection solver ',
     &                     'read from control file: ',inprec(21:30)
        write(iout,'(A)') 
     &          'This version of CAMx no longer supports Smolarkiewicz'
        write(iout,'(A)')'Please select BOTT or PPM'
        call camxerr()
      elseif( advflg .EQ. 'BOTT      ' ) then
        lbot = .TRUE.
        iadvct = 2
      elseif( advflg .EQ. 'PPM       ' ) then
        iadvct = 3
      else
        write(iout,'(//,a)') 'ERROR in STARTUP:' 
        write(iout,'(3A)') 'Incorrect horizontal advection solver ',
     &                     'read from control file: ',inprec(21:30)
        write(iout,'(A)')'CAMx supports BOTT and PPM'
        call camxerr()
      endif
c
      irec = irec + 1
      action = 'Reading the flag for chemistry solvers.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(A)',ERR=7002,END=7002) slvstr
      call jstlft( slvstr )
      call toupper( slvstr )
      if( slvstr .EQ. CDCMC ) then
          idsolv = IDCMC
      else if( slvstr .EQ. CDIEH ) then
          idsolv = IDIEH
      else
          write(iout,'(//,a)') 'ERROR in STARTUP:'
          write(iout,'(/,1X,6A)') 'Invalid option for chemistry ',
     &                'solver specified in control file: ',inprec(21:30)
          write(iout,'(1X,A)') 'Acceptable options are: '
          write(iout,'(10X,A)') CDCMC
          write(iout,'(10X,A)') CDIEH
          write(iout,'(1X,2A)') 'Make sure the control is updated for',
     &                        ' this version of CAMx.'  
          call camxerr()
      endif
c
      irec = irec + 1
      action = 'Reading the flag for restart.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) lrstrt 
c
      irec = irec + 1
      action = 'Reading the flag for chemstry.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) lchem
c
      irec = irec + 1
      action = 'Reading the flag for dry depostion.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) ldry 
c
      irec = irec + 1
      action = 'Reading the flag for wet depostion.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) lwet 
c
      irec = irec + 1
      action = 'Reading the flag for PiG treatment.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
c
c------------------------------------------------------------------------
c     The following block was commented out because the IRON-PIG 
c     code is not ready for distribution with version 4.00 of CAMx.
c             gwilson  Tue May 13 10:18:37 PDT 2003
c------------------------------------------------------------------------
c
cironpig      pigstr = inprec(21:30)
cironpig      call jstlft( pigstr )
cironpig      call toupper( pigstr )
cironpig      if( pigstr .eq. 'GREASD    ' ) then
cironpig        ipigflg = GRESPIG
cironpig      elseif( pigstr .EQ. 'IRON      ' ) then
cironpig        ipigflg = IRONPIG
cironpig      elseif( pigstr .EQ. 'NONE      ' ) then
cironpig        ipigflg = 0
cironpig      else
cironpig        write(iout,'(//,a)') 'ERROR in STARTUP:'
cironpig        write(iout,'(3A)') 'Incorrect response for PiG option ',
cironpig     &                     'read from control file: ',inprec(21:30)
cironpig        write(iout,'(A)')'This version of CAMx supports: ',
cironpig     &                                       'NONE, GREASD, and IRON'
cironpig        call camxerr()
cironpig      endif
c
c------------------------------------------------------------------------
c     The following block was put in to default to the old option
c     PiG on or off.  The IRON-PIG is not ready for distribution with
c     version 4.
c             gwilson  Tue May 13 10:18:37 PDT 2003
c------------------------------------------------------------------------
c
      read(inprec(21:),'(L10)',ERR=7002,END=7002) lpig
      ipigflg = 0
      if( lpig ) ipigflg = GRESPIG
c
c------------------------------------------------------------------------
c

      irec = irec + 1
      action = 'Reading the flag for staggered winds.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) lstagw
c
      irec = irec + 1
      action = 'Reading the flag for treating area emissions.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) larsrc
c
      irec = irec + 1
      action = 'Reading the flag for treating point emissions.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) lptsrc
      if( (ipigflg .EQ.  GRESPIG .OR. ipigflg .EQ. IRONPIG) .AND.
     &                                          .NOT. lptsrc) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,'(/,a,/)')
     &         'Point source emissions are required for PiG module.'
        write(iout,*)'Either turn off PiG module or turn on ',
     &                                      'point sources treatment.'
        call camxerr()
      endif
c
      irec = irec + 1
      action = 'Reading the flag for 1-day emissions files.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) le1day
c
      irec = irec + 1
      action = 'Reading the flag for 3-D average files.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) l3davg
c
c======================== Source Apportion Begin =======================
c
c  --- call routine to read the source apportionment specific user options ---
c
      lddm = .FALSE.
      lproca = .FALSE.
      irec = irec + 1
      action = 'Reading the flag for probing tools.'
      read(inp,'(A)',ERR=7001,END=7004) inprec 
      read(inprec(21:),'(L10)',ERR=7002,END=7002) ltrace
      if( ltrace ) then
c
c   --- Output file root name
c
         irec = irec + 1
         action = 
     &     'Reading the output filename root for probing tools.'
         read(inp,'(A)',ERR=7001,END=7004) inprec
         flrtsa = inprec(21:)
         call jstlft( flrtsa )
c
c   --- flag for probing tools technology type ---
c
         irec = irec + 1
         action = 
     &     'Reading the technology type for probing tools.'
         read(inp,'(A)',ERR=7001,END=7002) inprec
         tecsav = inprec(21:30)
         tectyp = tecsav
         call jstlft( tectyp )
         call toupper( tectyp )
         lfound = .FALSE.
         if( tectyp .EQ. OSAT  ) then
             lfound = .TRUE.
             verson = VEROSAT
         else if( tectyp .EQ. GOAT  ) then
             lfound = .TRUE.
             verson = VERGOAT
         else if( tectyp .EQ. APCA  ) then
             lfound = .TRUE.
             verson = VERAPCA
         else if( tectyp .EQ. RTRAC ) then
             lfound = .TRUE.
             verson = VERRTRAC
         else if( tectyp .EQ. DDM  ) then 
             ltrace = .FALSE.
             lddm = .TRUE.
             lfound = .TRUE.
             verson = VERDDM
         else if( tectyp .EQ. STRPA ) then
             ltrace = .FALSE.
             lddm = .FALSE.
             lfound = .TRUE.
             lproca = .TRUE.
             lipr = .TRUE.
             lirr = .TRUE.
             lsfcfl(IDXCRS) = .TRUE.
             if( ngrid .GT. 1 ) lsfcfl(IDXFIN) = .TRUE.
             verson = VERPA
         else if( tectyp .EQ. STRIPR ) then
             ltrace = .FALSE.
             lddm = .FALSE.
             lfound = .TRUE.
             lproca = .TRUE.
             lipr = .TRUE.
             lirr = .FALSE.
             verson = VERPA
         else if( tectyp .EQ. STRIRR ) then
             ltrace = .FALSE.
             lddm = .FALSE.
             lfound = .TRUE.
             lproca = .TRUE.
             lipr = .FALSE.
             lirr = .TRUE.
             lsfcfl(IDXCRS) = .TRUE.
             if( ngrid .GT. 1 ) lsfcfl(IDXFIN) = .TRUE.
             verson = VERPA
         endif
         if( .NOT. lfound ) goto 7003
c
c  --- call routine to call the rest of the options ---
c
         if( ltrace ) then
            if( tectyp .NE. RTRAC ) then
               call rdoptsa(inp,ctlfil,irec)
            endif
         else if( lddm ) then
            call rdoptddm(inp,ctlfil,irec)
         else if( lproca ) then
            call rdoptpa(inp,ctlfil,irec)
         endif
      endif
c
c========================= Source Apportion End ========================
c
c-----Perform checks on run control parameters
c
c-----Date/time
c
      if( enddate .LT. begdate ) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*) 'Simulation end date is less than ',
     &                'simulation start date.'
        write(iout,'(a,i10.5)') 'Simulation start: ',begdate
        write(iout,'(a,i10.5)') 'Simulation end:   ',enddate
        call camxerr()
      elseif( enddate .EQ. begdate ) then
        if( endtim .LE. begtim ) then
          write(iout,'(//,a)') 'ERROR in STARTUP:'
          write(iout,*) 'Simulation end time is less than ',
     &                  'simulation start time.'
          write(iout,*) 'Simulation start: ',begtim
          write(iout,*) 'Simulation end:   ',endtim
          call camxerr()
        endif
      endif
c
c-----Check I/O frequency
c
c-----The first check is to help User's migrate to specifying
c     frequencies in minutes starting with CAMx3.1.
c
      if( dtmax .LE. 1 .OR. dtinp .LE. 1. .or. dtems .LE. 1. 
     &                                     .or. dtout .LE. 1.) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,'(/,a,/)')
     &            'I/O frequencies must be greater than 1 minute'
        write(iout,*)'Maximum time step read  (DTMAX): ',dtmax
        write(iout,*)'Input interval read     (DTINP): ',dtinp
        write(iout,*)'Emissions interval read (DTEMS): ',dtems
        write(iout,*)'Output interval read    (DTOUT): ',dtout
        write(iout,'(/,2a)')
     &            'Starting with version 3.1 these values must be in ',
     &            'minutes rather than hours.',
     &            'Do you need to convert these inputs to minutes?'
        call camxerr()
      endif
c
c-----Really important checks follow
c
      if( dtmax .GT. MAXDT ) then
        write(iout,'(//,A)') 'WARNING in STARTUP:'
        write(iout,'(A)') 
     &               'Invalid value specified for maximum time step.'
        write(iout,'(A)') 'A default value will be used instead.'
        write(iout,*)'Maximum time step read (DTMAX): ',dtmax
        write(iout,*)'Default value to be used      : ',MAXDT
        dtmax = MAXDT
      endif
c 
      if( (dtmax .GT. 60. .AND. amod(dtmax,60.) .GT. 0. ) .OR.
     &    (dtinp .GT. 60. .AND. amod(dtinp,60.) .GT. 0. ) .OR.
     &    (dtems .GT. 60. .AND. amod(dtems,60.) .GT. 0. ) .OR.
     &    (dtout .GT. 60. .AND. amod(dtout,60.) .GT. 0. )) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*)'An Input/Output interval is > 60 minutes.'
        write(iout,*)'It must be an integer multiple of 60 minutes.'
        write(iout,*)'Maximum time step read  (DTMAX): ',dtmax
        write(iout,*)'Input interval read     (DTINP): ',dtinp
        write(iout,*)'Emissions interval read (DTEMS): ',dtems
        write(iout,*)'Output interval read    (DTOUT): ',dtout
        call camxerr()
      endif
      if( (dtmax .LT. 60. .AND. amod(60.,dtmax) .GT. 0. ) .OR.
     &    (dtinp .LT. 60. .AND. amod(60.,dtinp) .GT. 0. ) .OR.
     &    (dtems .LT. 60. .AND. amod(60.,dtems) .GT. 0. ) .OR.
     &    (dtout .LT. 60. .AND. amod(60.,dtout) .GT. 0. )) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*)'An Input/Output interval is < 60 minutes.'
        write(iout,*)'It must divide 60 minutes evenly.'
        write(iout,*)'Maximum time step read  (DTMAX): ',dtmax
        write(iout,*)'Input interval read     (DTINP): ',dtinp
        write(iout,*)'Emissions interval read (DTEMS): ',dtems
        write(iout,*)'Output interval read    (DTOUT): ',dtout
        call camxerr()
      endif
      if( amod(amax1(dtinp,dtems),amin1(dtinp,dtems)) .GT. 0. .OR.
     &    amod(amax1(dtinp,dtout),amin1(dtinp,dtout)) .GT. 0. .OR.
     &    amod(amax1(dtems,dtout),amin1(dtems,dtout)) .GT. 0. ) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*)'Input/Output intervals must be even multiples'
        write(iout,*)'of each other.'
        write(iout,*)'Input interval read (DTINP)    : ',dtinp
        write(iout,*)'Emissions interval read (DTEMS): ',dtems
        write(iout,*)'Output interval read (DTOUT)   : ',dtout
        call camxerr()
      endif
c
c-----Coarse grid parameters
c
      if( ncol(1) .GT. MXCOLA .OR. nrow(1) .GT. MXROWA .OR.
     &                                       nlay(1) .GT. MXLAYA ) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*) 'Input grid specifications greater than ',
     &                'internal model dimensions.'
        write(iout,'(10X,A,5X,A)') 'Parameters','Grid Definition'
        write(iout,'(A10,5X,I5,10X,I5)') 'Rows    :',MXROWA,nrow(1)
        write(iout,'(A10,5X,I5,10X,I5)') 'Columns :',MXCOLA,ncol(1)
        write(iout,'(A10,5X,I5,10X,I5)') 'Layers  :',MXLAYA,nlay(1)
        call camxerr()
      endif
c
c-----Fine grid parameters
c
      if( ngrid .gt. 1 ) then
        do n = 2,ngrid
          if( inst1(n) .LE. 1 .OR. inst1(n) .GE. ncol(1) .OR.
     &        inst2(n) .LE. 1 .OR. inst2(n) .GE. ncol(1) .OR.
     &        jnst1(n) .LE. 1 .OR. jnst1(n) .GE. nrow(1) .OR.
     &        jnst2(n) .LE. 1 .OR. jnst2(n) .GE. nrow(1) ) then
              write(iout,'(//,a)') 'ERROR in STARTUP:'
              write(iout,*) 'For fine grid # ',n
              write(iout,*) 'Starting/ending indices exceed extent',
     &                      ' of coarse grid.'
              write(iout,*) 'Nest column range   :',inst1(n),inst2(n)
              write(iout,*) 'Coarse column range :',1,ncol(1)
              write(iout,*) 'Nest row range      :',jnst1(n),jnst2(n)
              write(iout,*) 'Coarse row range    :',1,nrow(1)
              call camxerr()
          endif
c
c   --- calculate the dimensions for the grid and check for
c       array overflow ---
c
          ncol(n) = (inst2(n) - inst1(n) + 1)*meshold(n) + 2
          nrow(n) = (jnst2(n) - jnst1(n) + 1)*meshold(n) + 2
          if( ncol(n) .GT. MXCOLA .OR. nrow(n) .GT. MXROWA .OR.
     &                                        nlay(n) .GT. MXLAYA ) then
            write(iout,'(//,a)') 'ERROR in STARTUP:'
            write(iout,*) 'For fine grid # ',n
            write(iout,*) 'Input grid specifications greater than ',
     &                    'internal model dimensions.'
            write(iout,'(10X,A,5X,A)') 'Parameters','Grid Definition'
            write(iout,'(A10,5X,I5,10X,I5)') 'Rows    :',MXROWA,nrow(n)
            write(iout,'(A10,5X,I5,10X,I5)') 'Columns :',MXCOLA,ncol(n)
            write(iout,'(A10,5X,I5,10X,I5)') 'Layers  :',MXLAYA,nlay(n)
            call camxerr()
          endif
        enddo
      endif
c
c-----Call routine to set up the pointers for grid vectors
c
      call iniptr()
c 
c-----PiG requires Chemistry 
c 
      if( .NOT. lchem .AND. ipigflg .NE. 0 ) then
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*) 'PiG treatment requires Chemistry, but..' 
cironpig        write(iout,*) 'PiG flag read              : ',pigstr
        write(iout,*) 'PiG flag read       (LPIG) : ',lpig
        write(iout,*) 'Chemistry flag read (LCHEM): ',lchem 
        call camxerr()
      endif
c
c-----Echo run control parameters
c
      write(idiag,'(A,F10.0,I10.6,I10.5)')
     &         'Simulation start time/date : ',begtim,cbdate,begdate
      write(idiag,'(A,F10.0,I10.6,I10.5)')
     &         'Simulation end time/date   : ',endtim,cedate,enddate
      write(idiag,'(A,F10.0)')
     &         'Max timestep (min)         : ',dtmax
      write(idiag,'(A,F10.0)')
     &         'Met Input interval (min)   : ',dtinp
      write(idiag,'(A,F10.0)')
     &         'Emiss Input interval (min) : ',dtems
      write(idiag,'(A,F10.0)')
     &         'Output interval (min)      : ',dtout
      write(idiag,'(A,3I10)')
     &         'NCOL NROW NLAY             : ',ncol(1),nrow(1),nlay(1)
      write(idiag,'(A,A)')
     &         'Grid Projection Type       : ',grdflg
      write(idiag,'(A,4F10.3)')
     &         'X/Y Origin, cell size      : ',xorg,yorg,delx,dely
      if (lutm) then
        write(idiag,'(A,I10)')
     &         'UTM: zone number           : ',iuzon
      elseif (lpolar) then
        write(idiag,'(A,2F10.3)')
     &         'Polar: pole lon/lat        : ',psp_stdlon, psp_truelat1, psp_lon1, psp_lat1
      elseif (lambrt) then
        write(idiag,'(A,4F10.3)')
     &         'Lmbrt: lonc/latc, true lats: ',xlonc,ylatc,tlat1,tlat2
      endif
      write(idiag,'(A,I10)')
     &         'Time zone                  : ',itzon
cironpig      write(idiag,'(A,A)')
cironpig     &         'PiG submodel               : ',pigstr
      write(idiag,'(A,L10)')
     &         'PiG submodel               : ',lpig
      if( ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG ) then
        write(idiag,'(A,F10.0)')
     &         'Maximum PiG puff length (m): ',xlmax
        write(idiag,'(A,F10.0)')
     &         'Maximum PiG age (hr)       : ',agemax
        agemax = 3600.*agemax
      endif
      write(idiag,'(A,A)')
     &         'Advection Solver           : ',advflg
      write(idiag,'(A,A)')
     &         'Chemistry Solver           : ',slvstr
      write(idiag,*)
      if( navspc .GT. 0 ) then
        write(idiag,'(A,I3)')'Number of output average species: ',navspc
        write(idiag,'(8A10)') (spavg(l),l=1,navspc)
      else
        write(idiag,*)'Average concentrations will not be output.'
      endif
      write(idiag,*)
c
      if( nnest .GT. 0 ) then
        write(idiag,'(A,I3)')'Number of nested fine grids: ',nnest
        write(idiag,*)' Nest        x-range       ncol        y-range',
     &                '       nrow   mesh factor  nlay'      
        do n = 2,ngrid
          write(idiag,'(I5,2(5X,I5,1X,I5,5X,I5),5X,I5,5X,I5)')
     &          n,inst1(n),inst2(n),ncol(n),jnst1(n),jnst2(n),
     &          nrow(n),meshold(n),nlay(n)
        enddo
        write(idiag,*)
        write(idiag,*) '|',('-',i=1,74),'|'
        write(idiag,*) '|',(' ',i=1,74),'|'
        write(idiag,'(2A)') ' | NOTE:  The nest order listed above ',
     &       'is the original order specified in    |'
        write(idiag,'(2A)') ' |        the CAMx.in file.           ',
     &       '                                      |'
        write(idiag,'(2A)') ' |        CAMx may re-order the nests.',
     &      ' See the internal nest order provided |'
        write(idiag,'(2A)') ' |        in the table below.         ',
     &       '                                      |'
        write(idiag,*) '|',(' ',i=1,74),'|'
        write(idiag,*) '|',('-',i=1,74),'|'
      else
        write(idiag,*)'Fine grid nests are not specified.'
      endif
      write(idiag,*)
c
      write(idiag,*)'CAMx control flags'
      write(idiag,'(A,L10)') 'Restart             :',lrstrt 
      write(idiag,'(A,L10)') 'Chemistry           :',lchem
      write(idiag,'(A,L10)') 'Dry deposition      :',ldry 
      write(idiag,'(A,L10)') 'Wet deposition      :',lwet 
      write(idiag,'(A,L10)') 'Staggered winds     :',lstagw
      write(idiag,'(A,L10)') 'Area sources        :',larsrc
      write(idiag,'(A,L10)') 'Point sources       :',lptsrc
      write(idiag,'(A,L10)') '1-day emiss files   :',le1day
      write(idiag,'(A,L10)') '3-D average file    :',l3davg
c 
c======================== Source Apportion Begin =======================
c
c   --- echo the source apportionment specific flags ---
c
      if( .NOT. ltrace .AND. .NOT. lddm .AND. .NOT. lproca ) then
         write(idiag,'(A,L10)') 'Probing Tools       :',.FALSE.
      else
         write(idiag,'(A,L10)') 'Probing Tools       :',.TRUE.
         write(idiag,'(/,A,/)') '           Probing Tool Parameters'
         write(idiag,'(A,A10)') 'Technology type     :',tectyp
         if( lddm ) then
             write(idiag,'(A,I10)') 'Number of IC groups :',nicddm
             if( nicddm .GT. 0 ) then
                 write(idiag,'(A)') 'DDM IC species used :'
                 do i=1,nicddm
                    write(idiag,'(22X,A)') icddmsp(i)
                 enddo
             else
                 write(idiag,'(22X,A)') 'No IC species modeled.'
             endif 
             write(idiag,'(A,I10)') 'Number of BC groups :',nbcddm
             if( nbcddm .GT. 0 ) then
                 write(idiag,'(A)') 'DDM BC species used :'
                 do i=1,nbcddm
                    write(idiag,'(22X,A)') bcddmsp(i)
                 enddo
             else
                 write(idiag,'(22X,A)') 'No BC species modeled.'
             endif 
         endif     
         if( tectyp .NE. RTRAC ) write(idiag,'(A,L10)') 
     &                                'Stratify Boundary   :',lbndry
         if( lddm ) then
             write(idiag,'(A,I10)') 'Number of EM groups :',nemddm
             if( nemddm .GT. 0 ) then
                 write(idiag,'(A)') 'DDM EM species used :'
                 do i=1,nemddm
                    write(idiag,'(22X,A)') emddmsp(i)
                 enddo
             else
                 write(idiag,'(22X,A)') 'No emissions species modeled.'
             endif 
         endif
         if( tectyp .NE. RTRAC ) then
             write(idiag,'(A,I10)') 'No. of source areas :',nregin
c
c   --- set the number of emissions groups from the source groups ----
c
             if( leftovr ) then
                 nemiss = ngroup + 1
             else
                 nemiss = ngroup
             endif
             write(idiag,'(A,I10)') 'No. of source groups:',nemiss
             if( ltrace ) then
                write(idiag,'(A,L10)') 'Use Leftover group  :',leftovr
                write(idiag,'(A,I10)') 'No. of time releases:',ntrtim
             endif
         endif 
      endif
      write(idiag,*)
c
c========================= Source Apportion End ========================
c
c
c======================== Process Analysis Begin ====================================
c
c
c   --- echo the irmb specific flags ---
c
      if( lproca ) then
         write(idiag,'(A,L10,/)/')
     &                  'Integrated Process Rates Analysis :',lipr
         write(idiag,'(A,L10,//)')
     &                  'Integrated Reaction Rates Analysis:',lirr
      endif
c
c========================= Process Analysis End =====================================
c
c
c-----Open the remaining I/O files
c
      call openfils(inp,iifroot,filroot,ctlfil,irec)
c
c-----Initialize simulation clock
c
      time = begtim
      date = begdate
c
c-----Read chemistry files
c
      call chmprep
c
c----Prepare met file for netCDF
c
      call metprep
c             
c-----Check number of average output species 
c             
      if (navspc.gt.nspec) then 
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*)'Number of average species to be output ', 
     &               'is greater than number of species to model' 
        write(iout,*)'NAVSPC,NSPEC: ',navspc,nspec 
        call camxerr()
      endif
c  
c-----Map average species to model species
c  
      do 46 lav = 1,navspc
        do l = 1,nspec  
          if (spavg(lav).eq.spname(l)) then  
            lavmap(lav) = l
            write(idiag,'(2(a,i5,2x,a))')
     &                   'Average species ',lav,spavg(lav),
     &                   ' mapped to modeled species ',l,spname(l)
            goto 46  
          endif  
        enddo  
        write(iout,'(//,a)') 'ERROR in STARTUP:'
        write(iout,*)'Did not find average species: ',spavg(lav)  
        write(iout,*)'in species list'
        call camxerr()
 46   continue
      write(idiag,*)
c
c======================== Source Apportion Begin =======================
c
      if( ltrace ) then
c
c  --- check that switches are compatible ---
c
         if( (ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG) 
     &                                  .AND. tectyp .EQ. RTRAC ) then
             write(iout,'(//,a)') 'ERROR in STARTUP:'
             write(iout,*) 'Cannot use PiG submodel with RTRAC source ',
     &                                                 'apportionment.'
             write(iout,*) '  Turn off PiG switch and try again.'
             call camxerr()
         endif
c
c   --- call routine to initialize the source apportionment
c       data structures ---
c
          call initsa(version,ncol,nrow,
     &                             begdate,begtim,enddate,endtim)
c
c   --- if this is a restart, call routine to read the instantaneous files ---
c
          if( lrstrt ) then
             call rdinstsa(1,begdate,time,ncol(1),nrow(1),nlay(1),
     &                     ntotsp,ptconc(1))
          endif
c
c  ---- call routine to calculate the average reactivity
c       of boundary conditions ---
c
          if( tectyp .NE. RTRAC ) then
              call clcbwt(begdate,begtim,enddate,endtim,
     &                                        ncol(1),nrow(1),nlay(1))
          endif
c
c   --- call routine to write the header of the coarse grid
c       average surface tracer concentrations file ---
c
          if( lsfcfl(IDXCRS) ) then
             call hdrwsa(IOWSFC+IDXCRS,sfcfil(IDXCRS),'AVERAGE   ',
     &                        ntotsp,1,begdate,begtim,enddate,endtim)
          endif
      endif
c
c========================= Source Apportion End ========================
c
c
c======================== DDM Begin ====================================
c
      if( lddm ) then
c
c  --- check that switches are compatible ---
c
         if( naero .GT. 0 ) then
             write(iout,'(//,a)') 'ERROR in STARTUP:'
             write(iout,*) 'Cannot use DDM with Aerosol chemistry.'
             write(iout,*) 
     &          '  Turn off DDM or use a different chemical mechanism.'
             call camxerr()
         endif

         if( ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG ) then
             write(iout,'(//,a)') 'ERROR in STARTUP:'
             write(iout,*) 'Cannot use PiG submodel with DDM.'
             write(iout,*) '  Turn off PiG switch and try again.'
             call camxerr()
         endif
c
c  --- check that chemical mechanism is valid ---
c
         if( idmech .LT. 2 .OR. idmech .GT. 5 ) then
            write(iout,'(//,a)') 'ERROR in STARTUP:'
            write(iout,*) 'Chemical mechanism is invalid for ',
     &                                               'for use with DDM.'
            write(iout,'(A,I2)') ' Mechanism ID: ',idmech
            write(iout,*) '  Please supply a DDM compatible chemistry ',
     &                                             'file and try again.'
            call camxerr()
         endif
c
c  --- check that advection solver is correct ---
c
         if( .NOT. lbot ) then
            write(iout,'(//,a)') 'ERROR in STARTUP:'
            write(iout,*) 'Advection solver used with DDM ',
     &                                                   'must be BOTT.'
            write(iout,*)'Please change the control file and try again.'
            call camxerr()
          endif
c
c   --- call routine to initialize the source apportionment
c       data structures ---
c
         call initsa(version,ncol,nrow,begdate,begtim,enddate,endtim)
c
c   --- if this is a restart, call routine to read the instantaneous files ---
c
         if( lrstrt ) then
             call rdinstsa(1,begdate,time,ncol(1),nrow(1),nlay(1),
     &                     ntotsp,ptconc(1))
         endif
c
c   --- call routine to write the header of the coarse grid
c       average surface tracer concentrations file ---
c
         if( lsfcfl(IDXCRS) ) then
             call hdrwsa(IOWSFC+IDXCRS,sfcfil(IDXCRS),'AVERAGE   ',
     &                        ntotsp,1,begdate,begtim,enddate,endtim)
         endif
      endif
c
c======================== DDM End ====================================
c
c-----Read BC file header and set irregular boundary cells
c
      call bndprep(begtim,begdate,endtim,enddate)
c
c-----Calculate grid parameters for coarse grid
c
      call grdprep(ncol(1),nrow(1),cellon(1),cellat(1),mapscl(1))

c-----Pavan (Pavan_Nandan_Racherla@alumni.cmu.edu)
c-----Updated Oct 10 2008
      write(iout, *) "Pavan: print lon/lat's for outermost domain"
      do domlen_vec2d = 1, mxvec2d
         write(iout, '(2(E14.6,1X))') cellon(domlen_vec2d), cellat(domlen_vec2d)
      end do
c
c-----Calculate nested grid mapping parameters
c
      if (ngrid.gt.1) then
        call nstprep
      else
        do j=1,nrow(1)
          do i=1,ncol(1)
            n = (j-1)*ncol(1) + i
            idfin(n) = 0
          enddo
        enddo
        mapgrd(1) = 1
        nchdrn(1) = 0
        meshold(1) = 1
        nmesh(1) = 1
      endif
c
c======================== Source Apportion Begin =======================
c
      if( ltrace ) then
          if( tectyp .EQ. RTRAC ) then
c
c   --- get each emissions files to the correct place ----
c
              do i=1,ngrid
                  call empreprt(begdate,begtim,i)
              enddo
c
c   ---- call routine to read the receptor definition file ---
c
              call rdrcprt()
c
c  --- call routine to write the header of receptor average file ---
c
              call hdrcprt(begdate,begtim,enddate,endtim)
c
c  ---- call routine to calculate the average reactivity
c       of initial conditions ---
c
          else
              call clciwt(begdate,begtim,enddate,endtim,ly2k,
     &                ncol(1),nrow(1),nlay(1))
c
c   --- get each emissions files to the correct place ----
c
              do i=1,ngrid
                  call emprepsa(begdate,begtim,i)
              enddo
c
c  ---- call routine to calculate the average reactivity
c       of emissions ---
c
              call clcewt(enddate,endtim)
c
c   --- get each emissions files BACK to the correct place ----
c
              do i=1,ngrid
                  call emprepsa(begdate,begtim,i)
              enddo
c
c   ---- call routine to read the receptor definition file ---
c
              call rercp()
c
c  --- call routine to write the header of receptor average file ---
c
              call hdrrcp(begdate,begtim,enddate,endtim)
          endif
      endif
c
c========================= Source Apportion End ========================
c
c
c======================== DDM Begin ====================================
c
          if( lddm ) then
c
c   ---- call routine to read the receptor definition file ---
c
              call rercp()
c
c  --- call routine to write the header of receptor average file ---
c
              call hdrrcpddm(begdate,begtim,enddate,endtim)
          endif
c
c======================== DDM End ====================================
c
c
c============================= Process Analysis Begin ==============================
c
c-----Call routine to initialize Process Analysis subdomains
c
      if( lproca ) then
         call pagrids()
c
c-----Call routine to write the header to the otuput files ---
c
         if( lipr ) call wrtiprhdr(begdate,begtim,enddate,endtim)
         if( lirr ) then
            call pasetup
            call wrtirrhdr(begdate,begtim,enddate,endtim)
            if( lsfcfl(IDXCRS) ) then
                if( l3davg ) then
                    nlayer = nlay(1) 
                else
                    nlayer = 1
                endif
                call hdrwsa(IOWSFC+IDXCRS,sfcfil(IDXCRS),'AVERAGE   ',
     &                     ntotsp,nlayer,begdate,begtim,enddate,endtim)
            endif
         endif
      endif
c
c============================= Process Analysis End ==============================
c
c-----Assign albedo values for fine grids
c
      if (lchem) then
        do ip = 1,ngrid
          do ic = 1,nchdrn(ip)
            igrd = idchdrn(ic,ip)
            if( .NOT. lrdalb(igrd) ) call iassgn2d(
     &                    ncol(ip),nrow(ip),i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    icdalb(iptr2d(ip)),icdalb(iptr2d(igrd)))
          enddo
        enddo
      endif
c
c------Calculate grid parameters for fine grids
c
      do ip = 1,ngrid
        do ic = 1,nchdrn(ip)
          igrd = idchdrn(ic,ip)
          call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                  nmesh(igrd),ncol(igrd),nrow(igrd),
     &                  cellat(iptr2d(ip)),cellat(iptr2d(igrd)) )
          call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                  nmesh(igrd),ncol(igrd),nrow(igrd),
     &                  cellon(iptr2d(ip)),cellon(iptr2d(igrd)) )
          call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                  nmesh(igrd),ncol(igrd),nrow(igrd),
     &                  mapscl(iptr2d(ip)),mapscl(iptr2d(igrd)) )
        enddo
      enddo
c
c-----Read height/pressure for all grids, and winds & temps for all available
c     grids, to current time
c
      do igrd=1,ngrid
        call readpi0(igrd,ncmet,ncol(igrd),nrow(igrd),
     &               nlay(igrd),pi0(iptr3d(igrd)) )
        call readzpwt(igrd,ihtp(igrd),iwind(igrd),itemp(igrd),iout,
     &                ncol(igrd),nrow(igrd),nlay(igrd),
     &                height(iptr3d(igrd)),
     &                press(iptr3d(igrd)),depth(iptr3d(igrd)),
     &                windu(iptr3d(igrd)),
     &                windv(iptr3d(igrd)),tempk(iptr3d(igrd)),
     &                tsurf(iptr2d(igrd)) )
      enddo
c
c-----Determine vertical meshing number
c
      do ip = 1,ngrid
        do ic = 1,nchdrn(ip)
          ig = idchdrn(ic,ip)
c
c-----If ZP file is not provided, make sure that there is
c     no vertical meshing----
c
          if( ihtp(ig) .LE. 0 ) then
              if( nlay(ig) .NE. nlay(ip) ) then
                 write(iout,'(//,a)') 'ERROR in STARTUP:'
                 write(iout,*) 'A height/pressure file is required ',
     &                                          'for vertical nesting.'
                 write(iout,'(1X,A,I2,A)') 
     &                 'Set the number of layers for grid:',ig,
     &                        ' to the same value as its parent grid.'
                 write(iout,*) 'OR supply a height/pressure file.'
                 call camxerr()
               endif
          endif
c
          call vnmshcal(ig,ncol(ip),nrow(ip),nlay(ip),i1(ig),j1(ig),
     &                  ncol(ig),nrow(ig),nlay(ig),height(iptr3d(ip)),
     &                  height(iptr3d(ig)),nmshv(1,ig))
c
c-----If no ZP file is supplied call routines to 
c     interpolate ZP from parent grid---
c
          if( ihtp(ig) .LE. 0 ) then
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Assigning heights from parent grid',
     &                             time, date,' grid',ig
            call rassgn3d(ncol(ip),nrow(ip),nlay(ip),
     &           i1(ig),j1(ig),nmesh(ig),ncol(ig),nrow(ig),
     &                          height(iptr3d(ip)),height(iptr3d(ig)) )
            call rassgn3d(ncol(ip),nrow(ip),nlay(ip),
     &           i1(ig),j1(ig),nmesh(ig),ncol(ig),nrow(ig),
     &                          depth(iptr3d(ip)),depth(iptr3d(ig)) )
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating pressure from parent grid',
     &                             time, date,' grid',ig
            call interp2d(ncol(ip),nrow(ip),nlay(ip),
     &             i1(ig),j1(ig),nmesh(ig),ncol(ig),nrow(ig),
     &                            press(iptr3d(ip)),press(iptr3d(ig)) )
            call interpv(ncol(ip),nrow(ip),nlay(ip),ncol(ig),
     &                   nrow(ig),nlay(ig),nmesh(ig),
     &                   nmshv(1,ig),i1(ig),j1(ig),
     &                   height(iptr3d(ip)),
     &                   height(iptr3d(ig)),press(iptr3d(ig)) )
          endif
c
          write(idiag,*)
          write(idiag,'(a,i2)')'Vertical meshing factors for grid: ',ig
          kg1 = 1
          do kp = 1,nlay(ip)
            if (nmshv(kp,ig).lt.1) then
              write(iout,'(//,a)') 'ERROR in STARTUP:'
              write(iout,*) 'Vertical meshing factor < 1!'
              write(iout,*) 'In parent layer ',kp
              call camxerr()
            endif
            do kg = kg1,kg1+nmshv(kp,ig)-1
              if (nmshv(kp,ig).gt.1) then
                write(idiag,'(a,i2,a,i2)') 'Layer ',kg,
     &                                     '   is in parent layer ',kp
              else
                write(idiag,'(a,i2,a,i2)') 'Layer ',kg,
     &                                     ' matches parent layer ',kp
              endif
            enddo
            kg1 = kg1 + nmshv(kp,ig)
          enddo
          write(idiag,*)
        enddo
      enddo
c
c======================== Source Apportion Begin =======================
c
c
c   --- if this is a restart, call routine to read the
c       instantaneous files for fine grids ---
c
      if( (ltrace .OR. lddm) .AND. lrstrt .AND. ngrid .GT. 1 )
     &                                       call rdfgsa(begdate,time)
      
c
c======================== Source Apportion End =======================
c
c
c-----Interpolate wind fields that were not read to fine grids
c
      do ip = 1,ngrid 
        do ic = 1,nchdrn(ip) 
          igrd = idchdrn(ic,ip)
          iunit = iwind(igrd) 
          call finwind(iunit,ncol(ip),nrow(ip),nlay(ip),i1(igrd),
     &                 j1(igrd),nmesh(igrd),nmshv(1,igrd),ncol(igrd),
     &                 nrow(igrd),nlay(igrd),
     &                 windu(iptr3d(ip)),windv(iptr3d(ip)),
     &                 windu(iptr3d(igrd)),windv(iptr3d(igrd)),
     &                 igrd,date,time,iout ) 
c 
c-----Interpolate temperature fields to each grid
c 
          if (itemp(igrd).eq. 0) then
            call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    tempk(iptr3d(ip)),
     &                    tempk(iptr3d(igrd)) )
            call interpv(ncol(ip),nrow(ip),nlay(ip),ncol(igrd),
     &                   nrow(igrd),nlay(igrd),nmesh(igrd),
     &                   nmshv(1,igrd),i1(igrd),j1(igrd),
     &                   height(iptr3d(ip)),
     &                   height(iptr3d(igrd)),tempk(iptr3d(igrd)))
            call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd), 
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    tsurf(iptr2d(ip)),tsurf(iptr2d(igrd)) ) 
          endif
        enddo
      enddo
c
c-----Read surface files and initialize arrays
c
      if (ldry) then
        do igrd=1,ngrid
          call srfprep(igrd,ncol(igrd),nrow(igrd),fsurf(iptrlu(igrd)) )
        enddo
      endif
c
c-----Assign fine grid landuse fractions
c
      do ip = 1,ngrid
        do ic = 1,nchdrn(ip)
          igrd = idchdrn(ic,ip)
          if (isurf(igrd).eq.0) then
            call luassgn(ncol(ip),nrow(ip),NLU,i1(igrd),j1(igrd),
     &                   nmesh(igrd),ncol(igrd),nrow(igrd),
     &                   fsurf(iptrlu(ip)),fsurf(iptrlu(igrd)) )
          endif
        enddo
      enddo
c
c-----Read emission file headers
c
      if (lptsrc) call pntprep(begtim,begdate,endtim,enddate)
      if (larsrc) then
        do igrd = 1,ngrid
          dxmod = delx
          dymod = dely
          if (igrd.gt.1) then
            dxmod = delx/meshold(igrd)
            dymod = dely/meshold(igrd)
          endif
          call areaprep(igrd,begtim,begdate,endtim,enddate,
     &                  iarem(1,igrd),iout,idiag,dxmod,dymod)
        enddo
      endif
c
c======================== DDM Begin ====================================
c
      if( lddm ) then
c
c   --- get each emissions files to the correct place,
c       NOTE:  Done here for DDM because we need the point locations
c       from regular model ----
c
          do i=1,ngrid
              call emprepsa(begdate,begtim,i)
          enddo
      endif
c
c======================== DDM End ====================================
c
c-----Read IC or restart files headers and write output concentration
c     file headers 
c 
      call cncprep(endtim,enddate) 
c
c-----Write deposition output file headers
c
      if (ldry .or. lwet) call depprep(endtim,enddate)
c
c-----Initialize PiG submodel
c
      if( ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG ) 
     &                             call pigprep(begdate,begtim)
c
c-----Determine times and dates for next inputs/emissions/write
c
      inptim = time
      inpdate = date
      emstim = time 
      emsdate = date
      bndtim = time 
      bnddate = date
      haztim = time 
      hazdate = date
      ozntim = time 
      ozndate = date
      whr = aint(time/100.)
      wmn = amod(time,100.)
      wrttim = 100.*(whr + aint((wmn + dtout)/60.)) +
     &             amod((wmn + dtout),60.)
      wrtdate = date
      if (wrttim.ge.2400.) then 
        wrttim = wrttim - 2400. 
        wrtdate = wrtdate + 1 
        if( MOD(wrtdate,1000) .GT. 365 ) then
           if( MOD(INT(wrtdate/1000),4) .EQ. 0 ) then
              if( MOD(wrtdate,1000) .EQ. 367 )
     &                   wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
           else
              wrtdate = (INT(wrtdate/1000)+1)*1000 + 1
           endif
        endif
      endif
c
c  --- everything worked correctly, return to calling routine
c
      call flush(iout)
      call flush(idiag)
      return
c
c  --- Error messages
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in STARTUP:'
      write(iout,'(2A)') ' ERROR: Could not open file: ',
     &                                   filroot(:istrln(filroot))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in STARTUP:'
      write(iout,'(2A,I2)') ' ERROR: Reading control file ( ',
     &                     ctlfil(:istrln(ctlfil)),' ) at line: ',irec
      write(iout,'(A)') action(:istrln(action))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in STARTUP:'
      write(iout,'(3A,I2)') ' ERROR: Reading control file ( ',
     &                     ctlfil(:istrln(ctlfil)),' ) at line: ',irec
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') '     Last line read:'
      write(iout,'(A)') inprec(:istrln(inprec))
      if( irec .EQ. 3 ) then
         write(iout,'(/,2A,/)') '**** Are you certain you are using a ',
     &                 'control file for CAMx version 4.01 ****'
      endif
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in STARTUP:'
      write(iout,'(/,1X,6A)') 'ERROR:  Invalid technology type ',
     &                            'specified in control file: ',tecsav
      write(iout,'(10X,A)') 'Acceptable options are: '
      write(iout,'(10X,A)') OSAT
      write(iout,'(10X,A)') GOAT
      write(iout,'(10X,A)') APCA
      write(iout,'(10X,A)') RTRAC
      write(iout,'(10X,A)') DDM
      write(iout,'(10X,A)') STRPA
      write(iout,'(10X,A)') STRIPR
      write(iout,'(10X,A)') STRIRR
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in STARTUP:'
      write(iout,'(/,1X,2A)') 'Premature end-of-file reached in ',   
     &                  'CAMx.in file when reading input parameters.' 
      write(iout,'(A)') action(:istrln(action))
      call camxerr()   
c
 7005 continue
      write(*,'(//,A)') 'ERROR in STARTUP:'
      write(*,'(2A)') ' ERROR: Could not open file: ',
     &                                   filroot(:istrln(filroot))
      stop
c
 7006 continue
      write(*,'(//,A)') 'ERROR in STARTUP:'
      write(*,'(2A,I2)') ' ERROR: Reading control file ( ',
     &                     ctlfil(:istrln(ctlfil)),' ) at line: ',irec
      write(*,'(A)') action(:istrln(action))
      stop
c
c  --- Return point
c
      end
