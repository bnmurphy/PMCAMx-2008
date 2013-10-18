      subroutine openfils(inp,ii,filroot,ctlfil,irec)
c
c-----CAMx v4.03 031205
c
c     OPENFILS opens Fortran I/O files
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        7/5/02    Changed to account for new type of the PiG flag
c        8/30/02   Modified to read combined cloud/rain file, and now
c                  water vapor and cloud/rain files can be provided for each
c                  nest
c        01/30/02  Added code for RTRAC probing tool
c        1/10/03   Added open of deposition output file
c        4/2/03    Removed option for UAM-V type cloud file
c
c     Input arguments:
c        inp                 CAMx control file unit number
c        ii                  Length of fileroot string 
c        filroot             root name of CAMx output files
c        ctlfil              name of control file
c        irec                line number in control file
c
c     Output arguments:
c        irec                line number in control file
c
c     Routines Called:
c        none
c
c     Called by:
c        STARTUP
c
      include 'camx.prm'
      include 'camx.com'
      include 'camxfld.com'
      include 'filunit.com'
      include 'chmstry.com'
      include 'ahomap.com'
      include 'grid.com'
      include 'flags.com'
c
c======================== Source Apportion Begin =======================
c
      include 'tracer.com'
c
c-----External functions
c
      integer istrln
c
c========================= Source Apportion End ========================
c
c
c======================== Process Analysis Begin =======================
c
      include 'procan.com'
c
c======================== Process Analysis End   =======================
c
      character*200 ctlfil, filroot, filtmp
      character*80  action
      character*15  cldhdr
      logical       lexist
      integer       nopen
c
c-----Entry point
c
c-----Open output files for instantaneous concentration
c
      nopen = 0
      filroot(ii+1:) = '.inst.1'
      iconc(1) = 10 
      nopen = nopen + 1
      action = 'Opening output INST file for coarse grid.'
      open(unit=iconc(1),file=filroot(1:ii+7),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
      write(iout,9000)'1st output INST file coarse grid     (unit):',
     &                                                          iconc(1)
      write(iout,9002) '   File: ',filroot(1:ii+7)
c
      filroot(ii+1:) = '.inst.2'
      iconc(2) = 11 
      nopen = nopen + 1
      action = 'Opening output INST file for coarse grid.'
      open(unit=iconc(2),file=filroot(1:ii+7),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
      write(iout,9000)'2nd output INST file coarse grid     (unit):',
     &                                                          iconc(2)
      write(iout,9002) '   File: ',filroot(1:ii+7)
c
c-----Open optional output concentration files
c
      nfils = 12
      if (nnest.gt.0) then
        filroot(ii+1:) = '.finst.1'
        filtmp = filroot
        ifconc(1) = nfils
        nopen = nopen + 1
        action = 'Opening output INST file for fine grids.'
        open(unit=ifconc(1),file=filroot(1:ii+8),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
        write(iout,9000)'1st output INST file fine grids      (unit):',
     &                                                         ifconc(1)
        write(iout,9002) '   File: ',filroot(1:ii+8)
        nfils = nfils + 1
c
        filroot(ii+1:) = '.finst.2'
        filtmp = filroot
        ifconc(2) = nfils
        nopen = nopen + 1
        action = 'Opening output INST file for fine grids.'
        open(unit=ifconc(2),file=filroot(1:ii+8),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
        write(iout,9000)'2nd output INST file fine grids      (unit):',
     &                                                         ifconc(2)
        write(iout,9002) '   File: ',filroot(1:ii+8)
        nfils = nfils + 1
      endif
c
      if (navspc.gt.0) then
        write(filroot(ii+1:),'(A)') '.avrg'
        filtmp = filroot
        iavg = nfils
        nopen = nopen + 1
        action = 'Opening output AVERAGE file for coarse grid.'
        open(unit=iavg,file=filroot(1:ii+5),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
        write(iout,9000)'Output AVERAGE file coarse grid      (unit):',
     &                                                          iavg
        write(iout,9002) '   File: ',filroot(1:ii+5)
        nfils = nfils + 1
        if (ngrid.gt.1) then
          write(filroot(ii+1:),'(A)') '.favrg'
          filtmp = filroot
          ifavg = nfils
          nopen = nopen + 1
          open(unit=ifavg,file=filroot(1:ii+6),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
          action = 'Opening output AVERAGE file for fine grids.'
         write(iout,9000)'Output AVERAGE file fine grids       (unit):',
     &                                                          ifavg
          write(iout,9002) '   File: ',filroot(1:ii+6)
          nfils = nfils + 1
        endif
      endif
c
      if (ldry .or. lwet) then
        write(filroot(ii+1:),'(A)') '.depn'
        filtmp = filroot
        idep = nfils
        nopen = nopen + 1
        action = 'Opening output DEPOSITION file for coarse grid.'
        open(unit=idep,file=filroot(1:ii+5),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
        write(iout,9000)'Output DEPOSITION file coarse grid   (unit):',
     &                                                          idep
        write(iout,9002) '   File: ',filroot(1:ii+5)
        nfils = nfils + 1
        if (ngrid.gt.1) then
          write(filroot(ii+1:),'(A)') '.fdepn'
          filtmp = filroot
          ifdep = nfils
          nopen = nopen + 1
          open(unit=ifdep,file=filroot(1:ii+6),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
          action = 'Opening output DEPOSITION file for fine grids.'
         write(iout,9000)'Output DEPOSITION file fine grids    (unit):',
     &                                                          ifdep
          write(iout,9002) '   File: ',filroot(1:ii+6)
          nfils = nfils + 1
        endif
      endif
c
      if( ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG ) then
        write(filroot(ii+1:),'(A)') '.pig'
        filtmp = filroot
        ipig = nfils
        nopen = nopen + 1
        action = 'Opening output PiG diagnostic file.'
        open(unit=ipig,file=filroot(1:ii+4),form='UNFORMATTED',
     &                                       status= 'UNKNOWN',ERR=7005)
        write(iout,9000)'Output PiG diagnostics file          (unit):',
     &                                                          ipig
        write(iout,9002) '   File: ',filroot(1:ii+4)
        nfils = nfils + 1
      endif
c
c============================== Process Analysis Begin =============================
c
      if( lipr ) then
         jj = istrln(flrtsa)
         ipr_unit = IOWCN1
         filtmp = flrtsa
         write(filtmp(jj+1:),'(A)') '.ipr'
         open(unit=ipr_unit,file=filtmp(1:jj+4),form='UNFORMATTED',
     &                                    status= 'UNKNOWN',ERR=7000)
         write(iout,9000)'Cell Specific Process output         (unit):',
     &                                                         ipr_unit
         write(iout,9002) '   File: ',filtmp(1:jj+4)
      endif
      if( lirr ) then
         jj = istrln(flrtsa)
         irr_unit = IOWCN2
         filtmp = flrtsa
         write(filtmp(jj+1:),'(A)') '.irr'
         open(unit=irr_unit,file=filtmp(1:jj+4),form='UNFORMATTED',
     &                                     status= 'UNKNOWN',ERR=7000)
         write(iout,9000)'Cell Specific Rates output           (unit):',
     &                                                         irr_unit
         write(iout,9002) '   File: ',filtmp(1:jj+4)
         if( lsfcfl(IDXCRS) ) then
             filtmp = flrtsa
             write(filtmp(jj+1:),'(A)') '.grid.cpa'
             open(unit=IOWSFC+IDXCRS,file=filtmp(1:jj+9),
     &                  form='UNFORMATTED',status= 'UNKNOWN',ERR=7000)
             write(iout,9000)
     &       'Gridded CPA file for coarse grid     (unit):',
     &                                                     IOWSFC+IDXCRS
             write(iout,9002) '   File: ',filtmp(1:jj+9)
             sfcfil(IDXCRS) = filtmp
         endif
         if( lsfcfl(IDXFIN) ) then
             filtmp = flrtsa
             write(filtmp(jj+1:),'(A)') '.fgrid.cpa'
             open(unit=IOWSFC+IDXFIN,file=filtmp(1:jj+10),
     &                  form='UNFORMATTED',status= 'UNKNOWN',ERR=7000)
             write(iout,9000)
     &       'Gridded CPA file for fine grids      (unit):',
     &                                                     IOWSFC+IDXFIN
             write(iout,9002) '   File: ',filtmp(1:jj+10)
             sfcfil(IDXFIN) = filtmp
         endif
      endif
c
c=============================== Process Analysis End ==============================
c
c
c-----Read and open input file names
c
      filtmp=' '
      irec = irec + 1
      action = 'Reading CHEMPARAM filename.'
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp
      ichem = nfils
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      open(unit=ichem,file=filtmp,status='OLD',ERR=7000)
      write(iout,9000)'Chemistry parameters file            (unit):',
     &                                                          ichem
      write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
      action = 'Reading Photolysis rates filename.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( lchem ) then
        iphot = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        open(unit=iphot,file=filtmp,status='OLD',ERR=7000) 
        write(iout,9000)'Photolysis rates file                (unit):',
     &                                                           iphot
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        write(iout,9000)'Photolysis rates file                      :'
        write(iout,9002) '   Record ignored.' 
      endif
c
      action = 'Reading Landuse filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( ldry ) then
        isurf(1) = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        open(unit=isurf(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
        write(iout,9000)'Landuse file                         (unit):',
     &                                                         isurf(1)
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        write(iout,9000)'Landuse file                               :'
        write(iout,9002) '   Record ignored.' 
      endif
c
      action = 'Reading Height/pressure filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      ihtp(1) = nfils
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      open(unit=ihtp(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      write(iout,9000)'Height/pressure file                 (unit):',
     &                                                         ihtp(1)
      write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
      action = 'Reading Wind filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      iwind(1) = nfils
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      open(unit=iwind(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      write(iout,9000)'Wind file                            (unit):',
     &                                                        iwind(1)
      write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
      action = 'Reading Temperature filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      itemp(1) = nfils
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist ) goto 7002
      nopen = nopen + 1
      open(unit=itemp(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      write(iout,9000)'Temperature file                     (unit):',
     &                                                        itemp(1)
      write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
      action = 'Reading Wator vapor filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( lchem ) then
        ih2o(1) = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        open(unit=ih2o(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
        write(iout,9000)'Water vapor file                     (unit):',
     &                                                         ih2o(1)
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        write(iout,9000)'Water vapor file                           :'
        write(iout,9002) '   Record ignored.' 
      endif
c
      action = 'Reading Cloud/rain filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( filtmp .EQ. ' ' ) then
        write(iout,9000)'Cloud/rain file                            :'
        write(iout,9002) '   Record ignored.' 
      endif
      icld(1) = 0
      if( (lchem .or. lwet .or. ldry) .AND. filtmp .NE. ' ')  then
        icld(1) = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist ) goto 7002
        nopen = nopen + 1
        open(unit=icld(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
        write(iout,9000)'Cloud/rain file                      (unit):',
     &                                                         icld(1)
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1

        read(icld(1),ERR=7008) cldhdr
        call jstlft( cldhdr )
        call toupper( cldhdr )
        if( cldhdr .EQ. 'CAMX CLOUD_RAIN' ) then
          backspace(icld(1))
          read(icld(1),ERR=7008) cldhdr,nxcl,nycl,nzcl
          if (nxcl.ne.ncol(1) .or. nycl.ne.nrow(1) .or.
     &                                          nzcl.ne.nlay(1)) then
            write(iout,'(//,a)') 'ERROR in OPENFILS:'
            write(iout,'(2A)')'Cloud/rain file dimensions do not',
     &                        ' match coarse grid.'
            write(iout,*)'Cloud file: ',nxcl,nycl,nzcl
            write(iout,*)'CAMx      : ',ncol(1),nrow(1),nlay(1)
            call camxerr()
          endif
        else
          if( cldhdr(1:10) .EQ. 'CAMX CLOUD' ) goto 7011
          goto 7006
        endif
      endif
      if( lwet .and. filtmp .EQ. ' ' ) goto 7003
c
      action = 'Reading Vertical diffusivity filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      ikv(1) = nfils
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7009
      nopen = nopen + 1
      open(unit=ikv(1),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      write(iout,9000)'Vertical diffusivity file            (unit):',
     &                                                          ikv(1)
      write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
      action = 'Reading Initial conditions filename.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( .NOT. lrstrt ) then
        iic = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
        nopen = nopen + 1
        open(unit=iic,file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
        write(iout,9000)'Initial conditions file              (unit):',
     &                                                            iic
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        write(iout,9000)'Initial conditions file                    :'
        write(iout,9002) '   Record ignored.' 
      endif
c
      action = 'Reading Boundary conditions filename.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      ibc = nfils
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
      nopen = nopen + 1
      open(unit=ibc,file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
      write(iout,9000)'Boundary conditions file             (unit):',
     &                                                            ibc
      write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
      action = 'Reading Top concentrations filename.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      itopc = nfils
      inquire(file=filtmp,exist=lexist)
      if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
      nopen = nopen + 1
      open(unit=itopc,file=filtmp,status='OLD',ERR=7000) 
      write(iout,9000)'Top concentration file               (unit):',
     &                                           itopc
      write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
      nfils = nfils + 1
c
      action = 'Reading Albedo/haze/ozone filename.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( lchem ) then
        iaho = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
        nopen = nopen + 1
        open(unit=iaho,file=filtmp,status='OLD',ERR=7000) 
        write(iout,9000)'Albedo/haze/ozone file               (unit):',
     &                                                            iaho
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        write(iout,9000)'Albedo/haze/ozone file                     :'
        write(iout,9002) '   Record ignored.' 
      endif
c
      action = 'Reading Point source emissions filename.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( lptsrc ) then
        iptem = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
        nopen = nopen + 1
        open(unit=iptem,file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
        write(iout,9000)'Point source emissions file          (unit):',
     &                                                            iptem
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        write(iout,9000)'Point source emissions file                :'
        write(iout,9002) '   Record ignored.' 
      endif
c
      action = 'Reading Area source emissions filename for coarse grid.'
      irec = irec + 1
      read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
      if( larsrc ) then
        iarem(1,1) = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
        nopen = nopen + 1
        open(unit=iarem(1,1),file=filtmp,form='UNFORMATTED',
     &                                         status='OLD',ERR=7000)
        write(iout,9000)'Area source emissions file           (unit):',
     &                                                      iarem(1,1)
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
      else
        write(iout,9000)'Area source emissions file                 :'
        write(iout,9002) '   Record ignored.' 
      endif
c
c-----If nested grids are to be simulated, read and open fine grid files
c
      do 25 n = 2,ngrid
       isurf(n) = 0
       write(action,'(A,I4)') 'Reading Landuse filename for nest:',n-1
       irec = irec + 1
       read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
       if( ldry .AND. filtmp .NE. ' ')  then
         isurf(n) = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
         nopen = nopen + 1
         open(unit=isurf(n),file=filtmp,form='UNFORMATTED',
     &                                           status='OLD',ERR=7000)
         write(iout,9001)'Fine landuse file for grid # ',n,
     &                   '      (unit):',nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1
       else
         write(iout,9001)'Fine landuse file for grid # ',n,
     &                                                  '            :'
         write(iout,9002) '   Record ignored.' 
       endif
 25   continue
c
      do 30 n = 2,ngrid
        ihtp(n) = 0
        write(action,'(A,I4)') 
     &               'Reading Height/pressure filename for nest:',n-1
        irec = irec + 1
        read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
        if( filtmp .NE. ' ' ) then
          ihtp(n) = nfils
          inquire(file=filtmp,exist=lexist)
          if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
          nopen = nopen + 1
          open(unit=ihtp(n),file=filtmp,form='UNFORMATTED',
     &                                          status='OLD',ERR=7000)
          write(iout,9001)'Fine height/press file for grid # ',n,
     &                   ' (unit):',nfils
          write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
          nfils = nfils + 1
        else
          write(iout,9001)'Fine height/press file for grid # ',n,
     &                                                   '       :'
          write(iout,9002) '   Record ignored.' 
        endif
 30   continue
c
      do 35 n = 2,ngrid
       iwind(n) = 0
       write(action,'(A,I4)') 'Reading wind filename for nest:',n-1
       irec = irec + 1
       read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
       if( filtmp .NE. ' ' ) then
         iwind(n) = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
         nopen = nopen + 1
         open(unit=iwind(n),file=filtmp,form='UNFORMATTED',
     &                                        status='OLD',ERR=7000)
         write(iout,9001)'Fine wind file for grid # ',n,
     &                   '         (unit):',nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1
       else
         write(iout,9001)'Fine wind file for grid # ',n,
     &                                               '               :'
         write(iout,9002) '   Record ignored.' 
       endif
 35   continue
c
      do 36 n = 2,ngrid
       itemp(n) = 0
       write(action,'(A,I4)') 
     &               'Reading temperature filename for nest:',n-1
       irec = irec + 1
       read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
       if( filtmp .NE. ' ' ) then
         itemp(n) = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
         nopen = nopen + 1
         open(unit=itemp(n),file=filtmp,form='UNFORMATTED',
     &                                        status='OLD',ERR=7000)
         write(iout,9001)'Fine temp file for grid # ',n,
     &                   '         (unit):',nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1
       else
         write(iout,9001)'Fine temp file for grid # ',n,
     &                                               '               :'
         write(iout,9002) '   Record ignored.' 
       endif
 36   continue
c
      do 37 n = 2,ngrid
       ih2o(n) = 0
       write(action,'(A,I4)') 
     &               'Reading water vapor filename for nest:',n-1
       irec = irec + 1
       read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
       if( lchem .and. filtmp .NE. ' ' ) then
         ih2o(n) = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist ) goto 7002
         nopen = nopen + 1
         open(unit=ih2o(n),file=filtmp,form='UNFORMATTED',
     &                                        status='OLD',ERR=7000)
         write(iout,9001)'Fine water vapor file for grid # ',n,
     &                   '         (unit):',nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1
       else
         write(iout,9001)'Fine water vapor file for grid # ',n,
     &                                                  '        :'
         write(iout,9002) '   Record ignored.' 
       endif
 37   continue
c
      do 38 n = 2,ngrid
       icld(n) = 0
       write(action,'(A,I4)') 
     &               'Reading cloud/rain filename for nest:',n-1
       irec = irec + 1
       read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
       if( (lchem .or. lwet .or. ldry) .and. filtmp .NE. ' ' ) then
         icld(n) = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist ) goto 7002
         nopen = nopen + 1
         open(unit=icld(n),file=filtmp,form='UNFORMATTED',
     &                                        status='OLD',ERR=7000)
         write(iout,9001)'Fine cloud file for grid # ',n,
     &                   '         (unit):',nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1

         read(icld(n)) cldhdr
         call jstlft( cldhdr )
         call toupper( cldhdr )
         if( cldhdr .EQ. 'E   M   I   S  ' ) goto 7010
         if( cldhdr .EQ. 'CAMX CLOUD_RAIN' ) then
           backspace(icld(n))
           read(icld(n)) cldhdr,nxcl,nycl,nzcl
           if (nxcl.ne.ncol(n) .or. nycl.ne.nrow(n) .or.      
     &                                         nzcl.ne.nlay(n)) then
             write(iout,'(//,a)') 'ERROR in OPENFILS:'
             write(iout,'(2A)')'Cloud file dimensions do not',
     &                         ' match nested grid.'
             write(iout,*)'Cloud file: ',nxcl,nycl,nzcl
             write(iout,*)'CAMx      : ',ncol(n),nrow(n),nlay(n)
             call camxerr()   
           endif
         else
           if( cldhdr(1:10) .EQ. 'CAMX CLOUD') goto 7011
           goto 7006 
         endif
       else
         write(iout,9001)'Fine cloud file for grid # ',n,
     &                                                '              :'
         write(iout,9002) '   Record ignored.' 
       endif
 38   continue
c
      do 40 n = 2,ngrid
       ikv(n) = 0
       write(action,'(A,I4)') 
     &     'Reading Vertical diffusivity filename for nest:',n-1
       irec = irec + 1
       read(inp,'(20X,A)',ERR=7001,END=7012) filtmp 
       if( filtmp .NE. ' ' ) then
         ikv(n) = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
         nopen = nopen + 1
         open(unit=ikv(n),file=filtmp,form='UNFORMATTED',
     &                                      status='OLD',ERR=7000)
         write(iout,9001)'Fine Kv file for grid # ',n,
     &                   '           (unit):',nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1
       else
         write(iout,9001)'Fine Kv file for grid # ',n,
     &                                          '                 :'
         write(iout,9002) '   Record ignored.' 
       endif
 40   continue
c
      do 45 n = 2,ngrid
       iarem(1,n) = 0
       write(action,'(A,I4)') 
     &     'Reading Area source emissions filename for nest:',n-1
       irec = irec + 1
       read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
       if( .NOT. larsrc ) goto 45
       if( filtmp .NE. ' ' ) then
         iarem(1,n) = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
         nopen = nopen + 1
         open(unit=iarem(1,n),file=filtmp,form='UNFORMATTED',
     &                                    status='OLD',ERR=7000)
         write(iout,9001)'Fine area emiss file for grid # ',n,
     &                   '   (unit):',nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1
       else
         write(iout,9001)'Fine area emiss file for grid # ',n,
     &                                                  '         :'
         write(iout,9002) '   Record ignored.' 
       endif
 45   continue
c
c-----Open restart files if LRSTRT = T
c
      if( lrstrt ) then
        action = 'Reading restart filename for coarse grid.'
        irec = irec + 1
        read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
        irstc = nfils
        inquire(file=filtmp,exist=lexist)
        if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
        nopen = nopen + 1
        open(unit=irstc,file=filtmp,form='UNFORMATTED',
     &                                    status='OLD',ERR=7005)
        write(iout,9000)'Coarse grid restart file             (unit):',
     &                                                         nfils
        write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
        nfils = nfils + 1
c
        if( nnest .GT. 0 ) then
          action = 'Reading restart filename for fine grids.'
          irec = irec + 1
          read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
          if( filtmp .EQ. ' ' ) then 
             irstf = 0
             write(iout,9001)
     &               'Fine grid restart file                     :'
             write(iout,9002) '   Record ignored.' 
          else 
             irstf = nfils
             inquire(file=filtmp,exist=lexist)
             if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
             nopen = nopen + 1
             open(unit=irstf,file=filtmp,form='UNFORMATTED',
     &                                  status='OLD',ERR=7005)
             write(iout,9000)
     &          'Fine grid restart file               (unit):',nfils
             write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
             nfils = nfils + 1
          endif
        endif
c
        if( ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG ) then
         action = 'Reading PiG restart filename.'
         irec = irec + 1
         read(inp,'(20X,A)',ERR=7001,END=7004) filtmp 
         irstp = nfils
         inquire(file=filtmp,exist=lexist)
         if( .NOT. lexist .OR. filtmp .EQ. ' ' ) goto 7002
         nopen = nopen + 1
         open(unit=irstp,file=filtmp,form='UNFORMATTED',
     &                                  status='OLD',ERR=7005)
         write(iout,9000)'PiG restart file                     (unit):',
     &                                                            nfils
         write(iout,9002) '   File: ',filtmp(:istrln(filtmp))
         nfils = nfils + 1
       endif
      endif
c
c======================== Source Apportion Begin =======================
c
c  --- call routine to get the filenames for source apportionment data ---
c
      if( ltrace ) then
          if( tectyp .EQ. RTRAC ) then
              call startrt(inp,iout,nopen)
          else
              call startsa(inp,iout,nopen)
          endif
      endif
      if( lddm ) then
          call startddm(inp,iout,nopen)
      endif
c
c========================= Source Apportion End ========================
c
      write(iout,*)
      write(iout,*)'Finished reading control file.'
      write(iout,'(A,I2)')'Number of input files opened: ',nopen
      write(iout,*)
      goto 9999
c
c  --- Format statements ---
c
 9000 format(/,A,I2)
 9001 format(/,A,I2,A,I2)
 9002 format(2A)
c
c  --- Error messages
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Could not open file: ',
     &                                   filtmp(:istrln(filtmp))
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(3A,I2)') 'Reading control file ( ',
     &                     ctlfil(:istrln(ctlfil)),' ) at line: ',irec
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      if( filtmp .EQ. ' ' ) then
         write(iout,'(2A,I2)') 'Blank filename provided in ',
     &                         'control file at line: ',irec
      else
         write(iout,'(2A)') 'Input file does not exist: ',
     &                                       filtmp(:istrln(filtmp))
      endif
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(2A)') 'A cloud/rain file must be supplied if the ',
     &                   'wet deposition flag is set.'
      write(iout,'(2A,/,A)') 'Either turn off the wet dep flag or ',
     &            'supply a cloud/rain file in the','CAMx control file.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(/,1X,2A)') 'Premature end-of-file reached in ',
     &                  'CAMx.in file when reading filenames.'
      write(iout,'(1X,2A)') 'Please see the .out file for a history ',
     &             'of filenames read from the control file.'
      write(iout,'(1X,2A,/,1X,A)') 'If there is a mismatch, refer to ',
     &                     'the CAMx user guide for the proper order ',
     &                                       'of options and filenames.'
      call camxerr()
c
 7005 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'Could not open file: ',
     &                                   filtmp(:istrln(filtmp))
      write(iout,'(10X,2A)') 'Make sure the name for output files ',
     &                       'reflects the correct simulation day.'
      call camxerr()
c
 7006 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') 'CAMx cloud/rain file header is incorrect.'
      write(iout,'(2A)') 'This version of CAMx incorporates a ',
     &                     'different format for the Cloud/Rain File.'
      write(iout,'(A)') 'Do not include a line for the RAIN file.'
      write(iout,'(A)') 'Expecting: CAMx CLOUD/RAIN'
      write(iout,'(2A)')'Read: ',cldhdr
      call camxerr()
c
 7008 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)')'Could not read first record of cloud/rain file.'
      write(iout,'(A)') filtmp(:istrln(filtmp))
      call camxerr()
c
 7009 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A,I2)') 'Blank filename provided in ',
     &                         'control file at line: ',irec
      write(iout,'(A)') 
     &     'This version of CAMx uses a combined CLOUD and RAIN file.' 
      write(iout,'(A)') 'Do not include a line for the RAIN file.'
      call camxerr()
c
 7010 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'This version of CAMx requires a Cloud/Rain ',
     &                           'file and a Water Vapor file for each'
      write(iout,'(2A)') 'nest.  If you do not have a file, just ',
     &                    'leave a blank line for each nest after the '
      write(iout,'(2A)') 'temperature files but before the KV files.'
      call camxerr()
c
 7011 continue
      write(iout,'(//,a)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(A)') 
     &     'This version of CAMx uses a combined CLOUD and RAIN file.' 
      write(iout,'(2A)') 'It appears you are using a CLOUD ',
     &                   'file from a previous version of CAMx.'
      call camxerr()
c
 7012 continue
      write(iout,'(//,A)') 'ERROR in OPENFILS:'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(/,1X,2A)') 'Premature end-of-file reached in ',
     &                  'CAMx.in file when reading filenames.'
      write(iout,'(1X,2A)') 'This version of CAMx requires that you ',
     &                               'provide a water vapor file',
     &              'and a cloud/rain file for each nest.'
      write(iout,'(1X,2A)') 'If the filename is blank the ',
     &           'fields will be interpolated from the coarse grid.'
      write(iout,'(1X,2A)') 'Please see the .out file for a history ',
     &                    'of filenames read from the control file.'
      write(iout,'(1X,3A)') 'Please refer to ',
     &                     'the CAMx user guide for the proper order ',
     &                                       'of options and filenames.'
      call camxerr()
c
 9999 continue
      return
      end
