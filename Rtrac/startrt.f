c**** STARTRT
c
      subroutine startrt(iounit,iout,nopen)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c      This routine reads all of the files needed for the reactive 
c      tracer algorithm (RTRAC)
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c      Argument description:
c        Inputs:
c          iounit   I   unit number for input
c          iout     I   unit number for output
c        Outputs:   
c          nopen    I   number of files opened
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/16/02   --gwilson--   Original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'grid.com'
      include 'flags.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Argument declaration:
c-----------------------------------------------------------------------
c
      integer   iounit
      integer   iout
      integer   nopen
c
c-----------------------------------------------------------------------
c    External functions:
c-----------------------------------------------------------------------
c
      integer istrln
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*80 action
      integer      i, j
      logical      lexist
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. ltrace .OR. tectyp .NE. RTRAC ) goto 9999
      nchar = istrln( flrtsa )
c
c   --- intialize first time through to true ---
c
      lfirst = .TRUE.
      write(iout,'(/,A,/)') 
     &            '           **** Source Apportionment Files ****'
c
c   --- read the filenames from the job script files (standard input)
c       and open the files ---
c
c   --- chemistry parameters file ---
c
      action = 'Reading filename for chemistry definition file for '//
     &                                             'RTRAC.'
      read(iounit,'(20X,A)',END=7003) chmfil
      fname = chmfil
      inquire(file=chmfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      nopen =  nopen + 1
      write(iout,9000)'Chemistry definition file             (unit):',
     &                                                            IORCHM
      write(iout,9002) '   File: ',chmfil(:istrln(rcpfil))
c
c   --- receptor definition file ----
c
      action = 'Reading filename for receptor definition for '//
     &                                             'RTRAC.'
      read(iounit,'(20X,A)',END=7003) rcpfil
      fname = rcpfil
      inquire(file=rcpfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      nopen =  nopen + 1
      write(iout,9000)'Receptor definition file             (unit):',
     &                                                            IORRCP
      write(iout,9002) '   File: ',rcpfil(:istrln(rcpfil))
c
c   --- instantaneous file used for initialization ---
c
      if( lrstrt ) then
         action = 'Reading coarse grid reatsrt filename for '//
     &                                             'RTRAC.'
         read(iounit,'(20X,A)',END=7003) inifil(IDXCRS)
         fname = inifil(IDXCRS)
         inquire(file=inifil(IDXCRS),exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         write(iout,9000)
     &         'Instantaneous file for coarse grid   (unit):',IORINI
         write(iout,9002) '   File: ',
     &                          inifil(IDXCRS)(:istrln(inifil(IDXCRS)))
         if( ngrid .GT. 1 ) then
            action = 'Reading fine grid restart filename for '//
     &                                             'RTRAC.'
            read(iounit,'(20X,A)',END=7003) inifil(IDXFIN)
            if( inifil(IDXFIN) .NE. ' ' ) then
                write(iout,9002) '   File: ',
     &                          inifil(IDXFIN)(:istrln(inifil(IDXFIN)))
                
                fname = inifil(IDXFIN)
                inquire(file=inifil(IDXFIN),exist=lexist)
                if( .NOT. lexist ) goto 7000
                nopen =  nopen + 1
                write(iout,9001)
     &             'Instantaneous file for fine grids    (unit):',
     &                                                        IORINI
                write(iout,9002) '   File: ',
     &                          inifil(IDXFIN)(:istrln(inifil(IDXFIN)))
            endif
         endif
      else
c
c   --- initial conditions file ---	
c
         action = 'Reading initial conditions filename for '//
     &                                             'RTRAC.'
         read(iounit,'(20X,A)',END=7003) icfil
         fname = icfil
         inquire(file=icfil,exist=lexist)
         if( .NOT. lexist ) goto 7000
         nopen =  nopen + 1
         write(iout,9000)
     &             'Initial conditions file for RTRAC    (unit):',IORIC
         write(iout,9002) '   File: ',icfil(:istrln(icfil))
      endif
c
c   --- boundary conditions file ---
c
      action = 'Reading boundary conditions filename for '//
     &                                             'RTRAC.'
      read(iounit,'(20X,A)',END=7003) bcfil
      fname = bcfil
      inquire(file=bcfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      nopen =  nopen + 1
      write(iout,9000)
     &           'Boundary conditions file for RTRAC    (unit):',IORBC
      write(iout,9002) '   File: ',bcfil(:istrln(bcfil))
c
c   
c   --- boundary conditions file ---
c
      action = 'Reading top concentrations filename for '//
     &                                             'RTRAC.'
      read(iounit,'(20X,A)',END=7003) tcfil
      fname = tcfil
      inquire(file=tcfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      nopen =  nopen + 1
      write(iout,9000)
     &           'Top concentrations file for RTRAC     (unit):',IORTC
      write(iout,9002) '   File: ',tcfil(:istrln(bcfil))
c
c   --- only one group for RTRAC ----
c
      ngroup = 1
c
c   --- emissions files for toxics ---
c
c   --- surface emissions file --- 
c
      do 10 j=1,ngrid
          action = 'Reading area emissions filename for '//
     &                                             'RTRAC.'
          read(iounit,'(20X,A)',END=7003) temfil(j,1)
          fname = temfil(j,1)
          if( fname .EQ. ' ' ) then
             goto 7004
          else
             ltemfl(j,1) = .TRUE.
             inquire(file=temfil(j,1),exist=lexist)
             if( .NOT. lexist ) goto 7000
             nopen = nopen + 1
             open(unit=IORTEM+j,file=temfil(j,1),
     &                   ERR=7001,form='UNFORMATTED',status='UNKNOWN')
             if( j .GT. 0 ) then
                write(iout,9005)
     &             'Surface emissions file grid#         (unit):',
     &                                                      j,IORTEM+j
                write(iout,9002) '   File: ',
     &                                 temfil(j,1)(:istrln(temfil(j,1)))
             else
                write(iout,9005)
     &             'Surface emissions file coarse grid   (unit):',
     &                                                        IORTEM+j
                write(iout,9002) '   File: ',
     &                                 temfil(j,1)(:istrln(temfil(j,1)))
             endif
          endif
   10 continue
c
c   --- elevated point source emissions ----
c
      read(iounit,'(20X,A)',END=7003) tptfil(1)
      fname = tptfil(1)
      action = 'Reading point emissions filename for '//
     &                                             'RTRAC.'
      if( fname .EQ. ' ' ) then
         goto 7005
      else
          ltptfl(1) = .TRUE.
          inquire(file=tptfil(1),exist=lexist)
          if( .NOT. lexist ) goto 7000
          nopen = nopen + 1
          open(unit=IORTPT,file=tptfil(1),
     &                    ERR=7001,form='UNFORMATTED',status='UNKNOWN')
          write(iout,9005)
     &                 'Elevated emissions file for RTRAC    (unit):',
     &                                                          IORTPT
          write(iout,9002) '   File: ',tptfil(1)(:istrln(tptfil(1)))
      endif
c
c   --- automatically output average concs for RTRAC ---
c
      lsfcfl(IDXCRS) = .TRUE.
      lsfcfl(IDXFIN) = .TRUE.
c
c   --- output filenames ---
c
      do 20 i=IDXCRS,IDXFIN
         if( i .EQ. IDXFIN .AND. ngrid .EQ. 1 ) goto 20
c
c   --- first output instantaneous file ---
c
         cn1fil(i) = flrtsa 
         if( i .EQ. IDXCRS ) then
            cn1fil(i)(nchar+1:) = '.rt.inst.1'
         else
            cn1fil(i)(nchar+1:) = '.rt.finst.1'
         endif
         fname = cn1fil(i)
         nopen = nopen + 1
         open(unit=IOWCN1+i,file=cn1fil(i),ERR=7002,
     &                          form='UNFORMATTED',status='UNKNOWN')
         if( i .EQ. IDXCRS ) then
            write(iout,9000)
     &             'First output inst. file coarse grid  (unit):',
     &                                                          IOWCN1+i
            write(iout,9002) '   File: ',cn1fil(i)(:istrln(cn1fil(i)))
         else
            write(iout,9000)
     &             'First output inst. file fine grids   (unit):',
     &                                                          IOWCN1+i
            write(iout,9002) '   File: ',cn1fil(i)(:istrln(cn1fil(i)))
         endif
c
c    ---- second output instantaneous file ---
c
         cn2fil(i)(1:) = flrtsa(1:nchar) 
         if( i .EQ. IDXCRS ) then
            cn2fil(i)(nchar+1:) = '.rt.inst.2'
         else
            cn2fil(i)(nchar+1:) = '.rt.finst.2'
         endif
         fname = cn2fil(i)
         nopen = nopen + 1
         open(unit=IOWCN2+i,file=cn2fil(i),ERR=7002,
     &                              form='UNFORMATTED',status='UNKNOWN')
         if( i .EQ. IDXCRS ) then
            write(iout,9000)
     &             'Second output inst. file coarse grid (unit):',
     &                                                          IOWCN2+i
            write(iout,9002) '   File: ',cn2fil(i)(:istrln(cn2fil(i)))
         else
            write(iout,9000)
     &             'Second output inst. file fine grids  (unit):',
     &                                                          IOWCN2+i
            write(iout,9002) '   File: ',cn2fil(i)(:istrln(cn2fil(i)))
         endif
c
c    ---- surface concentrations file ----
c
         if( .NOT. lsfcfl(i) ) then
             if( i .EQ. IDXCRS ) then
                 write(iout,9006)
     &             'Output surface conc. file coarse grid(unit):',
     &                                                  ' Not supplied'
             else
                 write(iout,9006)
     &             'Output surface conc. file fine grids (unit):',
     &                                                  ' Not supplied'
             endif
          else
             sfcfil(i)(1:) = flrtsa(1:nchar) 
             if( i .EQ. IDXCRS ) then
                sfcfil(i)(nchar+1:) = '.rt.surf'
             else
                sfcfil(i)(nchar+1:) = '.rt.fsurf'
             endif
             fname = sfcfil(i)
             nopen = nopen + 1
             open(unit=IOWSFC+i,file=sfcfil(i),ERR=7002,
     &                             form='UNFORMATTED',status='UNKNOWN')
             if( i .EQ. IDXCRS ) then
               write(iout,9000)
     &             'Output surface conc. file coarse grid(unit):',
     &                                                          IOWSFC+i
               write(iout,9002)'   File: ',sfcfil(i)(:istrln(sfcfil(i)))
             else
               write(iout,9000)
     &             'Output surface conc. file fine grids (unit):',
     &                                                          IOWSFC+i
               write(iout,9002)'   File: ',sfcfil(i)(:istrln(sfcfil(i)))
             endif
          endif
   20 continue
c
c    ---- receptor output file ----
c
      avgfil(1:) = flrtsa(1:nchar) 
      avgfil(nchar+1:) = '.rt.receptor'
      fname = avgfil
      nopen = nopen + 1
      open(unit=IOWAVG,file=avgfil,ERR=7002,
     &                             form='FORMATTED',status='UNKNOWN')
      write(iout,9000)
     &             'Output receptor decay rate file   (unit):',IOWAVG
      write(iout,9002)'   File: ',avgfil(:istrln(avgfil))
c
c  --- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
 9000 format(/,A,I2)
 9001 format(/,A,I2,1X,I2)
 9002 format(2A)
 9005 format(/,A,I2,1X,I2,1X,I2)
 9006 format(/,2A)
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in STARTRT:'
      write(iout,'(/,1X,2A)') 'Input file does not exist: ',
     &                                           fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in STARTRT:'
      write(iout,'(/,1X,2A)') 'Cannot open input file: ',
     &                                           fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in STARTRT:'
      write(iout,'(/,1X,2A)') 'Cannot open output file: ',
     &                                           fname(:istrln(fname))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in STARTRT:'
      write(iout,'(/,1X,2A)') 'Premature end-of-file reached in ',
     &                  'CAMx.in file when reading RTRAC filenames.'
      write(iout,'(A)') action(:istrln(action))
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in STARTRT:'
      write(iout,'(/,1X,2A,I5)') 'A filename must be supplied ',
     &                                             'for grid: ',j
      write(iout,'(A)') action(:istrln(action))
      call camxerr()
c
 7005 continue
      write(iout,'(//,A)') 'ERROR in STARTRT:'
      write(iout,'(/,1X,A)') 'A filename must be supplied.'
      write(iout,'(A)') action(:istrln(action))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
