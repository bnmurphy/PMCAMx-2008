c**** STARTSA
c
      subroutine startsa(iounit,iout,nopen)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c      This routine reads all of the files needed for the source 
c      apportionment algorithm.  The filenames are read and the output
c      files are opened.  The input files are opened as needed.
c
c      Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
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
c     12/16/99   --gwilson--   Fixed a bug which caused the model to try
c                              and read a restart file for each find grid
c     07/19/02   --gwilson--   Added code for source area map for nests
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
      integer      igrd, i, j
      logical      lexist
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- if option is turned off, just return ---
c
      if( .NOT. ltrace ) goto 9999
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
c       source area mapping file ---
c
      if( nregin .GT. 0 ) then
         do 10 igrd=1,ngrid
            write(action,'(2A,I2)') 'Reading filename source area map ',
     &                          'for source apportionment.  Grid: ',igrd
            read(iounit,'(20X,A)',END=7003) mapfil(igrd)
            fname = mapfil(igrd)
            if( fname .EQ. ' ' ) then
               if( igrd .EQ. 1 ) goto 7004
               lmapfl(igrd) = .FALSE.
               write(iout,9001)
     &            'No Source area mapping file for Grid #        :',igrd
               goto 10
            endif
            inquire(file=mapfil(igrd),exist=lexist)
            if( .NOT. lexist ) goto 7000
            nopen =  nopen + 1
            write(iout,9001)
     &        'Source area mapping file for Grid #  (unit):',igrd,IORMAP
            write(iout,9002) '   File: ',
     &                               mapfil(igrd)(:istrln(mapfil(igrd)))
             lmapfl(igrd) = .TRUE.
   10    continue
      endif
c
c   --- receptor definition file ----
c
      action = 'Reading filename for receptor definition for source '//
     &                                             'apportionment.'
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
         action = 'Reading coarse grid reatsrt filename for source '//
     &                                             'apportionment.'
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
     &                                          'source apportionment.'
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
      endif
c
c   --- emissions files for the source groupings ---
c
      do 20 i=1,ngroup
c
c   --- surface emissions file --- 
c
         do 30 j=1,ngrid
             action = 'Reading area emissions filename for '//
     &                                          'source apportionment.'
             read(iounit,'(20X,A)',END=7003) temfil(j,i)
             fname = temfil(j,i)
             if( fname .EQ. ' ' ) then
                ltemfl(j,i) = .FALSE.
                if( j .GT. 0 ) then
                    write(iout,9003)
     &                  'Surface emissions file for grid#/group#    :',
     &                                            j,i,' Not supplied'
                else
                    write(iout,9004)
     &                  'Surface emissions file for coarse/group#   :',
     &                                                i,' Not supplied'
                endif
             else
                 ltemfl(j,i) = .TRUE.
                 inquire(file=temfil(j,i),exist=lexist)
                 if( .NOT. lexist ) goto 7000
                 nopen = nopen + 1
                 open(unit=IORTEM+(j-1)*ngroup+i,file=temfil(j,i),
     &                   ERR=7001,form='UNFORMATTED',status='UNKNOWN')
                 if( j .GT. 0 ) then
                    write(iout,9005)
     &                 'Surface emissions file grid#/group#  (unit):',
     &                                     j,i,IORTEM+(j-1)*ngroup+i
                    write(iout,9002) '   File: ',
     &                                 temfil(j,i)(:istrln(temfil(j,i)))
                 else
                    write(iout,9005)
     &                 'Surface emissions file coarse/group# (unit):',
     &                                         i,IORTEM+(j-1)*ngroup+i
                    write(iout,9002) '   File: ',
     &                                 temfil(j,i)(:istrln(temfil(j,i)))
                 endif
             endif
   30    continue
c
c   --- elevated point source emissions ----
c
         read(iounit,'(20X,A)',END=7003) tptfil(i)
         fname = tptfil(i)
         action = 'Reading point emissions filename for '//
     &                                          'source apportionment.'
         if( fname .EQ. ' ' ) then
             ltptfl(i) = .FALSE.
             write(iout,9004)
     &                  'Elevated emissions file for group#         :',
     &                                                i,' Not supplied'
         else
             ltptfl(i) = .TRUE.
             inquire(file=tptfil(i),exist=lexist)
             if( .NOT. lexist ) goto 7000
             nopen = nopen + 1
             open(unit=IORTPT+i,file=tptfil(i),
     &                    ERR=7001,form='UNFORMATTED',status='UNKNOWN')
             write(iout,9005)
     &                 'Elevated emissions file for group#   (unit):',
     &                                                        i,IORTPT+i
             write(iout,9002) '   File: ',tptfil(i)(:istrln(tptfil(i)))
         endif
   20 continue
c
c   --- output filenames ---
c
      do 40 i=IDXCRS,IDXFIN
         if( i .EQ. IDXFIN .AND. ngrid .EQ. 1 ) goto 40
c
c   --- first output instantaneous file ---
c
         cn1fil(i) = flrtsa 
         if( i .EQ. IDXCRS ) then
            cn1fil(i)(nchar+1:) = '.sa.inst.1'
         else
            cn1fil(i)(nchar+1:) = '.sa.finst.1'
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
            cn2fil(i)(nchar+1:) = '.sa.inst.2'
         else
            cn2fil(i)(nchar+1:) = '.sa.finst.2'
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
         if( .not.lsfcfl(i) ) then
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
                sfcfil(i)(nchar+1:) = '.sa.surf'
             else
                sfcfil(i)(nchar+1:) = '.sa.fsurf'
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
   40 continue
c
c    ---- tracer average file ----
c
      avgfil(1:) = flrtsa(1:nchar) 
      avgfil(nchar+1:) = '.sa.receptor'
      fname = avgfil
      nopen = nopen + 1
      open(unit=IOWAVG,file=avgfil,ERR=7002,
     &                             form='FORMATTED',status='UNKNOWN')
      write(iout,9000)
     &             'Output average receptor conc. file   (unit):',IOWAVG
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
 9003 format(/,A,1X,I2,1X,I2,A)
 9004 format(/,A,1X,I2,A)
 9005 format(/,A,I2,1X,I2,1X,I2)
 9006 format(/,2A)
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in STARTSA:'
      write(iout,'(/,1X,2A)') 'Input file does not exist: ',
     &                                             fname(:istrln(fname))
      call camxerr()
c
 7001 continue
      write(iout,'(//,A)') 'ERROR in STARTSA:'
      write(iout,'(/,1X,2A)') 'Cannot open input file: ',
     &                                             fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in STARTSA:'
      write(iout,'(/,1X,2A)') 'Cannot open output file: ',
     &                                             fname(:istrln(fname))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in STARTSA:'
      write(iout,'(/,2A)') 'Premature end-of-file reached in ',
     &                  'CAMx.in file when reading OSAT filenames.'
      write(iout,'(A)') action(:istrln(action))
      write(iout,'(2A)') 'This version of CAMx requires that you ',
     &           'provide a source area map for each nest.'
      write(iout,'(2A)') 'If you do not have a fine grid source area ',
     &                           'map, just leave the filename blank.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in STARTSA:'
      write(iout,'(/,A)') 'A source area mapping file must be ',
     &                'supplied for the coarse grid, Grid #1'
      call camxerr()
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
