*** EMPREPRT
c
      subroutine empreprt(idate,btim,igrid)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads the emissions files for the RTRAC and moves the 
c   file pointer to the first hour to be used in this run.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c      01/16/02   --gwilson-- original development
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'filunit.com'
      include 'grid.com'
      include 'flags.com'
      include 'chmstry.com'
      include 'tracer.com'
      include 'rtracchm.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
      integer   igrid
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*200 fname
      character*10  namein
      character*4   ifname(10), ifnote(60), iname(10,MXSPEC)
      integer       nseg, nspect
      integer       ispc, ibgdat, iendat
      integer       idateb, idatee, nxcelt, nycelt, nzcelt, nzlowt
      integer       nzuppt, idum, izonet, iounit, isegm, npoint
      integer       iseg, i, j, ndate 
      real          tutmx, tutmy, tdxcel, tdycel, thtsur, thtlow, thtupp
      real          seg4(4), tbegti, tendti, begtim, endtim, dum
      real          txorig, tyorig, ttime
      real          txloc(MXPTSRC), tyloc(MXPTSRC), tdiam(MXPTSRC)
      logical       luse
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- set the date and time ---
c
      ndate = idate
      ttime = btim/100.0
      if( ttime .EQ. 24.0 ) then
         ndate = ndate + 1
         ttime = 0.
         if( MOD(ndate,1000) .GT. 365 ) then
            if( MOD(INT(ndate/1000),4) .EQ. 0 ) then
               if( MOD(ndate,1000) .EQ. 367 )
     &                     ndate = (INT(ndate/1000)+1)*1000 + 1
            else
               ndate = (INT(ndate/1000)+1)*1000 + 1
            endif
         endif
      endif
c
c   --- initialize the flag for PiG sources ---
c
      if( (ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG) .AND. 
     &                                               igrid .EQ. 1 ) then
         do 15 i=1,MXPTSRC
            lpigsa(i) = .FALSE.
   15    continue
      endif 
c
c   --- only read if surface emissions file is supplied ---
c
      iounit = IORTEM + igrid
      fname = temfil(igrid,1)
c
c   --- rewind the file ----
c
      rewind(iounit,ERR=7005)
c
c   --- read the header of the emissions file ---
c
      read(iounit,ERR=7000) ifname, ifnote, nseg, nspect, idateb, 
     &                                        tbegti, idatee, tendti
      read(iounit,ERR=7000) tutmx, tutmy, izonet, txorig,tyorig, 
     &                tdxcel, tdycel, nxcelt, nycelt, nzcelt, nzlowt, 
     &                                   nzuppt, thtsur, thtlow, thtupp
      read(iounit,ERR=7000) seg4
      read(iounit,ERR=7000) ((iname(i,ispc),i=1,10),ispc=1,nspect)
c
c   --- set up index for species names ---
c
      nspcem(igrid,1) = nspect
      do 20 ispc=1,nspcem(igrid,1)
         write(namein,'(10A1)') (iname(i,ispc),i=1,10)
         luse = .FALSE.
         do 30 j=1,ntotsp
           if(namein .EQ. ptname(j) ) then
               idxems(igrid,1,ispc) = j
               luse = .TRUE.
           endif
  30     continue
c
c  --- if not found in model species list, write a message ----     
c   
         if( .NOT. luse ) then
              write(idiag,'(1X,4A)') 'Species in RTRAC ',
     &              'surface emissions file: ',namein(:istrln(namein)),
     &                      ' not found in species list ... Skipping.'
         endif
  20  continue
c
c   --- read the data until the correct position in the file is reached ---
c
  111 continue
      read(iounit,ERR=7002,END=7003) ibgdat, begtim, iendat, endtim
      if( le1day ) then
         if( INT(begtim) .EQ. INT(ttime) ) goto 222
      else
         if( ibgdat .EQ. ndate .AND. 
     &                    INT(begtim) .EQ. INT(ttime) ) goto 222
      endif
      do 40 ispc=1,nspcem(igrid,1)
          read(iounit,ERR=7002,END=7003) idum, 
     &              (dum,i=1,10),((dum,i=1,ncol(igrid)),j=1,nrow(igrid))
   40 continue
      goto 111
  222 continue
      backspace(iounit)
c
c   --- point sources are only for the coarse grid -- grid #1 ---
c
      if( igrid .GT. 1 ) goto 9999
c
c   --- only read if elevated point source emissions file is supplied ---
c
      iounit = IORTPT
      fname = tptfil(1)
c
c   --- rewind the file ----
c
      rewind(iounit,ERR=7005)
c
c   --- read the header of the emissions file ---
c
      read(iounit,ERR=7000) ifname, ifnote, nseg, nspect, idateb, 
     &                                         tbegti, idatee, tendti
      read(iounit,ERR=7000) tutmx, tutmy, izonet, txorig, tyorig, 
     &               tdxcel, tdycel, nxcelt, nycelt, nzcelt, nzlowt, 
     &                                   nzuppt, thtsur, thtlow, thtupp
      read(iounit,ERR=7000) seg4
      read(iounit,ERR=7000) ((iname(i,ispc),i=1,10),ispc=1,nspect)
c
c   --- set up index for species names ---
c
      nspcpt(1) = nspect
      do 50 ispc=1,nspcpt(1)
         write(namein,'(10A1)') (iname(i,ispc),i=1,10)
         luse = .FALSE.
         do 60 j=1,ntotsp
           if(namein .EQ. ptname(j) ) then
               idxpts(1,ispc) = j
               luse = .TRUE.
           endif
  60     continue
c
c  --- if not found in model species list, write a message ----
c  
         if( .NOT. luse ) then
              write(idiag,'(1X,4A)') 'Species in RTRAC ',
     &              'point emissions file: ',namein(:istrln(namein)),
     &                      ' not found in species list ... Skipping.'
         endif
  50  continue
      do 70 iseg=1,nseg
          read(iounit,ERR=7000) isegm, npoint
          if( npoint .LE. 0 ) goto 70
          read(iounit,ERR=7002) (txloc(i),tyloc(i),idum,tdiam(i),
     &                                             dum,dum,i=1,npoint)
   70 continue
      do 80 i=1,npoint
         xlocpt(i) = txloc(i)
         ylocpt(i) = tyloc(i)
         if( (ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG) 
     &                      .AND. tdiam(i) .LT. 0. ) lpigsa(i) = .TRUE.
   80 continue
c
c   --- read the data until the correct position in the file is reached ---
c
  333 continue
      read(iounit,ERR=7002,END=7003) ibgdat, begtim, iendat, endtim
      if( le1day ) then
         if( begtim .EQ. ttime ) goto 444
      else
         if( ibgdat .EQ. ndate .AND. begtim .EQ. ttime ) goto 444
      endif
      read(iounit,ERR=7002) isegm, npoint
      read(iounit,ERR=7002) (idum,idum,idum,dum,dum,i=1,npoint)
      do 90 ispc=1,nspcpt(1)
          read(iounit,ERR=7002,END=7003) idum,(dum,i=1,10),
     &                                                (dum,i=1,npoint)
   90 continue
      goto 333
  444 continue
      backspace(iounit)
c
c  ---- return to calling routine ---
c
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,3A)') 'Reading header of ',
     &                     'emissions file: ',fname(:istrln(fname))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,2A)') 'Reading emissions grid in file: ',
     &                                        fname(:istrln(fname))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,3A)') 'Premature end-of-file reading ',
     &                     'emissions file: ',fname(:istrln(fname))
      call camxerr()
c
 7005 continue
      write(iout,'(//,A)') 'ERROR in EMPREPRT:'
      write(iout,'(/,1X,3A)') 'Cannot access emissions file',
     &     ' provided for OSAT processing: ',fname(:istrln(fname))
      call camxerr()
c
c-----------------------------------------------------------------------
c    Format statements:
c-----------------------------------------------------------------------
c
c-----------------------------------------------------------------------
c    Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
