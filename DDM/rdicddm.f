c**** RDICDDM.F
c
      subroutine rdicddm(nox,noy,noz,nspsa,saconc)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine fils one hour of initial conditions for the DDM
c   species. It then places these concentrations in the  appropriate 
c   place in the gridded array used for tracer concentrations.  
c   The values are left in PPM so that they can be interpolated to
c   the nests in a later step.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c      Argument description:
c       Outputs:
c           saconc   R  tracer concentrations
c       Inputs:
c           nox      I  number of X cells in the grid
c           noy      I  number of Y cells in the grid
c           noz      I  number of layers in the grid
c           nspsa    I  number of tracer species
c       
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     03/31/99   --gwilson--    Original development
c     10/24/01  Removed BSWAP and converted integer strings to character*4
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'filunit.com'
      include 'chmstry.com'
      include 'grid.com'
      include 'flags.com'
      include 'bndary.com'
      include 'tracer.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   nox
      integer   noy
      integer   noz
      integer   nspsa
      real      saconc(nox,noy,noz,nspsa)
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      character*10 aqfil, cfile, cspec
      character*4  ifile(10), note(60), ispec(10)
      integer      ispc, ifirst, iddm, i, j, n
      integer      iseg, nspc, idat1, idat2
      integer      izone, nx, ny, nz, ilay, idum, imod
      real         tim1, tim2, orgx, orgy, utmx, utmy, dx, dy
      real         concic(MXCOLA,MXROWA)
      logical      lexist, luse
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data aqfil /'AIRQUALITY'/

c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- open the IC ---
c
      inquire(file=icfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      open(file=icfil,unit=IORIC,status='UNKNOWN',
     &                                   form='UNFORMATTED',ERR=7001)
c
c  --- call routine to zero out the array ---
c
      call zeros(saconc,nox*noy*noz*nspsa)
c
c  --- read 1st IC header record and check inputs ---
c
      read(IORIC,ERR=7002) ifile,note,iseg,nspc,idat1,tim1,idat2,tim2
      if( INT(tim2) .EQ. 24 ) then
        idat2 = idat2 + 1
        tim2 = 0.
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      write(cfile,'(10A1)') (ifile(i),i=1,10)
      if( cfile .NE. aqfil ) goto 7003
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if( idat1 .GT. begdate .AND. .NOT. le1day ) goto 7004
      if( idat1 .EQ. begdate .AND. tim1 .GT. begtim
     &                              .AND. .NOT. le1day ) goto 7004
c
c  --- read 2nd IC header record and check inputs ---
c
      read(IORIC,ERR=7002) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz
      if( .NOT. llatlon ) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if( nx .NE. ncol(1) .OR. ny .NE. nrow(1) 
     &                                .OR. nz .NE. nlay(1) ) goto 7005
c
c  --- skip the next 2 records ---
c
      read(IORIC,ERR=7002)
      read(IORIC,ERR=7002) 
c
c  --- read the date and time and make sure it is the correct hour ---
c
  111 continue
      read(IORIC,ERR=7007,END=7006) idat1, tim1, idat2, tim2
      tim1 = 100*tim1
      tim2 = 100*tim2
      if( INT(tim2) .EQ. 24 ) then
        idat2 = idat2 + 1
        tim2 = 0.
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                     idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      if( (idat1 .LT. begdate .OR. 
     &       (idat1 .EQ. begdate .AND. tim1 .LE. begtim)) .AND.
     &            (idat2 .GT. begdate .OR. 
     &                (idat2 .EQ. begdate .AND. tim2 .GT. begtim))) then
c
c  --- read the concentrations for each species ---
c
         do ispc=1,nspc
c
c   --- loop over layers ---
c
            do ilay=1,nz
               read(IORIC) idum,(ispec(n),n=1,10),
     &                                 ((concic(i,j),i=1,nx),j=1,ny)
               write(cspec,'(10A1)') ispec
               if(cspec .EQ. 'HNO2      ') cspec = 'HONO      '
               if(cspec .EQ. 'HCHO      ' .and. kHCHO.eq.nspec+1) 
     &                                            cspec = 'FORM      '
c
c  --- find this species in the modeled species list ---
c
               do imod=1,nspec
                  if( cspec .EQ. spname(imod) ) then
                      ifirst = iptddm(imod) 
                      do iddm=1,nicddm
                         luse = .FALSE.
                         if( icddmsp(iddm) .EQ. cspec ) luse = .TRUE. 
                         if( icddmsp(iddm) .EQ. NAMALL  ) luse = .TRUE. 
                         if( icddmsp(iddm) .EQ. NAMVOC 
     &                               .AND. lvocsp(imod) ) luse = .TRUE. 
                         if( icddmsp(iddm) .EQ. NAMNOX 
     &                               .AND. lnoxsp(imod) ) luse = .TRUE. 
c
c  --- if this DDM species matches this modeled species, load it ---
c
                         if( luse ) then
                            do 10 j=2,ny-1
                               icel1 = 2
                               icel2 = nx - 1
                               if( ibeg(j) .EQ. -999 ) goto 10
                               icel1 = ibeg(j)
                               icel2 = iend(j)
                               do i=icel1,icel2
                                  if( concic(i,j) .GE. bdnl(imod) )
     &                               saconc(i,j,ilay,ifirst+iddm-1) = 
     &                                                       concic(i,j)
                               enddo
   10                       continue
                         endif
c
c  --- check the next IC DDM species ---
c
                      enddo
c
c  --- check the next modeled species ---
c
                  endif
               enddo
c
c  --- next layer ---
c
            enddo
c
c  --- read nect species worth of data ---
c
         enddo
c
c  --- if not the right hour, read through it ---
c
      else
         do ispc=1,nspc
            do ilay=1,nz
               read(IORIC,ERR=7007) 
            enddo
         enddo
         goto 111
      endif
c
c  --- close the file and return to the calling routine ---
c
      close(IORIC)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue 
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  IC file for DDM does not exist: ',
     &                                             icfil(:istrln(icfil))
      call camxerr()
c
 7001 continue 
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  Opening IC file for DDM: ',
     &                                             icfil(:istrln(icfil))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  Reading header of IC file for DDM: ',
     &                                             icfil(:istrln(icfil))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  IC input file for DDM is not labelled ',
     &                                                 'AIRQUALITY.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  IC input file for DDM is not for correct ',
     &                                                   'time period.'
      write(iout,*) '   ***  Episode ***'
      write(iout,'(a,i10.5)') 'Date   : ',begdate
      write(iout,'(a,f10.1)') 'Time   : ',begtim
      write(iout,*) '   ***  IC File ***'
      write(iout,'(a,i10.5)') 'Date   : ',idat1
      write(iout,'(a,f10.1)') 'Time   : ',tim1
      call camxerr()
c
 7005 continue
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  IC input file does not appear to be for ',
     &                ' the correct grid.' 
      write(iout,*) '   ***  Coarse Grid ***'
      write(iout,*) 'No. of Cells : (',ncol(1),',',nrow(1),')'
      write(iout,*) 'No. of layers:  ',nlay(1)
      write(iout,*) '   ***  IC File ***'
      write(iout,*) 'No. of Cells : (',nx,',',ny,')'
      write(iout,*) 'No. of layers:  ',nz
      call camxerr()
c
 7006 continue
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  Premature end-of-file reached in IC file ',
     &                                                      'for DDM.'
      write(iout,*) 'Make sure the file contains the correct date/time.'
      call camxerr()
c
 7007 continue
      write(iout,'(//,A)') 'ERROR in RDICDDM:'
      write(iout,*) 'ERROR:  Reading IC file for DDM: ',
     &                                             icfil(:istrln(icfil))
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
