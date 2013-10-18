c**** RDBCRT.F
c
      subroutine rdbcrt(nsteps,nox,noy,noz,nspas,saconc,tpgrd,prgrd)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine fills one hour of boundary conditions for the RTRAC
c   species. It then places these concentrations in the  appropriate 
c   place in the gridded array used for tracer concentrations.  
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c      Argument description:
c       Outputs:
c           saconc   R  tracer concentrations
c       Inputs:
c           nsteps   I  number of steps already performed
c           nox      I  number of X cells in the grid
c           noy      I  number of Y cells in the grid
c           noz      I  number of layers in the grid
c           nspas    I  number of tracer species
c           tpgrd    I  3-D temperature field
c           prgrd    I  3-D pressure field
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     01/21/02   --gwilson--    Original development
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
      include 'rtracchm.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   nsteps
      integer   nox
      integer   noy
      integer   noz
      integer   nspas
      real      saconc(nox,noy,noz,nspas)
      real      tpgrd(nox,noy,noz)
      real      prgrd(nox,noy,noz)
c
c-----------------------------------------------------------------------
c    Local variables:
c----------------------------------------------------------------------
c
      character*10 bnfil, cfile, cspec
      character*4  ifile(10), note(60), ispec(10)
      integer      icl, jcl, ispc, iedge
      integer      iseg, nspc, idat1, idat2
      integer      izone, nx, ny, nz, imod
      integer      i, j, k
      real         tim1, tim2, orgx, orgy, utmx, utmy, dx, dy
      real         conwst(MXLAYA,MXROWA), conest(MXLAYA,MXROWA)
      real         consth(MXLAYA,MXCOLA), connth(MXLAYA,MXCOLA)
      logical      lexist, luse
c
c-----------------------------------------------------------------------
c    Data statements:
c-----------------------------------------------------------------------
c
      data bnfil /'BOUNDARY'/
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c  --- initialize the boundary conditions to lower bound ---
c
      do ispc=1,ntotsp
        do k=1,noz
          do j=1,noy
            icl = ibeg(j) - 1
            if( icl .GT. 0 .AND. icl .LE. nox ) saconc(icl,j,k,ispc) = 
     &                                                      rtlbnd(ispc)
            icl = iend(j) + 1
            if( icl .GT. 0 .AND. icl .LE. nox ) saconc(icl,j,k,ispc) = 
     &                                                      rtlbnd(ispc)
          enddo
          do i=1,nox
            jcl = jbeg(i) - 1
            if( jcl .GT. 0 .AND. jcl .LE. noy ) saconc(i,jcl,k,ispc) =
     &                                                      rtlbnd(ispc)
            jcl = jend(i) + 1
            if( jcl .GT. 0 .AND. jcl .LE. noy ) saconc(i,jcl,k,ispc) =
     &                                                      rtlbnd(ispc)
          enddo
        enddo
      enddo
c
c  --- open the BC ---
c
      inquire(file=bcfil,exist=lexist)
      if( .NOT. lexist ) goto 7000
      open(file=bcfil,unit=IORBC,status='UNKNOWN',
     &                                   form='UNFORMATTED',ERR=7001)
c
c  --- read 1st BC header record and check inputs ---
c
      read(IORBC,ERR=7002) ifile,note,iseg,nspc,idat1,tim1,idat2,tim2
      if( INT(tim2) .EQ. 24 ) then
        idat2 = idat2 + 1
        tim2 = 0.
        if( MOD(idat2,1000) .GT. 365 ) then
           if( MOD(INT(idat2/1000),4) .EQ. 0 ) then
              if( MOD(idat2,1000) .EQ. 367 )
     &                       idat2 = (INT(idat2/1000)+1)*1000 + 1
           else
              idat2 = (INT(idat2/1000)+1)*1000 + 1
           endif
        endif
      endif
      write(cfile,'(10A1)') (ifile(i),i=1,10)
      if( cfile .NE. bnfil ) goto 7003
      tim1 = 100.*tim1
      tim2 = 100.*tim2
      if( idat1 .GT. date .AND. .NOT. le1day ) goto 7004
      if( idat1 .EQ. date .AND. tim1 .GT. time
     &                              .AND. .NOT. le1day ) goto 7004
c
c  --- read 2nd BC header record and check inputs ---
c
      read(IORBC,ERR=7002) orgx,orgy,izone,utmx,utmy,dx,dy,nx,ny,nz
      if( .NOT. llatlon ) then
        dx = dx/1000.
        dy = dy/1000.
      endif
      if( nx .NE. ncol(1) .OR. ny .NE. nrow(1) 
     &                                .OR. nz .NE. nlay(1) ) goto 7005
c
c  --- skip the next 6 records ---
c
      read(IORBC,ERR=7002)
      read(IORBC,ERR=7002) 
      do i=1,4
        read(IORBC,ERR=7002) 
      enddo
c
c  --- read the date and time and make sure it is the correct hour ---
c
  111 continue
      read(IORBC,ERR=7007,END=7006) idat1, tim1, idat2, tim2
      tim1 = tim1*100
      tim2 = tim2*100
      if( INT(tim2) .EQ. 24 ) then
        tim2 = 0.
        idat2 = idat2 + 1
      endif
      if( (idat1 .LT. date .OR. 
     &       (idat1 .EQ. date .AND. tim1 .LE. time)) .AND.
     &            (idat2 .GT. date .OR. 
     &                (idat2 .EQ. date .AND. tim2 .GT. time))) then
c
c  --- read the concentrations for each species ---
c
         do ispc=1,nspc
            read(IORBC,ERR=7007,END=7006) iseg, (ispec(i),i=1,10),iedge,
     &                                     ((conwst(k,j),k=1,nz),j=1,ny)
            read(IORBC,ERR=7007,END=7006) iseg, (ispec(i),i=1,10),iedge,
     &                                     ((conest(k,j),k=1,nz),j=1,ny)
            read(IORBC,ERR=7007,END=7006) iseg, (ispec(i),i=1,10),iedge,
     &                                     ((consth(k,i),k=1,nz),i=1,nx)
            read(IORBC,ERR=7007,END=7006) iseg, (ispec(i),i=1,10),iedge,
     &                                     ((connth(k,i),k=1,nz),i=1,nx)
            write(cspec,'(10A1)') ispec
c
c  --- find this species in the modeled species list ---
c
            luse = .FALSE.
            do imod=1,ntotsp
               if( cspec .EQ. ptname(imod) ) then
                  luse = .TRUE.
                  do k=1,nz
                     do j=1,ny
                        icl = ibeg(j) - 1 
                        if( icl .GT. 0 .AND. icl .LE. nx ) then
c
c  --- get the conversion factor to umol/m3 ---
c
                           if( imod .LE. nrtgas ) then
                              convfac = densfac*273./
     &                               tpgrd(icl,j,k)*prgrd(icl,j,k)/1013.
                           else
                              convfac = 1.0
                           endif
c
c  --- load the west boundary ---
c
                           saconc(icl,j,k,imod) = conwst(k,j) * convfac
                           if( saconc(icl,j,k,imod) .LT. rtlbnd(imod) )
     &                              saconc(icl,j,k,imod) = rtlbnd(imod)
                        endif
c
c  --- load the east boundary ---
c
                        icl = iend(j) + 1
                        if( icl .GT. 0 .AND. icl .LE. nx ) then
c
c  --- get the conversion factor to umol/m3 ---
c
                           if( imod .LE. nrtgas ) then
                                convfac = densfac*273./
     &                               tpgrd(icl,j,k)*prgrd(icl,j,k)/1013.
                           else
                                convfac = 1.0
                           endif
                           saconc(icl,j,k,imod) = conest(k,j) * convfac
                           if( saconc(icl,j,k,imod) .LT. rtlbnd(imod) )
     &                               saconc(icl,j,k,imod) = rtlbnd(imod)
                        endif
                     enddo
c
c  --- load the south boundary ---
c
                     do i=1,nx
                       jcl = jbeg(i) - 1 
                       if( jcl .GT. 0 .AND. jcl .LE. ny ) then
c
c  --- get the conversion factor to umol/m3 ---
c
                           if( imod .LE. nrtgas ) then
                              convfac = densfac*273./
     &                               tpgrd(i,jcl,k)*prgrd(i,jcl,k)/1013.
                           else
                              convfac = 1.0
                           endif
                           saconc(i,jcl,k,imod) = consth(k,j) * convfac
                           if( saconc(i,jcl,k,imod) .LT. rtlbnd(imod) )
     &                               saconc(i,jcl,k,imod) = rtlbnd(imod)
                       endif
c
c  --- load the east boundary ---
c
                       jcl = jend(i) + 1
                       if( jcl .GT. 0 .AND. jcl .LE. ny ) then
c
c  --- get the conversion factor to umol/m3 ---
c
                          if( imod .LE. nrtgas ) then
                              convfac = densfac*273./
     &                              tpgrd(i,jcl,k)*prgrd(i,jcl,k)/1013.
                          else
                              convfac = 1.0
                          endif
                          saconc(i,jcl,k,imod) = connth(k,j) * convfac
                          if( saconc(i,jcl,k,imod) .LT. rtlbnd(imod) )
     &                               saconc(i,jcl,k,imod) = rtlbnd(imod)
                       endif
                     enddo
c
c  --- get the next layer ---
c
                  enddo
c
c  --- check the next modeled species ---
c
               endif
            enddo
c
c  --- if not found in model species list, write a message ----   
c   
            if( .NOT. luse .AND. nsteps .EQ. 1 ) then
                write(idiag,'(1X,4A)') 'Species in RTRAC ',     
     &              'boundary conditions file: ',cspec(:istrln(cspec)),
     &                      ' not found in species list ... Skipping.'
            endif
c
c  --- read next species worth of data ---
c
         enddo
c
c  --- if not the right hour, read through it ---
c
      else
         do ispc=1,nspc
            do i=1,4
               read(IORBC,ERR=7007) 
            enddo
         enddo
         goto 111
      endif
c
c  --- finally got the right hour, close the file ---
c
      close(IORBC)
c
c  --- open the TOPCONC file ---
c
      inquire(file=tcfil,exist=lexist)
      if( .NOT. lexist ) goto 7008
      open(file=tcfil,unit=IORTC,status='UNKNOWN',ERR=7009)
c
c  --- load the TOP concentrations into the appropriate place ---
c
      do imod=1,ntotsp
         ptloft(imod) = rtlbnd(imod)
      enddo
 222  continue
      read(IORTC,'(A10,F10.0)',END=333) cspec, ctin
      imod = 0
      do i=1,nspec
         if( cspec .EQ. ptname(i) ) imod = i
      enddo
      if( imod .LE. 0 ) then 
         if( nsteps .EQ. 1 ) then
            write(idiag,'(1X,4A)') 'Species in RTRAC ',
     &              'topconc file: ',cspec(:istrln(cspec)),
     &                      ' not found in species list ... Skipping.'
         endif
         goto 222
      endif
      if( ctin .GT. rtlbnd(imod) ) then
          ptloft(imod) = ctin
      endif
      goto 222
c
c  --- close the file and return to the calling routine ---
c
  333 continue
      close(IORTC)
      goto 9999
c
c-----------------------------------------------------------------------
c    Error messages:
c-----------------------------------------------------------------------
c
 7000 continue 
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  BC file for RTRAC does not exist: ',
     &                                             bcfil(:istrln(bcfil))
      call camxerr()
c
 7001 continue 
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  Opening BC file for RTRAC: ',
     &                                            bcfil(:istrln(bcfil))
      call camxerr()
c
 7002 continue
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  Reading header of BC file for RTRAC: ',
     &                                            bcfil(:istrln(bcfil))
      call camxerr()
c
 7003 continue
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  BC input file for RTRAC is not labelled ',
     &                                                 'BOUNDARY.'
      call camxerr()
c
 7004 continue
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  BC input file for RTRAC is not for ',
     &                                            'corect time period.'
      write(iout,*) '   ***  Episode ***'
      write(iout,'(a,i10.5)') 'Date   : ',date
      write(iout,'(a,f10.1)') 'Time   : ',time
      write(iout,*) '   ***  BC File ***'
      write(iout,'(a,i10.5)') 'Date   : ',idat1
      write(iout,'(a,f10.1)') 'Time   : ',tim1
      call camxerr()
c
 7005 continue
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  BC input file does not appear to be for ',
     &                ' the correct grid.' 
      write(iout,*) '   ***  Coarse Grid ***'
      write(iout,*) 'No. of Cells : (',ncol(1),',',nrow(1),')'
      write(iout,*) 'No. of layers:  ',nlay(1)
      write(iout,*) '   ***  BC File ***'
      write(iout,*) 'No. of Cells : (',nx,',',ny,')'
      write(iout,*) 'No. of layers:  ',nz
      call camxerr()
c
 7006 continue
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  Premature end-of-file reached in BC file ',
     &                                                      'for RTRAC.'
      write(iout,*) 'Make sure the file contains the correct date/time.'
      call camxerr()
c
 7007 continue
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  Reading BC file for RTRAC: ',
     &                                          bcfil(:istrln(bcfil))
      call camxerr()
c
 7008 continue 
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  TC file for RTRAC does not exist: ',
     &                                          tcfil(:istrln(tcfil))
      call camxerr()
c
 7009 continue 
      write(iout,'(//,A)') 'ERROR in RDBCRT:'
      write(iout,*) 'ERROR:  Opening TopConc file for RTRAC: ',
     &                                          tcfil(:istrln(tcfil))
      call camxerr()
c
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
