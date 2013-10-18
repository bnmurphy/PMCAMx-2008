c**** READARSA.F
c
      subroutine readarsa(igrid,idate,btim,nox,noy,nspsa,dx,dy,
     &                                                        emisar)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions from each of the surface
c   emissions files and calculates the NOx and VOC levels.  It then
c   places these emissions in the appropriate place in the gridded
c   array used for tracer emissions.  The tracer to which the NOx and
c   VOC emissions is assigned depends on the source group and the
c   source region.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Argument declarations:
c        Outputs:
c          emisar R  array to store the tracer emissions
c        Inputs:
c          igrid  I  grid number for this grid
c          idate  I  date of current hour (YYJJJ)
c          btim   R  time of current hour
c          nox    I  number of X cells in grid
c          noy    I  number of Y cells in grid
c          nspsa  I  number of Y cells in grid
c          dx     R  cell width in X direction
c          dy     R  cell width in Y direction
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     06/06/96   --gwilson--    Original development
c     01/09/97   --gwilson--    Fixed a bug in fineding the source
c                               region in the fine grid
c     01/12/97   --gwilson--    Now checks for negative emissions in
c                               leftover group and exits if found
c     02/01/97   --gwilson--    Put fuzz factor of 0.1 ppb to determine
c                               if emissions are truly negative.
c     02/03/97   --gwilson--    Put code to ignore emissions in the
c                               boundary. Now includes 'bndary.com'.
c     11/06/01   --cemery--     Input dates are now Julian
c     07/19/02   --gwilson--    Seperate source area map for each grids
c     12/04/02   --gyarwood--   Improved message for negative leftover
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'grid.com'
      include 'bndary.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   igrid
      integer   idate
      real      btim
      integer   nox
      integer   noy
      integer   nspsa
      real*4    dx(noy)
      real*4    dy
      real      emisar(nox,noy,nspsa)
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   FUZZ  R  fuzz factor used to determine if emissions are truly < 0.
c
      real   FUZZ
c
      parameter( FUZZ = -0.0001 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      ndate, imap, ivoc, inox, itim
      integer      ibegcl, iendcl, i, j, k, l
      real         ttime, lefnox, lefvoc
      real         emsvoc(0:MXTEMF,MXCOLA,MXROWA)
      real         emsnox(0:MXTEMF,MXCOLA,MXROWA)
c
c-----------------------------------------------------------------------
c    Entry point:
c-----------------------------------------------------------------------
c
c   --- set the date and times ---
c
      ndate = idate
      ttime = btim/100.0
c
c   --- if up to time for releasing a new timing tracer, bump up
c       counter for number of tracer species ---
c
      if( igrid .EQ. 1 .AND. ntrtim .GT. 0 ) then
          if( MOD( INT(ttime), INT(24/ntrtim) ) .EQ. 0 ) then
              nreles = nreles + 1
              nsaspc = nsaspc + 2 * nregin
          endif
      endif
c
c   --- initialize emissions to zero ---
c
      do 71 i=1,nox
         do 81 j=1,noy
            do 91 l=1,nspsa
               emisar(i,j,l) = 0.
 91         continue
 81      continue
 71   continue
      do 12 j=1,MXROWA
         do 22 i=1,MXCOLA
            do 32 l=0,MXTEMF
              emsvoc(l,i,j) = 0.
              emsnox(l,i,j) = 0.
   32       continue
   22    continue
   12 continue
c
      call rdargrp(igrid,ndate,ttime,nox,noy,deltax(1,igrid),
     &             deltay(igrid),emsvoc,emsnox)
c
c  --- all files processed, if doing source groups calculate the 
c      left-over and put emissions into proper place in the emissions
c      array ---
c
      if( ngroup .GT. 0 ) then
         do 60 j=2,noy-1
            ibegcl = 2
            iendcl = nox-1
            if( igrid .EQ. 1 ) then
               if( ibeg(j) .EQ. -999 ) goto 60
               ibegcl = ibeg(j)
               iendcl = iend(j)
            endif
            do 70 i=ibegcl,iendcl
c
c  --- get the region for this cell from mapping array ----
c
                imap = igrmap(igrid,i,j)
                if( imap .LT. 0 .OR. imap .GT. nregin ) goto 70
c
c  --- calculate the leftover emissions to use in last source group ---
c
                lefvoc = emsvoc(0,i,j)
                lefnox = emsnox(0,i,j)
                do 80 k=1,ngroup
                   lefvoc = lefvoc - emsvoc(k,i,j) 
                   lefnox = lefnox - emsnox(k,i,j) 
c
c  --- put emissions into position in gridded tracer emissions array ---
c
                   inox = iemnox - 1 + imap + (k-1)*nregin
                   emisar(i,j,inox) = emsnox(k,i,j)
                   ivoc = iemvoc - 1 + imap + (k-1)*nregin
                   emisar(i,j,ivoc) = emsvoc(k,i,j)
   80           continue
c
c  --- put leftover emissions in last group ----
c
                if( leftovr ) then
c
c   ---- make sure leftover group is not negative ---
c
                   if( lefnox .LT. 0. ) then
                      if( lefnox .GT. FUZZ ) then
                         lefnox = 0.
                      else
                         write(iout,'(//,a)')'ERROR in READARSA:'
                         write(iout,'(/,2A,I3,A,I3,2A,I3)') 
     &                     'Negative NOx emissions calculated ',
     &                      'in leftover group in cell: (',i,',',j,') ',
     &                       'in Grid: ', igrid
                         write(iout,'(A,F20.6)') '   Value = ',lefnox
                         write(iout,'(A,I2,A,F20.6)') (' group ', k,
     &                               ' = ', emsnox(k,i,j), k=1,ngroup)
                         call camxerr()
                      endif
                   endif
                   if( lefvoc .LT. 0. ) then
                      if( lefvoc .GT. FUZZ ) then
                         lefvoc = 0.
                      else
                         write(iout,'(//,a)')'ERROR in READARSA:'
                         write(iout,'(/,2A,I3,A,I3,2A,I3)') 
     &                     'Negative VOC emissions calculated ',
     &                      'in leftover group in cell: (',i,',',j,') ',
     &                       'in Grid: ', igrid
                         write(iout,'(A,F20.6)') '    Value = ',lefvoc
                         write(iout,'(A,I2,A,F20.6)') (' group ', k,
     &                               ' = ', emsvoc(k,i,j), k=1,ngroup)
                         call camxerr()
                      endif
                   endif
                   inox = iemnox - 1 + imap + ngroup * nregin
                   emisar(i,j,inox) = lefnox
                   ivoc = iemvoc - 1 + imap + ngroup * nregin
                   emisar(i,j,ivoc) = lefvoc
                endif
c
c  --- next cell ---
c
   70       continue
   60    continue
c
c  --- only one group, just load emissions into arrays from
c      postion 0 in gridded array ----
c
      else
         do 11 j=2,noy-1
            ibegcl = 2
            iendcl = nox-1
            if( igrid .EQ. 1 ) then
               if( ibeg(j) .EQ. -999 ) goto 11
               ibegcl = ibeg(j)
               iendcl = iend(j)
            endif
            do 21 i=ibegcl,iendcl
c
c  --- get the region for this cell from mapping array ----
c
                imap = igrmap(igrid,i,j)
                if( imap .LT. 0 .OR. imap .GT. nregin ) goto 21
c
c   --- put emissions in array at correct offset ---
c
                inox = iemnox - 1 + imap
                emisar(i,j,inox) = emsnox(0,i,j)
                ivoc =  iemvoc - 1 + imap
                emisar(i,j,ivoc) = emsvoc(0,i,j)
c
c  --- next cell ---
c
   21       continue
   11    continue
      endif
c
c  --- put in the timing tracers ---
c
      do 31 k=1,nreles
         do 41 i=1,nox
            do 51 j=1,noy
c
c  --- get the region for this cell from mapping array ----
c
                imap = igrmap(igrid,i,j)
                if( imap .LT. 0 .OR. imap .GT. nregin ) goto 51
c
c   --- put emissions in array at correct offset ---
c
                itim = iemtim + (imap-1)*2 + (k-1)*2*nregin
                if( k .NE. nreles ) then
                    emisar(i,j,itim) = 0.
                    emisar(i,j,itim+1) = 0.
               else
                    emisar(i,j,itim) = TIMEMS * dx(j) * dy
                    emisar(i,j,itim+1) = TIMEMS * dx(j) * dy
               endif
   51          continue
   41      continue
   31 continue
c
      return
      end
