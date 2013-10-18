c**** READPTSA
c
      subroutine readptsa(idate,btim)
c
c-----CAMx v4.02 030709
c
c-----------------------------------------------------------------------
c    Description:
c-----------------------------------------------------------------------
c
c   This routine reads one hour of emissions from each of the elevated
c   emissions files and calculates the NOx and VOC levels.  It then
c   places these emissions in the appropriate place in the point source
c   array used for tracer emissions.  The tracer to which the NOx and
c   VOC emissions is assigned depends on the source group and the
c   source region.
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Inputs:
c          idate  I  date of current hour (YYJJJ)
c          btim   R  time of current hour
c
c-----------------------------------------------------------------------
c    LOG:
c-----------------------------------------------------------------------
c
c     06/06/96   --gwilson--    Original development
c     12/10/96   --gwilson--    Added switch to skip emissions if PiG
c     02/01/97   --gwilson--    Put fuzz factor of 0.1 ppb to catch if
c                               emissions of leftover group is negative.
c     02/03/97   --gwilson--    Put code to ignore emissions on boundary.
c     01/20/99   --cemery---    Grid cell size from file should be meters
c                               for all cartesian projections (UTM, LCP,
c                               PSP)
c     11/06/01   --cemery--     Input dates are now Julian
c     07/19/02   --gwilson--    Added seperate source area map for each grids.
c
c-----------------------------------------------------------------------
c    Include files:
c-----------------------------------------------------------------------
c
      include 'camx.prm'
      include 'camx.com'
      include 'flags.com'
      include 'grid.com'
      include 'ptemiss.com'
      include 'tracer.com'
      include 'filunit.com'
c
c-----------------------------------------------------------------------
c    Arguement declarations:
c-----------------------------------------------------------------------
c
      integer   idate
      real      btim
c
c-----------------------------------------------------------------------
c    Local parameters:
c-----------------------------------------------------------------------
c
c   FUZZ   R    fuzz factor used to determine if emissions are truly < 0
c
      real   FUZZ
c
      parameter( FUZZ = -0.0001 )
c
c-----------------------------------------------------------------------
c    Local variables:
c-----------------------------------------------------------------------
c
      integer      imap, ivoc, inox, i, j, k, ii, jj
      integer      ndate, icel, jcel, izcel(MXPTSRC), xloccrs, yloccrs
      real         lefnox, lefvoc, xloctmp, yloctmp
      real         ttime
      real         emsvoc(0:MXTEMF,MXPTSRC), emsnox(0:MXTEMF,MXPTSRC)
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
c   --- initialize emissions to zero ---
c
      do 11 i=1,MXPTSRC
         do 21 j=1,nsaspc
            sapnts(i,j) = 0.
 21      continue
         do 31 j=0,MXTEMF
            emsvoc(j,i) = 0. 
            emsnox(j,i) = 0. 
 31      continue
 11   continue
c
      call rdptgrp(ndate,ttime,emsvoc,emsnox,izcel)
c
c  --- all files processed, if doing source groups calculate the 
c      left-over and put emissions into proper place in the emissions
c      array ---
c
      if( ngroup .GT. 0 ) then
         do 60 i=1,nptsrc
c
c  --- get the region for this cell from mapping array ----
c
            if( .NOT.llatlon ) then
               icel = INT( (xlocpt(i)/1000. - xorg) / delx ) + 1
               jcel = INT( (ylocpt(i)/1000. - yorg) / dely ) + 1
               xloccrs = xlocpt(i)/1000. - xorg
               yloccrs = ylocpt(i)/1000. - yorg
            else
               icel = INT( (xlocpt(i) - xorg) / delx ) + 1
               jcel = INT( (ylocpt(i) - yorg) / dely ) + 1
               xloccrs = xlocpt(i) - xorg
               yloccrs = ylocpt(i) - yorg
            endif
            if(icel .LE. 0 .OR. icel .GT. ncol(1)) goto 60
            if(jcel .LE. 0 .OR. jcel .GT. nrow(1)) goto 60
c  
c-----The following is specific to BAAQMD application  
c  
c     elseif (lpolar) then  
c       do i = 1,nptsrc  
c         call pspgeo(0,polelon,polelat,xstk,ystk,xlocpt(i),ylocpt(i))
c         icel = INT((xstk/1000. - xorg)/delx) + 1
c         jcel = INT((ystk/1000. - yorg)/dely) + 1
c       enddo 
c
c
c   --- find out if a nest contains this source  ---
c
            igrd = 1
            do ig = 2,ngrid
              xloctmp = xloccrs - (inst1(ig)-1)*delx
              yloctmp = yloccrs - (jnst1(ig)-1)*dely
              if( xloctmp .GT. 0. .AND. yloctmp .GT. 0. ) then
                  ii = 2 + INT( xloctmp/delx * FLOAT( meshold(ig) ) )
                  jj = 2 + INT( yloctmp/dely * FLOAT( meshold(ig) ) )
                  if( ii .GT. 1 .AND. jj .GT. 1 .AND. ii .LT. ncol(ig) 
     &                                    .AND. jj .LT. nrow(ig) ) then
                     igrd = ig
                     icel = ii
                     jcel = jj
                   endif
              endif
            enddo
            imap = igrmap(igrd,icel,jcel)
            if( izcel(i) .LT. 0 ) imap = ABS( izcel(i) )
            if( imap .LT. 0 .OR. imap .GT. nregin ) goto 60
c
c  --- calculate the leftover emissions to use in last source group ---
c
            lefvoc = emsvoc(0,i)
            lefnox = emsnox(0,i)
            do 70 k=1,ngroup
               lefvoc = lefvoc - emsvoc(k,i) 
               lefnox = lefnox - emsnox(k,i) 
c
c  --- put emissions into position in gridded tracer emissions array ---
c
               inox = iemnox - 1 + imap + (k-1)*nregin
               sapnts(i,inox) = emsnox(k,i)
               ivoc = iemvoc - 1 + imap + (k-1)*nregin
               sapnts(i,ivoc) = emsvoc(k,i)
   70       continue
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
                    write(iout,'(//,a)') 'ERROR in READPTSA:'
                    write(iout,'(/,2A,I5,A)') 
     &                     'Negative NOx emissions calculated ',
     &                     'in leftover group in point source: ',i,' :'
                    write(iout,'(A,F20.6)') '   Value = ',lefnox
                    call camxerr()
                  endif
               endif
               if( lefvoc .LT. 0. ) then
                  if( lefvoc .GT. FUZZ ) then
                    lefvoc = 0.
                  else
                    write(iout,'(//,a)') 'ERROR in READPTSA:'
                    write(iout,'(/,2A,I5,A)') 
     &                   'Negative VOC emissions calculated ',
     &                   'in leftover group in point source: ',i,' :'
                    write(iout,'(A,F20.6)') '   Value = ',lefvoc
                    call camxerr()
                  endif
               endif
               inox = iemnox - 1 + imap + ngroup * nregin
               sapnts(i,inox) = lefnox
               ivoc = iemvoc - 1 + imap + ngroup * nregin
               sapnts(i,ivoc) = lefvoc
             endif
c
c  --- next cell ---
c
   60    continue
c
c  --- only one group, just load emissions into arrays from
c      position 0 in gridded array ----
c
      else
         do 90 i=1,nptsrc
c
c  --- get the region for this cell from mapping array ----
c
            if( .NOT.llatlon ) then
               icel = INT( (xlocpt(i)/1000. - xorg) / delx ) + 1
               jcel = INT( (ylocpt(i)/1000. - yorg) / dely ) + 1
               xloctmp = xlocpt(i)/1000. - xorg
               yloctmp = ylocpt(i)/1000. - yorg
            else
               icel = INT( (xlocpt(i) - xorg) / delx ) + 1
               jcel = INT( (ylocpt(i) - yorg) / dely ) + 1
               xloctmp = xlocpt(i) - xorg
               yloctmp = ylocpt(i) - yorg
            endif
            if(icel .LE. 0 .OR. icel .GT. ncol(1)) goto 90
            if(jcel .LE. 0 .OR. jcel .GT. nrow(1)) goto 90
c
c   --- find out if a nest contains this source  ---
c
            igrd = 1
            do ig = 2,ngrid
              xloctmp = xloctmp - (inst1(ig)-1)*delx
              yloctmp = yloctmp - (jnst1(ig)-1)*dely
              ii = 2 + INT( xloctmp/delx * FLOAT( meshold(ig) ) )
              jj = 2 + INT( yloctmp/dely * FLOAT( meshold(ig) ) )
              if( ii .GT. 1 .AND. jj .GT. 1 .AND. ii .LT. ncol(ig) .AND.
     &                                           jj .LT. nrow(ig) ) then
                 igrd = ig
                 icel = ii
                 jcel = jj
               endif
            enddo
            imap = igrmap(igrd,icel,jcel)
            if( izcel(i) .LT. 0 ) imap = ABS( izcel(i) )
            if( imap .LT. 0 .OR. imap .GT. nregin ) goto 90
c
c   --- put emissions in array at correct offset ---
c
            inox = iemnox - 1 + imap
            sapnts(i,inox) = emsnox(0,i)
            ivoc = iemvoc - 1 + imap
            sapnts(i,ivoc) = emsvoc(0,i)
c
c  --- next cell ---
c
   90   continue
      endif
c
      return
      end
