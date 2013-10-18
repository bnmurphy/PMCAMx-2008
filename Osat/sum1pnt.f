      subroutine sum1pnt(igroup,idx,nptsrc,emsbas,emsoth,
     &                   emslft,emstot,emspnt,emssum)
c
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c
c-----CAMx v4.02 030709
c
c     SUM1PNT sums up the point emission of one species for a given group
c
c       07/19/02  --gwilso-- Added seperate source area map for each grids.
c
c
c     Input argument:
c        igroup            group ID
c        idx               specie ID
c        emspnt            the species emission for the group
c
c     Output arguments:
c        emssum            emission summed over grid
c        emsbas            base emission
c        emsoth            "otherwise" emission
c        emslft            leftover emission
c        emstot            total emission
c
      include "camx.prm"
      include "grid.com"
      include "flags.com"
      include "tracer.com"
c
c
      real emspnt(nptsrc),xloctmp(MXPTSRC),yloctmp(MXPTSRC)
      real emsbas(MXSPEC,MXTRSP), emsoth(MXSPEC,MXTRSP)
      real emssum(MXSPEC,MXTRSP)
      real emslft(MXCOLA,MXROWA),emstot(MXCOLA,MXROWA)
      real xlocnst, ylocnst
c
c   --- sum up the emissions for each point ---
c
c     if (igroup.eq.0) then
        if (llatlon) then
          do n = 1,nptsrc
            xloctmp(n) = xlocpt(n) - xorg
            yloctmp(n) = ylocpt(n) - yorg
          enddo
        else
          do n = 1,nptsrc
            xloctmp(n) = xlocpt(n)/1000. - xorg
            yloctmp(n) = ylocpt(n)/1000. - yorg
          enddo
c 
c-----The following is specific to BAAQMD application 
c 
c     elseif (lpolar) then 
c       do n = 1,nptsrc 
c         call pspgeo(0,polelon,polelat,xstk,ystk,xlocpt(n),ylocpt(n)) 
c         xloctmp(n) = xstk/1000. - xorg 
c         yloctmp(n) = ystk/1000. - yorg 
c       enddo
        endif
c     endif
c
      do 70 i=1,nptsrc
         icel = 1 + INT( xloctmp(i)/delx )
         jcel = 1 + INT( yloctmp(i)/dely )
         if(icel .LE. 0 .OR. icel .GT. ncol(1)) goto 70
         if(jcel .LE. 0 .OR. jcel .GT. nrow(1)) goto 70
c
c   --- find out if a nest contains this source  ---
c
         igrd = 1
         do ig = 2,ngrid
           xlocnst = xloctmp(i) - (inst1(ig)-1)*delx
           ylocnst = yloctmp(i) - (jnst1(ig)-1)*dely
           ii = 2 + INT( xlocnst/delx * FLOAT( meshold(ig) ) )
           jj = 2 + INT( ylocnst/dely * FLOAT( meshold(ig) ) )
           if( ii .GT. 1 .AND. jj .GT. 1 .AND. ii .LT. ncol(ig) .AND.
     &                                           jj .LT. nrow(ig) ) then
              igrd = ig
              icel = ii
              jcel = jj
            endif
         enddo
c
c  --- get the region for this cell from mapping array ----
c
         imap = igrmap(igrd,icel,jcel)
         if( imap .LT. 0 .OR. imap .GT. nregin ) goto 70
c
c  --- calculate the index into the tracer species for this gruoup/region ---
c
         if( ngroup .GT. 0 ) then
c
c   --- if group is base emissions, add to "leftover" group ----
c
            if( igroup .EQ. 0 ) then
               if( leftovr ) then
                  inox = iemnox - 1 + imap + ngroup*nregin
                  emsbas(idx,inox) = emsbas(idx,inox) + emspnt(i)
c
c  --- if doing PiG and source is a PIG source, set the PiG species index ---
c
                  if( (ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG) 
     &                          .AND. lpigsa(i) .AND. lnoxsp(idx) .AND. 
     &                                              emspnt(i) .GT. 0. )
     &                        ipigsp(i) = imap + ngroup*nregin
c
c  --- add emissions to totals array ----
c
                  ivoc = iemvoc - 1 + imap + ngroup*nregin
                  emsbas(idx,ivoc) = emsbas(idx,ivoc) + emspnt(i)
               endif
               emstot(icel,jcel) = emstot(icel,jcel) +  emspnt(i)
c
c   --- otherwise, add to this group/region and subtract from "leftover" ---
c
            else
               inox = iemnox - 1 + imap + (igroup-1)*nregin
               emssum(idx,inox) = emssum(idx,inox) + emspnt(i)
c
c  --- if doing PiG and source is a PIG source, set the PiG species index ---
c
               if( (ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG) 
     &                         .AND. lpigsa(i) .AND. lnoxsp(idx) .AND. 
     &                                              emspnt(i) .GT. 0. )
     &                     ipigsp(i) = imap + (igroup-1)*nregin
c
c  --- add emissions to totals array ----
c
               ivoc = iemvoc - 1 + imap + (igroup-1)*nregin
               emssum(idx,ivoc) = emssum(idx,ivoc) + emspnt(i)
               if( leftovr ) then
                  inox = iemnox - 1 + imap + ngroup*nregin
                  emsoth(idx,inox) = emsoth(idx,inox) + emspnt(i)
                  ivoc = iemvoc - 1 + imap + ngroup*nregin
                  emsoth(idx,ivoc) = emsoth(idx,ivoc) + emspnt(i)
               endif
               emslft(icel,jcel)= emslft(icel,jcel) + emspnt(i)
            endif
c
c   --- only using regular model emissions ---
c
         else
            inox = iemnox - 1 + imap
            emssum(idx,inox) = emssum(idx,inox) + emspnt(i)
c
c  --- if doing PiG and source is a PIG source, set the PiG species index ---
c
            if( (ipigflg .EQ. GRESPIG .OR. ipigflg .EQ. IRONPIG) 
     &       .AND. lpigsa(i) .AND. lnoxsp(idx) .AND. emspnt(i) .GT. 0. )
     &                                                  ipigsp(i) = imap
c
c  --- add emissions to totals array ----
c
            ivoc = iemvoc - 1 + imap
            emssum(idx,ivoc) = emssum(idx,ivoc) + emspnt(i)
         endif
  70  continue
c
      return
      end
