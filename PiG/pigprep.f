      subroutine pigprep(begdat,begtim)
c
c-----CAMx v4.02 030709
c
c     PIGPREP prepares the PiG submodel
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        6/27/00   Added some diagnostics to ensure that the point source PiG
c                  list is consistent on a restart with that used in previous 
c                  run.
c                - Now writing PiG point source information to PiG restart
c                  file, and echoing this information to diagnostic file on 
c                  model startup.
c       11/06/01   Input dates are now Julian
c        7/05/02   Added code to accommodate IRON-PiG
c
c     Input arguments:
c        begdat              model beginning date (YYJJJ)
c        begtim              model beginning time
c
c     Output arguments:
c        none
c
c     Routines Called:
c        none
c
c     Called by:
c        STARTUP
c
      include "camx.prm"
      include "grid.com"
      include "flags.com"
      include "pigsty.com"
      include "filunit.com"
      include "camxfld.com"
      include "ptemiss.com"
      include "chmstry.com"
c
      integer begdat,nsrc,irec,idpsrc(MXPTSRC)
      real xx(MXPTSRC),yy(MXPTSRC),hh(MXPTSRC),dd(MXPTSRC),tt(MXPTSRC),
     &     vv(MXPTSRC)
      real xsrc(MXPTSRC),ysrc(MXPTSRC),hsrc(MXPTSRC),dsrc(MXPTSRC),
     &     tsrc(MXPTSRC),vsrc(MXPTSRC)
c
c-----Entry point
c
c-----Check FLEAK param
c
      if( ipigflg .EQ. IRONPIG .AND. FLEAK .GT. 1.0) then
         write(iout,'(//,a)') 'ERROR in PIGPREP:'
         write(iout,*) 'FLEAK must be .LE. 1.0'
         write(iout,*) 'FLEAK = ',FLEAK
         call camxerr()
      endif
c
c-----Set overlay failure counter for IRONPIG to zero
c
      kopfail=0
c 
c-----Read PiG information if it is a restart run; check against point 
c     source PiG list to ensure consistency between runs
c
      irec = 0
      if (lrstrt) then
        write(iout,*)
        read(irstp) nsrc
        read(irstp) (idpsrc(n),xsrc(n),ysrc(n),hsrc(n),dsrc(n),tsrc(n),
     &               vsrc(n),n=1,nsrc)
        rewind(iptem)
        do n = 1,5
          read(iptem)
        enddo
        read(iptem) (xx(n),yy(n),hh(n),dd(n),tt(n),vv(n),n=1,nptsrc)
        nn = 0
        do n = 1,nptsrc
          if (dd(n).lt.0.) then
            nn = nn + 1
            if (nn.gt.nsrc) then
              write(iout,'(//,a)') 'ERROR in PIGPREP:'
              write(iout,*) 
     &                'The number of PiG sources on the point source'
              write(iout,*) 
     &                'file exceeds the number of PiG sources on the'
              write(iout,*) 'PiG restart file'
              write(iout,*) 'Point source file: ',nn
              write(iout,*) ' PiG restart file: ',nsrc
              call camxerr()
            endif
            if (idpsrc(nn).ne.n .or. xsrc(nn).ne.xx(n) .or. 
     &          ysrc(nn).ne.yy(n) .or. hsrc(nn).ne.hh(n) .or. 
     &          dsrc(nn).ne.dd(n) .or. tsrc(nn).ne.tt(n) .or. 
     &          vsrc(nn).ne.vv(n) ) then
              write(iout,'(//,a)') 'ERROR in PIGPREP:'
              write(iout,*) 'Point source list on PiG restart file'
              write(iout,*) 'is inconsistent (not in same order) as'
              write(iout,*) 'the PiG sources on point source file'
              write(iout,*)
              write(iout,*) 'Last point source read from point file'
              write(iout,*) 'is PiG source # ',nn
              write(iout,*) '   Point ID: ',n
              write(iout,*) ' x-location: ',xx(n)
              write(iout,*) ' y-location: ',yy(n)
              write(iout,*) ' stack hght: ',hh(n)
              write(iout,*) ' stack diam: ',dd(n)
              write(iout,*) ' stack temp: ',tt(n)
              write(iout,*) ' stack exit: ',vv(n)
              write(iout,*)
              write(iout,*) 
     &              'Should match this source read from PiG file:'
              write(iout,*) '   Point ID: ',idpsrc(nn)
              write(iout,*) ' x-location: ',xsrc(nn)
              write(iout,*) ' y-location: ',ysrc(nn)
              write(iout,*) ' stack hght: ',hsrc(nn)
              write(iout,*) ' stack diam: ',dsrc(nn)
              write(iout,*) ' stack temp: ',tsrc(nn)
              write(iout,*) ' stack exit: ',vsrc(nn)
              call camxerr()
            endif
          endif
        enddo
        if (nn.lt.nsrc) then
          write(iout,'(//,a)') 'ERROR in PIGPREP:'
          write(iout,*) 'The number of PiG sources on the point source'
          write(iout,*) 'file is less than the number of PiG sources on'
          write(iout,*) 'the PiG restart file'
          write(iout,*) 'Point source file: ',nn
          write(iout,*) ' PiG restart file: ',nsrc
          call camxerr()
        endif
c
c-----Read hourly data to current date/time
c
 100    continue
        irec = irec + 1
        if( ipigflg .EQ. GRESPIG ) then
           read(irstp,ERR=7000,END=7001) idatpig,timpig,npig
        else if( ipigflg .EQ. IRONPIG ) then
           read(irstp,ERR=7000,END=7001) idatpig,timpig,npig,nreactr
        endif
        irec = irec + 1
        if( ipigflg .EQ. GRESPIG ) then
           read(irstp,ERR=7002,END=7003) (ingrd(n),idpig(n),iipig(n),
     &            jjpig(n),xpig(n),ypig(n),zpig(n),xlength(n),axisy(n),
     &            axisz(n),sigy(n),sigz(n),(puffmass(i,n),i=1,4),
     &            fmspig(n),agepig(n),thetapig(n),n=1,npig)
        else if( ipigflg .EQ. IRONPIG ) then
           read(irstp,ERR=7002,END=7003) (ingrd(n),idpig(n),iipig(n),
     &            jjpig(n),xpig(n),ypig(n),zpig(n),xlength(n),axisy(n),
     &            axisz(n),sigy(n),sigz(n),pufftop(n),puffbot(n),
     &            ((puffpol(i,nr,n),i=1,nspec),nr=1,nreactr),
     &            fmspig(n),agepig(n),thetapig(n),n=1,npig)
        endif
        if (timpig.ge.2400.) then
          timpig = timpig - 2400.
          idatpig = idatpig + 1
          if( MOD(idatpig,1000) .GT. 365 ) then
             if( MOD(INT(idatpig/1000),4) .EQ. 0 ) then
                if( MOD(idatpig,1000) .EQ. 367 )
     &                         idatpig = (INT(idatpig/1000)+1)*1000 + 1
             else
                idatpig = (INT(idatpig/1000)+1)*1000 + 1
             endif
          endif
        endif
        write(iout,'(a,F10.1,i10.5,i10)') 'Read PiG file at  ',
     &                                timpig,idatpig,npig
        if (idatpig.lt.begdat .or. 
     &     (idatpig.eq.begdat .and. timpig.lt.begtim - 0.01)) goto 100
        if (idatpig.gt.begdat .or. 
     &     (idatpig.eq.begdat .and. timpig.gt.begtim + 0.01)) then
           write(iout,'(//,a)') 'ERROR in PIGPREP:'
           write(iout,*) 'Date or time in PiG file > beginning time'
           write(iout,'(2i10.5)') idatpig, timpig
           call camxerr()
        endif
c
c-----Check first and last puff contents
c
        if( ipigflg .EQ. IRONPIG ) then
           write(iout,*) 'pigprp2...npig= ',npig
           write(iout,*) '          nreactr= ',nreactr
c
           do n=1,npig
              if( ingrd(n) .gt. 0) then
                 write(iout,*) 'pigprp2 restart...first live puff is ',n
                 write(iout,*) 'puffpol in first reactor..',
     &                        ( puffpol(is,1,n),is=1,nspec)
                 go to 10
              endif
           enddo
 10        continue
           write(iout,*) ' last puff is ',npig
           write(iout,*) 'ingrd = ',ingrd(npig)
           write(iout,*) 'puffpol in first reactor is ',
     &                     (puffpol(is,1,npig),is=1,nspec)
c
c-----Fill radical array for puffs
c
           do n=1,npig
            do nr=1,nreactr
              do is=1,MXRADCL
                    cncradp(n,nr,is)=1.e-9
               enddo
            enddo
           enddo
         endif
c
c-----If not a restart, find PiG sources on input point source file
c
      else
c            
c-----Set nreactr from .prm file
c
        if( ipigflg .EQ. IRONPIG ) nreactr = MXRECTR
        npig = 0
        nsrc = 0
        rewind(iptem)
        do n = 1,5
          read(iptem)
        enddo
        read(iptem) (xx(n),yy(n),hh(n),dd(n),tt(n),vv(n),n=1,nptsrc)
        do n = 1,nptsrc
          if (dd(n).lt.0.) then
            nsrc = nsrc + 1
            idpsrc(nsrc) = n
            xsrc(nsrc) = xx(n)
            ysrc(nsrc) = yy(n)
            hsrc(nsrc) = hh(n)
            dsrc(nsrc) = dd(n)
            tsrc(nsrc) = tt(n)
            vsrc(nsrc) = vv(n)
          endif
        enddo
      endif
c
c-----Write source information to header of PiG file
c
      write(ipig) nsrc
      write(ipig) (idpsrc(n),xsrc(n),ysrc(n),hsrc(n),dsrc(n),tsrc(n),
     &             vsrc(n),n=1,nsrc)
c
c-----Echo PiG point source information to diagnostic file
c
      write(idiag,'(//,a)') 'PiG source information'
      write(idiag,'(80a)') ('-',m=1,80)
      write(idiag,'(a,a)') ' Pig Src   Pt Src    Xloc       Yloc   ',
     &                     ' Height   Diameter  Temperature  Velocity'
      write(idiag,'(a,a)') '   No.       No.  (m or deg) (m or deg)  ',
     &                     '(m)       (m)         (K)       (m/hr)'
      write(idiag,'(80a)') ('-',m=1,80)
      do n = 1,nsrc
        write(idiag,'(i5,1x,i10,f10.1,f11.1,f9.1,f9.2,f12.1,f13.0)')
     &      n,idpsrc(n),xsrc(n),ysrc(n),hsrc(n),dsrc(n),tsrc(n),vsrc(n)
      enddo
c
      goto 9999
c
c-----------------------------------------------------------------------
c   Error messages:
c-----------------------------------------------------------------------
c
 7000 continue
      write(iout,'(//,a)') 'ERROR in PIGPREP:'
      write(iout,*) 'Cannot read PiG restart file at record: ',irec
      call camxerr()
c
 7001 continue
      write(iout,'(//,a)') 'ERROR in PIGPREP:'
      write(iout,*) 'Premature end-of-file found in PiG ',
     &              'restart file at record: ',irec
      call camxerr()
c
 7002 continue
      write(iout,'(//,a)') 'ERROR in PIGPREP:'
      write(iout,'(a,i10.5,f10.1)') 
     &      'Cannot read data in PiG restart file at hour: ',
     &      idatpig,timpig
      call camxerr()
c
 7003 continue
      write(iout,'(//,a)') 'ERROR in PIGPREP:'
      write(iout,*) 'ERROR: Premature end-of-file found when ',
     &             'reading data in PiG restart file at record: ',irec
      call camxerr()
c
c-----------------------------------------------------------------------
c   Return point:
c-----------------------------------------------------------------------
c
 9999 continue
      return
      end
