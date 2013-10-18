      subroutine intrpdat
c
c-----CAMx v4.02 030709
c
c     INTRPDAT generates fine grid data fields from parent grids, through
c     assignment or linear interpolation.  For wind and temperature fields,
c     also determines time-dependence for each fine grid
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        8/30/02    Modifications for combined cloud/rain file, and to allow
c                   water vapor and cloud/rain files to be read for each nest
c
c     Input arguments:
c        none
c
c     Output arguments:
c        none
c
c     Subroutine called:
c        INTERP2D
c        EXPNDLAY
c        FINWIND
c        INTERPV
c        TIMRATES
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "camx.com"
      include "camxfld.com"
      include "filunit.com"
      include "flags.com"
      include "grid.com"
c
c-----Entry point
c
      do 100 ip = 1,ngrid
        do 99 ic = 1,nchdrn(ip)
          igrd = idchdrn(ic,ip)
c
c-----Vertical diffusion coefficient
c
          iunit = ikv(igrd)
          if (iunit.eq.0) then
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating KV from parent grid',
     &                             time, date,' grid',igrd
            call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    rkv(iptr3d(ip)),rkv(iptr3d(igrd)) )
            call expndlay(ncol(igrd),nrow(igrd),nlay(igrd),nlay(ip),
     &                    nmshv(1,igrd),rkv(iptr3d(igrd)) )
          endif
c
c-----Wind field
c
          iunit = iwind(igrd)
          call finwind(iunit,ncol(ip),nrow(ip),nlay(ip),i1(igrd),
     &                 j1(igrd),nmesh(igrd),nmshv(1,igrd),ncol(igrd),
     &                 nrow(igrd),nlay(igrd),unxt(iptr3d(ip)),
     &                 vnxt(iptr3d(ip)),
     &                 unxt(iptr3d(igrd)),vnxt(iptr3d(igrd)),
     &                 igrd,date,time,iout )
          call timrates(ncol(igrd),nrow(igrd),nlay(igrd),
     &                  windu(iptr3d(igrd)),unxt(iptr3d(igrd)),
     &                  pupt(iptr3d(igrd)) )
          call timrates(ncol(igrd),nrow(igrd),nlay(igrd),
     &                  windv(iptr3d(igrd)),vnxt(iptr3d(igrd)),
     &                  pvpt(iptr3d(igrd)))
c
c-----Temperature
c
          iunit = itemp(igrd)
          if (iunit.eq.0) then
            write(iout,'(a40,f7.0,i8.5,a,i3)')
     &                 'Interpolating temps from parent grid',
     &                             time, date,' grid',igrd
            call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    tnxt(iptr3d(ip)),tnxt(iptr3d(igrd)) )
            call interpv(ncol(ip),nrow(ip),nlay(ip),ncol(igrd),
     &                   nrow(igrd),nlay(igrd),nmesh(igrd),
     &                   nmshv(1,igrd),i1(igrd),j1(igrd),
     &                   hnxt(iptr3d(ip)),
     &                   hnxt(iptr3d(igrd)),tnxt(iptr3d(igrd)) )
            call timrates(ncol(igrd),nrow(igrd),nlay(igrd),
     &                    tempk(iptr3d(igrd)),tnxt(iptr3d(igrd)),
     &                    ptpt(iptr3d(igrd)) )
c
c-----Surface temperature
c
            call interp2d(ncol(ip),nrow(ip),1,i1(igrd),j1(igrd),
     &                    nmesh(igrd),ncol(igrd),nrow(igrd),
     &                    tsnxt(iptr2d(ip)),tsnxt(iptr2d(igrd)) )
            call timrates(ncol(igrd),nrow(igrd),1,
     &                    tsurf(iptr2d(igrd)),tsnxt(iptr2d(igrd)),
     &                    pspt(iptr2d(igrd)) )
          endif
c
c-----Water vapor
c
          if (lchem) then
            iunit = ih2o(igrd)
            if (iunit.eq.0) then
              call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),
     &                      j1(igrd),nmesh(igrd),ncol(igrd),nrow(igrd),
     &                      water(iptr3d(ip)),water(iptr3d(igrd)) )
              call expndlay(ncol(igrd),nrow(igrd),nlay(igrd),nlay(ip),
     &                      nmshv(1,igrd),water(iptr3d(igrd)) )
            endif
          endif
c
c-----Cloud/rain
c
          if (lchem .or. lwet .or. ldry) then
            iunit = icld(igrd)
            if (iunit.eq.0) then
              call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),
     &                      j1(igrd),nmesh(igrd),ncol(igrd),nrow(igrd),
     &                      fcloud(iptr3d(ip)),fcloud(iptr3d(igrd)) )
              call expndlay(ncol(igrd),nrow(igrd),nlay(igrd),nlay(ip),
     &                      nmshv(1,igrd),fcloud(iptr3d(igrd)) )
c
              call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),
     &                      j1(igrd),nmesh(igrd),ncol(igrd),nrow(igrd),
     &                      cldtrns(iptr3d(ip)),cldtrns(iptr3d(igrd)) )
              call expndlay(ncol(igrd),nrow(igrd),nlay(igrd),nlay(ip),
     &                      nmshv(1,igrd),cldtrns(iptr3d(igrd)) )

              call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),
     &                      j1(igrd),nmesh(igrd),ncol(igrd),nrow(igrd),
     &                      cwc(iptr3d(ip)),cwc(iptr3d(igrd)) )
              call expndlay(ncol(igrd),nrow(igrd),nlay(igrd),nlay(ip),
     &                      nmshv(1,igrd),cwc(iptr3d(igrd)) )
c
              call interp2d(ncol(ip),nrow(ip),nlay(ip),i1(igrd),
     &                      j1(igrd),nmesh(igrd),ncol(igrd),nrow(igrd),
     &                      pwc(iptr3d(ip)),pwc(iptr3d(igrd)) )
              call expndlay(ncol(igrd),nrow(igrd),nlay(igrd),nlay(ip),
     &                      nmshv(1,igrd),pwc(iptr3d(igrd)) )
            endif
          endif
c
  99    continue
 100  continue
c
      return
      end
