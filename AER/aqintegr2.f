      subroutine aqintegr2(neqa,y,stime,stout,temp,p)
      external aqfex2,aqjex
      include 'aerpar.inc'
      include 'droppar.inc'
      include 'math.inc'

      integer iwork(liw2)
      real*4 y(meqn2), rwork(lrw2), ysav(meqn2)
c
      common/time1/ tnow               ! information passes to aqrates
c  
      tnow = stime
c
      lrw = lrw2                        ! dimension of real work vector
      liw = liw2                        ! dimension of integer work vector
      do i=5, 10                        ! zero all optional inputs
        iwork(i) = 0
        rwork(i) = 0.0
      enddo
      iwork(6)=worki                    ! steps allowed
      rwork(6)=workr                    ! maximum step in seconds
      istate=1
      itask=1
      atol=tola
      rtol=tolr
c
c     SET UP SULFUR BALANCE CHECK
c
      do i=1,meqn2
        ysav(i) = y(i)
      enddo
      savtime = stime
      savtout = stout
      sulfbef = y(11)*32./96. + y(12)*32./113. + y(13)*32./111.
     &        + y(19)*32./96. + y(20)*32./113. + y(21)*32./111.
     &        + (y(2)+(y(6)+y(14))*0.7901)*32./64.
c
c     READY FOR THE CALL TO SVODE
c
      call svode(aqfex2,neqa,y,stime,stout,itol,rtol,atol,itask,
     &           istate,iopt,rwork,lrw,iwork,liw,aqjex,mf,
     &           rpar,ipar)    
c
c     CHECK SULFUR BALANCE
c
      sulfaft = y(11)*32./96. + y(12)*32./113. + y(13)*32./111.
     &        + y(19)*32./96. + y(20)*32./113. + y(21)*32./111.
     &        + (y(2)+(y(6)+y(14))*0.7901)*32./64.
      sbal = sulfaft/sulfbef
c
      if (sulfbef.gt.0.1 .and.
     &    (sbal.lt.0.999 .or. sbal.gt.1.001) ) then
        do i=1,meqn2
          y(i) = ysav(i)
        enddo
        stime = savtime
        stout = savtout
        istate = 1
        rtol = tolr*10.
        atol = tola*10.
        call slsode(aqfex2,neqa,y,stime,stout,itol,rtol,atol,itask,
     &             istate,iopt,rwork,lrw,iwork,liw,aqjex,mf)
        sulfaft = y(11)*32./96. + y(12)*32./113. + y(13)*32./111.
     &        + y(19)*32./96. + y(20)*32./113. + y(21)*32./111.
     &        + (y(2)+(y(6)+y(14))*0.7901)*32./64.
        sbal = sulfaft/sulfbef
      endif
      do i=1, meqn2
        y(i) = amax1(y(i), 1.0e-20)
      enddo
c
      return
      end
