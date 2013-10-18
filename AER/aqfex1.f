      subroutine aqfex1(neqa,t,y,f,rpar,ipar)
      include 'aerpar.inc'
      include 'droppar.inc'
      real*4 a(meqn1),b(meqn1),y(meqn1),f(meqn1),rpar
      common/time1/ tnow
c
      tnow=t
c
      call aqrates1(meqn1,y,a,b,f)
c
      return
      end
