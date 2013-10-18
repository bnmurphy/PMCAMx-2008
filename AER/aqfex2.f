      subroutine aqfex2(neqa,t,y,f,rpar,ipar)
      include 'aerpar.inc'
      include 'droppar.inc'
      real*4 a(meqn2),b(meqn2),y(meqn2),f(meqn2),rpar
      common/time1/ tnow
c
      tnow=t
c
      call aqrates2(meqn2,y,a,b,f)
c
      return
      end      
