      subroutine state(n,x,f,ipar)

c      X  : THE SOLUTION VECTOR
c      F  : VECTOR OF FUNCTIONS THAT MUST BE MADE ZERO
c      N  : NUMBER OF EQUATIONS TO BE SOLVED
c      IPAR: PARAMETER PASSED TO THE USER SUPPLIED SUBROUTINE
c           NOT USED IN THIS CASE
      

ckf       dimension x(7),f(7)
       dimension x(6),f(6)
       dimension c(46),cmet(4),con(28),gcon(22),scon(28),arr(120)
       dimension ar(28),akeq(17),akhen(21),akre(120)
       dimension afgl(21),aflg(21),agfgl(21),agflg(21)
       integer ipar, n
       real*4 temp, pres
       common /aqrate2/ akeq, akhen, akre
       common/sstate/gcon,con,cmet,rad,wl,cph
       common/tprof/temp,pres
c       rideal=0.08206d0
       do 1600 i=1,28
       scon(i)=con(i)
1600   continue

c      WE RENEW THE CON MATRIX USING THE NEW STEADY STATE CONCENTRATIONS
       scon(16)=x(1)
       scon(17)=x(2)
ckf       scon(18)=x(3)
       scon(23)=x(3)
       scon(24)=x(4)
       scon(25)=x(5)
       scon(28)=x(6)

c      WE RENEW THE C MATRIX
c       call equil(scon,cmet,akeq,xsol)
       xsol=cph
       call values(xsol,scon,cmet,akeq,c)

c      WE CALCULATE THE REACTION RATES (arr(i))
       call react(c,cmet,scon,akre,arr,arytm)

c      WE CALCULATE THE PRODUCTION OR LOSS RATES FOR THE 7 STEADY STATE
c      SPECIES. BADDIT IS A SMALL VERSION OF ADDIT GIVING ONLY THE
c      RATES FOR THE 7 SPECIES
       call baddit(arr,ar)
c      THE SAME AS BEFORE FOR MASS TRANFER
       call bmass(wl,rad,temp,gcon,scon,c,akeq,akhen,afgl,aflg,
     & agfgl,agflg)

       f(1)=1.0d13*(ar(16)+afgl(16)-aflg(16))
       f(2)=1.0d13*(ar(17)+afgl(17)-aflg(17))
ckf       f(3)=1.0d13*(ar(18)+afgl(18)-aflg(18))
       f(3)=1.0d13*ar(23)
       f(4)=1.0d13*ar(24)
       f(5)=1.0d13*ar(25)
       f(6)=1.0d13*ar(28)
       return
       end
