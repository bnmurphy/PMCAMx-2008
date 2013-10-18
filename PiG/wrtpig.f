      subroutine wrtpig(idatpig,timpig,idum,rdum)
c
c-----CAMx v4.02 030709
c
c     WRTPIG writes pig parameters and state variables for restart 
c                          
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c          
c     Modifications:
c        07/05/02   gwilson  Added code to accomodate IRON-PiG
c
c     Input arguments:
c        idatpig             model date (YYJJJ)
c        timpig              model time
c        idum                dummy integer variable
c        rdum                dummy real variable
c
c     Output arguments:
c        none
c
c     Subroutines called:
c        none
c
c     Called by:
c        CAMx
c
      include "camx.prm"
      include "filunit.com"
      include "pigsty.com"
      include "flags.com"
      include "chmstry.com"
c
c-----Entry point
c
      if( ipigflg .EQ. GRESPIG ) then
          write(ipig) idatpig,timpig,npig
c
          write(ipig) (ingrd(n),idpig(n),iipig(n),jjpig(n),xpig(n),
     &            ypig(n),
     &            zpig(n),xlength(n),axisy(n),axisz(n),sigy(n),sigz(n),
     &            (puffmass(i,n),i=1,4),fmspig(n),agepig(n),
     &            thetapig(n),n=1,npig)
c
      else if( ipigflg .EQ. IRONPIG ) then
          write(ipig) idatpig,timpig,npig,nreactr
c
          write(ipig) (ingrd(n),idpig(n),iipig(n),jjpig(n),xpig(n),
     &            ypig(n),
     &            zpig(n),xlength(n),axisy(n),axisz(n),sigy(n),sigz(n),
     &            pufftop(n), puffbot(n),
     &            ((puffpol(i,nr,n),i=1,nspec),nr=1,nreactr),
     &            fmspig(n),agepig(n),thetapig(n),n=1,npig)
      endif
c
      return
      end
