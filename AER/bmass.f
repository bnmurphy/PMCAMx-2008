      subroutine bmass(wl,radius,temp,gcon,con,c,akeq,akhen,fgl,flg,
     &gfgl,gflg)
      double precision kn,n,ikn,kmt
         real*4 temp
      dimension gcon(22),con(28),akeq(17),akhen(21),fgl(21),flg(21)
      dimension c(46),gfgl(21),gflg(21)
c     EKHEN(I) IS THE EFFECTIVE HENRY'S LAW CONSTANT
      dimension ekhen(21)

      DO 1360 I=1,21
      EKHEN(I)=0.0d0
      FGL(I)=0.0d0
      FLG(I)=0.0d0
      GFGL(I)=0.0d0
      GFLG(I)=0.0d0
1360  CONTINUE

      ekhen(16)=akhen(16)
      ekhen(17)=akhen(17)*(1.d0+akeq(15)/c(46))
      ekhen(18)=akhen(18)

c     AIRL IS THE MEAN FREE PATH OF AIR. LATER WE HAVE TO SUBSTITUTE
c     THE NUMERICAL VALUE GIVEN HERE BY A FUNCTION OF TEMPERATURE
c     AIRL=65x10-9  m
      airl=65.d-9
c     KN IS THE KNUDSEN NUMBER
      kn=airl/radius
      ikn=1.d0/kn
c     ACC IS THE ACCOMODATION COEFFICIENT ASSUMED THE SAME HERE FOR
c     ALL THE SPECIES
      acc=0.1D0
c     N IS THE COEFFICIENT ENTERING THE FLUX EXPRESSION
      n=1.d0/(1.d0+((1.33d0+0.71d0*ikn)/(1.d0+ikn)+4.d0*(1.d0-acc)
     &/(3.d0*acc))*kn)
c     DG IS THE GAS PHASE DIFFUSIVITY ASSUMED HERE THE SAME FOR ALL
c     THE GASES.WE SHALL PROBABLY HAVE TO CHANGE IT LATER.
c     DG=1.x10-5 m**2/sec
      dg=1.0d-5
c     RIDEAL IS THE GAS CONSTANT =0.082 (LT.ATM/MOL K)
      rideal=0.082058d0
      kmt=(3.0d0*n*dg)/(radius*radius)

      do 1500 i=16,18
      fgl(i)=kmt*gcon(i)
      flg(i)=(kmt*con(i))/(ekhen(i)*rideal*temp)
      gfgl(i)=fgl(i)*wl
      gflg(i)=flg(i)*wl
1500  continue

      return
      end
