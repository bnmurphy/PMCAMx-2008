      subroutine steady(radius,temp,c,con,gcon,akeq,akhen,akre)
c      
c..INPUTS:
c     RADIUS : DROPLET RADIUS IN M
c     TEMP : TEMPERATURE (IN K)
c     C(46) : THE CONCENTRATIONS OF THE REST OF THE AQUEOUS-PHASE SPECIES
c     GCON(22) : GAS-PHASE CONCENTRATIONS
c     AKEQ,AKHEN,AKRE : REACTION CONSTANTS
c..OUTPUTS:
c     X(8) THE VALUES OF THE STEADY STATE SPECIES CONCENTRATIONS
c
      real*4 kn,n,ikn,kmt
      real*4 temp
ckf      real*4 x(8)
      real*4 x(6)
      real*4 c(46),gcon(22),akeq(17),akhen(21),akre(120)
      real*4 con(28)
c
c     AIRL IS THE MEAN FREE PATH OF AIR. LATER WE HAVE TO SUBSTITUTE
c     THE NUMERICAL VALUE GIVEN HERE BY A FUNCTION OF TEMPERATURE
c     AIRL=65x10-9  m
      airl=65.e-9
c     KN IS THE KNUDSEN NUMBER
      kn=airl/radius
      ikn=1.0/kn
c     ACC IS THE ACCOMODATION COEFFICIENT ASSUMED THE SAME HERE FOR
c     ALL THE SPECIES
      acc=0.01
c     N IS THE COEFFICIENT ENTERING THE FLUX EXPRESSION
      n=1.0/(1.+((1.33+0.71*ikn)/(1.+ikn)+4.*(1.-acc)
     &/(3.*acc))*kn)
c     DG IS THE GAS PHASE DIFFUSIVITY ASSUMED HERE THE SAME FOR ALL
c     THE GASES. DG=1.x10-5 m**2/sec
      dg=1.0e-5
c     RIDEAL IS THE GAS CONSTANT =0.082 (LT.ATM/MOL K)
      rideal=0.082058
      kmt=(3.0*n*dg)/(radius*radius)
c
c     ITERATION LOOP
c
      do icount=1,2

c
c     NO3 CALCULATION
c      
ckf      x(3)=(kmt*gcon(18))/(akre(43)*c(8)+akre(45)+akre(46)*c(29)+
ckf     &akre(47)*c(30) +akre(48)*c(14)+
ckf     &akre(49)*c(27)+akre(54)*c(18)+akre(59)*c(19)+akre(71)*c(35)+
ckf     &akre(109)*c(2)+kmt/(akhen(18)*rideal*temp))
ckf      con(18)=x(3)
c
c     SO4-  CALCULATION
c      
ckf      x(5)=(akre(109)*c(2)*x(3)+2.*akre(86)*c(40)*c(40))
      x(4)=(akre(109)*c(2)*gcon(18)+2.*akre(86)*c(40)*c(40))
     & /(akre(89)*c(41)+akre(92)*c(2)+
     &akre(93)*c(3)+akre(94)*c(29)+akre(95)*c(30)+
     &akre(96)*c(45)+akre(97)*c(14)+akre(98)*c(8)+
     &akre(99)*c(12)+akre(100)*c(19)+akre(101)*c(27)+
     &akre(102)*c(18)+akre(108)*c(35))
ckf      c(39)=x(5)
      c(39)=x(4)
      con(24)=c(39)
c
      a1=c(46)/(akeq(15)+c(46))
      a2=akeq(15)/(akeq(15)+c(46))
c
c     HO2 CALCULATION
c      
      x(2)=((akre(48)*c(14)+akre(54)*c(18)+akre(59)*c(19)+
ckf     &      akre(71)*c(35)) * x(3)+
     &      akre(71)*c(35)) * gcon(18)+
     &  (akre(97)*c(14)+akre(100)*c(19)+akre(102)*c(18)+
ckf     &  akre(108)*c(35))*x(5)+
     &  akre(108)*c(35))*x(4)+
     & 2.0*akre(14)*c(45)*c(22)+
     &akre(28)*c(14)*c(36)+akre(29)*c(14)*c(37)+akre(55)*c(18)*c(22)+
     &akre(56)*c(18)*c(36)+akre(61)*c(19)*c(36)+akre(69)*c(35)*c(36)+
     &akre(65)*c(25)+akre(15)*c(22)*c(15)+akre(58)*c(19)*c(22)+
     &kmt*gcon(17) +akre(5)*c(14)*c(28)+akre(11)*c(22)*c(28)+
     &akre(20)*c(14)*c(44)+akre(50)*c(17)*c(28)+akre(52)*c(18)*c(28)+
     &akre(57)*c(19)*c(28)+akre(60)*c(19)*c(44)+akre(67)*c(35)*c(28)+
     &akre(68)*c(35)*c(44)+akre(84)*c(18)*c(40)+
     &akre(85)*c(19)*c(40))/
     &(a1*(akre(3)*c(28)+2.*akre(6)*c(29)+2.*akre(7)*c(30)+
     &akre(9)*c(14)+akre(12)*c(22)+akre(25)*c(36)+
ckf     &akre(27)*c(37)+akre(46)*x(3)+akre(63)*c(34)+
     &akre(27)*c(37)+akre(46)*gcon(18)+akre(63)*c(34)+
     &akre(94)*c(39)+akre(107)*c(3))+
     &a2*(akre(4)*c(28)+2.*akre(8)*c(30)+akre(10)*c(14)+
     &akre(13)*c(22)+akre(18)*c(12)+akre(19)*c(44)+akre(26)*c(36)+
ckf     &akre(47)*x(3)+akre(64)*c(34)+akre(83)*c(40)+akre(95)*c(39))
     &akre(47)*gcon(18)+akre(64)*c(34)+akre(83)*c(40)+akre(95)*c(39))
     &+(kmt*a1)/(akhen(17)*rideal*temp))
c
      ho2=(x(2)*c(46))/(akeq(15)+c(46))
      o2=(x(2)*akeq(15))/(akeq(15)+c(46))
      c(29)=ho2
      c(30)=o2
      con(17)=c(29)+c(30)
c
      a3=(akre(21)*akre(22)*c(27))/(akre(22)+akre(23)*c(46))
      a4=(akre(22)*akre(24)*c(37))/(akre(22)+akre(23)*c(46))
c
c     OH CALCULATION
c      
      x(1)=(2.*akre(1)*c(14)+akre(15)*c(15)*c(22)+akre(30)*c(45)*
     &c(36)+
     &akre(35)*c(7)+akre(36)*c(8)+akre(44)*c(10)+akre(55)*c(18)*c(22)+
     &akre(58)*c(19)*c(22)+akre(65)*c(25)+kmt*gcon(16)+a4+
     &(akre(9)*c(14)+akre(12)*c(22)+akre(107)*c(3))*ho2+
ckf     &(akre(10)*c(14)+akre(13)*c(22))*o2+akre(96)*c(45)*x(5))/
     &(akre(10)*c(14)+akre(13)*c(22))*o2+akre(96)*c(45)*x(4))/
     &(akre(3)*ho2+akre(4)*o2+akre(5)*c(14)+akre(11)*c(22)+
     &akre(17)*c(12)+akre(21)*c(27)+akre(33)*c(20)+akre(34)*c(21)+
     &akre(37)*c(7)+akre(38)*c(8)+akre(50)*c(17)+akre(52)*c(18)+
     &akre(57)*c(19)+akre(66)*c(25)+akre(67)*c(35)+akre(80)*c(3)+
     &akre(81)*c(2)+akre(88)*c(41)+akre(115)*c(42)+
     &kmt/(akhen(16)*rideal*temp)-a3)
      c(28)=x(1)
      con(16)=c(28)
c
c     ClOH- CALCULATION
c      
ckf      x(4)=(akre(21)*c(27)*x(1)+akre(24)*c(37))/(akre(22)+
      x(3)=(akre(21)*c(27)*x(1)+akre(24)*c(37))/(akre(22)+
     &akre(23)*c(46))
ckf      c(38)=x(4)
      c(38)=x(3)
      con(23)=c(38)
c
c     CO3- CALCULATION
c      
ckf      x(7)=(akre(17)*c(12)*x(1)+akre(99)*c(12)*x(5)+akre(18)*c(12)*o2)/
      x(6)=(akre(17)*c(12)*x(1)+akre(99)*c(12)*x(4)+akre(18)*c(12)*o2)/
     &(akre(19)*o2+akre(20)*c(14)+akre(41)*c(8)+akre(60)*c(19)+
     &akre(68)*c(35))
ckf      c(44)=x(7)
      c(44)=x(6)
      con(28)=c(44)
c
c     SO5- CALCULATION
c      
ckf      x(6)=(akre(116)*c(2)*c(36)+(akre(80)*c(3)+akre(81)*c(2)+
      x(5)=(akre(116)*c(2)*c(36)+(akre(80)*c(3)+akre(81)*c(2)+
     &akre(88)*c(41)+akre(115)*c(42))*x(1)+(akre(89)*c(41)+
ckf     &akre(92)*c(2)+akre(93)*c(3))*x(5))/
     &akre(92)*c(2)+akre(93)*c(3))*x(4))/
     &(akre(83)*o2+akre(84)*c(18)+akre(85)*c(19)+2.*akre(86)*c(40))
ckf      c(40)=x(6)
      c(40)=x(5)
      con(25)=c(40)

      enddo
c
      
      
      return
      end
