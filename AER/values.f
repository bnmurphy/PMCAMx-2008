        subroutine values(x,con,cmet,akeq,cc)

        real*4 con(28),cmet(4),akeq(17),cc(46)
        real*8 bparam,cparam,diak,cl,hcl
c
c    	species in the aqueous phase mechanism
c    	cc (1 - 46)
c	1.)     SO2*H2O		24.)	CH3C(O)OOH 
c	2.)	HSO3(-)		25.)	CH3OOH
c    	3.)	SO3(2-)		26.)	HCl
c    	4.)	H2SO4		27.)	Cl(-)
c    	5.)	HSO4(-)		28.)	OH
c    	6.)	SO4(2-)		29.)	HO2
c    	7.)	HNO2		30.)	O2(-)
c	8.)	NO2(-)		31.)	NO3
c	9.)	HNO3		32.)	NH4OH
c	10.)	NO3(-)		33.)	NH4(+)
c	11.)	CO2*H2O		34.)	CH3O2
c	12.)	HCO3(-)		35.)	CH3OH
c	13.)	CO3(2-)		36.)	Cl2(-)
c	14.)	H2O2		37.)	Cl
c	15.)	HO2(-)		38.)	ClOH(-)
c	16.)	HCHO		39.)	SO4(-)
c	17.)	H2C(OH)2	40.)	SO5(-)
c	18.)	HCOOH		41.)	HSO5(-)
c	19.)	HCOO(-)		42.)	HOCH2SO3(-)
c	20.)	NO		43.)	OCH2SO3(2-)
c	21.)	NO2 		44.)	CO3(-)
c	22.)	O3		45.)	OH(-)
c	23.)	PAN		46.)	H(+)
c	
c	con(1-28)
c
c	1.)	SO2(g)		15.)	HCl(g)
c	2.)	H2SO4(g)	16.)	OH(g)
c	3.)	HNO2(g)		17.)	HO2(g)
c	4.)	HNO3(g)		18.)	NO3(g)
c	5.)	CO2(g)		19.)	NH3(g)
c	6.)	H2O2(g)		20.)	CH3O2(g)
c	7.)	HCHO(g)		21.)	CH3OH(g)
c	8.)	HCOOH(g)	22.)	Cl2(-), Cl
c	9.)	NO(g)		23.)	ClOH(-)
c	10.)	NO2(g)		24.)	SO4(-)
c	11.)	O3(g)		25.)	SO5(-)
c	12.)	PAN(g)		26.)	HSO5(-)
c	13.)	CH3C(O)OOH(g)	27.)	HOCH2SO3(-),OCH2SO3(2-)
c	14.)	CH3OOH(g)	28.)	CO3(-)
c	
         bparam=akeq(16)+con(15)-con(22)
         cparam=akeq(16)*con(22)
         diak=bparam*bparam+4.0d0*cparam
         if (diak .le. 0.0d0) diak=1.0d-30
         cl=(-bparam+(diak)**0.5d0)/2.0d0
         hcl=(x*(con(15)-con(22)+cl))/(x+akeq(14))

       cc(1)=(con(1)*x*x)/(x*x+akeq(1)*x+akeq(1)*akeq(2))
       cc(2)=(akeq(1)*con(1)*x)/(x*x+akeq(1)*x+akeq(1)*akeq(2))
       cc(3)=(akeq(1)*akeq(2)*con(1))/(x*x+akeq(1)*x+akeq(1)*akeq(2))

       cc(4)=(con(2)*x*x)/(x*x+akeq(3)*x+akeq(3)*akeq(4))
       cc(5)=(akeq(3)*con(2)*x)/(x*x+akeq(3)*x+akeq(3)*akeq(4))
       cc(6)=(akeq(3)*akeq(4)*con(2))/(x*x+akeq(3)*x+akeq(3)*akeq(4))

       cc(7)=(x*con(3))/(x+akeq(7))
       cc(8)=(akeq(7)*con(3))/(x+akeq(7))

       cc(9)=(x*con(4))/(x+akeq(6))
       cc(10)=(akeq(6)*con(4))/(x+akeq(6))

       cc(11)=(x*x*con(5))/(x*x+akeq(8)*x+akeq(8)*akeq(9))
       cc(12)=(akeq(8)*con(5)*x)/(x*x+akeq(8)*x+akeq(8)*akeq(9))
       cc(13)=(akeq(8)*akeq(9)*con(5))/(x*x+akeq(8)*x+akeq(8)*akeq(9))

       cc(14)=(x*con(6))/(x+akeq(5))
       cc(15)=(akeq(5)*con(6))/(x+akeq(5))

       cc(16)=con(7)/(1.0d0+akeq(12))
       cc(17)=(akeq(12)*con(7))/(1.0d0+akeq(12))

       cc(18)=(x*con(8))/(x+akeq(13))
       cc(19)=(akeq(13)*con(8))/(x+akeq(13))

       cc(20)=con(9)

       cc(21)=con(10)

       cc(22)=con(11)

       cc(23)=con(12)

       cc(24)=con(13)

       cc(25)=con(14)

       cc(26)=hcl
       cc(27)=(akeq(14)*hcl)/x

       cc(28)=con(16)

       cc(29)=(x*con(17))/(x+akeq(15))
       cc(30)=(akeq(15)*con(17))/(x+akeq(15))

       cc(31)=con(18)

       cc(32)=(akeq(11)*con(19))/(akeq(11)+akeq(10)*x)
       cc(33)=(akeq(10)*x*con(19))/(akeq(11)+akeq(10)*x)

       cc(34)=con(20)

       cc(35)=con(21)

       cc(36)=(akeq(14)*cl*hcl)/(akeq(16)*x)
       cc(37)=cl

       cc(38)=con(23)
       cc(39)=con(24)
       cc(40)=con(25)
       cc(41)=con(26)
       cc(42)=(x*con(27))/(x+akeq(17))
       cc(43)=(akeq(17)*con(27))/(x+akeq(17))
       cc(44)=con(28)
       cc(45)=akeq(11)/x
       cc(46)=x

       return
       end
