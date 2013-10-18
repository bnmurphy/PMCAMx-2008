       subroutine electro(x,f,con,spres,cmet,akeq,akhen,wv,temp)
c       
       implicit real*8(a-h,o-z)
       real*8 con(28),spres(21),cmet(4),akeq(17),akhen(21),cc(46)
       real*8 bparam,cparam,diak,cl,hcl
c       
       cc(2)=(akeq(1)*con(1)*x)/(x*x+akeq(1)*x+akeq(1)*akeq(2))       ! HSO3-
       cc(3)=(akeq(1)*akeq(2)*con(1))/(x*x+akeq(1)*x+akeq(1)*akeq(2)) ! SO3--
       cc(5)=(akeq(3)*con(2)*x)/(x*x+akeq(3)*x+akeq(3)*akeq(4))       ! HSO4-
       cc(6)=(akeq(3)*akeq(4)*con(2))/(x*x+akeq(3)*x+akeq(3)*akeq(4)) ! SO4--
c  
c      ** NO2- CALCULATION FROM EQUILIBRIUM **
       dfac=8.314e-2*temp*wv*akhen(3)*(1.+akeq(7)/x)
       hno2 = spres(3)/(1.+dfac)                        ! New HNO2(g) in ppm
       cc(8)=akhen(3)*1.e-6*(akeq(7)/x)*hno2
c
       cc(10)=(akeq(6)*con(4))/(x+akeq(6))                            ! NO3-
c  
c      ** CO2 EQUILIBRIUM (CONSTANT GAS CO2 CONCENTRATION) **
       cc(12)= akeq(8)*akhen(5)*spres(5)*1.e-6/x
       cc(13)= akeq(9)*cc(12)/x
c       
       cc(15)=(akeq(5)*con(6))/(x+akeq(5))                            ! HO2-
c  
c      ** HCOO- EQUILIBRIUM (PARTITIONING WITH GAS-PHASE) **
       dform=8.314e-2*temp*wv*akhen(8)*(1.+akeq(13)/x)
       form=spres(8)/(1.+dform)                          ! New HCOOH
       cc(19)=akhen(8)*1.e-6*(akeq(13)/x)*form
c       
       cc(30)=(akeq(15)*con(17))/(x+akeq(15))                         ! O2-
       cc(38)=con(23)                                                 ! ClOH-
       cc(39)=con(24)                                                 ! SO4-
       cc(40)=con(25)                                                 ! SO5-
       cc(41)=con(26)                                                 ! HSO5-
       cc(42)=(x*con(27))/(x+akeq(17))                                ! HOCH2SO3-
       cc(43)=(akeq(17)*con(27))/(x+akeq(17))                         ! -OCH2SO3-
       cc(44)=con(28)                                                 ! CO3-
       cc(45)=akeq(11)/x                                              ! OH-

         bparam=akeq(16)+con(15)-con(22)
         cparam=akeq(16)*con(22)
         diak=bparam*bparam+4.0*cparam
         if (diak .le. 0.) diak=1.0e-20
         cl=(-bparam+(diak)**0.5)/2.0
         hcl=(x*(con(15)-con(22)+cl))/(x+akeq(14))
         cc(27)=(akeq(14)*hcl)/x                                      ! Cl-
         cc(36)=(akeq(14)*cl*hcl)/(akeq(16)*x)                        ! Cl2-

        cc(33)=(akeq(10)*x*con(19))/(akeq(11)+akeq(10)*x)             ! NH4+
        cc(46)=x                                                      ! H+

        f1=cc(2)+2.0*cc(3)+cc(5)+2.0*cc(6)+cc(8)+cc(10)
        f2=cc(12)+2.0*cc(13)+cc(15)+cc(19)+cc(27)+cc(30)
        f3=cc(36)+cc(38)+cc(39)+cc(40)+cc(41)+cc(42)
        f4=2.0*cc(43)+cc(44)+cc(45)-cc(33)-cc(46)
        f5=-3.0*cmet(1)-2.0*cmet(2)-cmet(3)-2.0*cmet(4)

        f=f1+f2+f3+f4+f5

        return
        end
