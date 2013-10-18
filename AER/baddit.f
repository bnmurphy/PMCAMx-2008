       subroutine baddit(rr,r)
c      PADD IS A WORK MATRIX
       dimension rr(120),r(28),padd(44)
       DO 1310 I=1,28
       R(I)=0.0d0
1310   CONTINUE
       padd(28)=rr(1)-rr(3)-rr(4)-rr(5)+rr(9)+rr(10)-rr(11)+rr(12)+
     & rr(1)+
     & rr(13)+rr(15)-rr(17)-rr(21)+rr(22)+rr(30)-rr(33)-rr(34)+rr(35)+
     & rr(36)-rr(37)-rr(38)+rr(44)-rr(50)-rr(52)+rr(55)-rr(57)+rr(58)+
     & rr(65)-rr(66)-rr(67)-rr(75)-rr(76)-rr(83)+rr(91)+rr(101)-rr(108)
       r(16)=padd(28)

       padd(29)=-rr(3)+rr(5)-rr(6)-rr(7)-rr(9)+rr(11)-rr(12)+rr(14)
     & -rr(6)+
     & rr(20)-rr(25)-rr(27)+rr(28)+rr(29)-rr(46)+rr(48)+rr(50)+rr(52)+
     & rr(54)+rr(55)+rr(56)+rr(57)+rr(59)+rr(60)+rr(61)-rr(63)+rr(65)+
     & rr(67)+rr(68)+rr(69)+rr(71)+rr(79)-rr(89)+rr(92)+rr(95)+rr(97)-
     & rr(101)+rr(102)
       padd(30)=-rr(4)-rr(7)-2.d0*rr(8)-rr(10)-rr(13)+rr(14)+rr(15)-
     & rr(18)-rr(19)-rr(26)-rr(47)+rr(58)-rr(64)-rr(78)+rr(80)-rr(90)
       r(17)=padd(29)+padd(30)

       padd(31)=-rr(43)-rr(45)-rr(46)-rr(47)-rr(48)-rr(49)-rr(54)-
     & rr(59)-rr(71)-rr(103)
       r(18)=padd(31)


       padd(38)=rr(21)-rr(22)-rr(23)+rr(24)
       r(23)=padd(38)

       padd(39)=2.d0*rr(81)-rr(84)-rr(87)-rr(88)-rr(89)-rr(90)-rr(91)-
     & rr(92)-rr(93)-rr(94)-rr(95)-rr(96)-rr(97)-rr(102)+rr(103)
       r(24)=padd(39)

       padd(40)=rr(75)+rr(76)-rr(78)-rr(79)-rr(80)-2.d0*rr(81)+rr(83)+
     & rr(84)+rr(87)+rr(88)+rr(108)+rr(109)
       r(25)=padd(40)


       padd(44)=rr(17)+rr(18)-rr(19)-rr(20)-rr(41)-rr(60)-rr(68)+rr(94)
       r(28)=padd(44)

       return
       end
