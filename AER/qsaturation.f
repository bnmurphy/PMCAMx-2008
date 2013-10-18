       subroutine qsaturation(temp,qsat)
c
c      THIS SUBROUTINE CALCULATES THE SATURATION MASS CONCENTRATION (in ug/m3)
c      OVER LIQUID WATER FOR A TEMPERATURE TEMP (K)
c
       real*8 rideal,a0,a1,a2,a3,a4,a5,a6,esat,csat
       real*4 temp,qsat
c
       t = temp-273.15                             ! in C
       rideal = 0.08206d0
       a0 = 6.107799961d-0
       a1 = 4.436518521d-1
       a2 = 1.428945805d-2
       a3 = 2.650648471d-4
       a4 = 3.031240396d-6
       a5 = 2.034080948d-8
       a6 = 6.136820929d-11
c
       esat=a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t))))) ! in mb
       psat = esat/(1000.0*1.01325)                    ! in atm
       csat = psat/(rideal*temp)                       ! in mole/l
       qsat = 18000.0d0*csat*1.e6                      ! in ug/m3
c       write(6,*)t,esat/1000.,qsat
       return
       end
