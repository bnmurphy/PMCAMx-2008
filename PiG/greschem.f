      subroutine greschem(dt,ldark,bdnlo3,bdnlno,h2o,pk,
     &                    convfac,conpig,o3dif)
c
c-----CAMx v4.02 030709
c
c     GRESCHEM updates NOx and O3 concentrations in the puff
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c
c     Modifications:
c        07/05/02   gwilson    Changed name of routine to accommodate IRON-PIG
c
c     Input arguments:
c        dt                  time step (s)
c        ldark               darkness flag (logical)
c        bdnlo3              lower bound ozone concentration (ppm)
c        bdnlno              lower bound NO concentration (ppm)
c        h2o                 water vapor concentration (ppm)
c        pk(1)               rate constant for   NO2 + O3 -> NO3
c        pk(2)               rate constant for         O3 -> O(1D)
c        pk(3)               rate constant for      O(1D) -> O(3P)
c        pk(4)               rate constant for      O(1D) -> 2 OH
c        pk(5)               rate constant for  NO3 + NO2 -> NO + NO2
c        pk(6)               rate constant for  NO3 + NO2 -> N2O5
c        pk(7)               rate constant for N2O5 + H2O -> 2 HNO3
c        pk(8)               rate constant for       N2O5 -> NO3 + NO2
c        pk(9)               rate constant for    NO + NO -> 2 NO2
c        convfac             conversion factor from umol/m3 to ppm
c        conpig              concentrations of NO, NO2, HNO3, O3 (umol/m3)
c  
c     Output arguments:
c        conpig              concentrations of NO, NO2, HNO3, O3 (umol/m3)
c        o3dif               ozone reduction due to chemistry (umol/m3)
c
c     Routines Called:
c        none
c
c     Called by:
c        GRESDRIVE
c
      dimension conpig(4), pk(9)
      real no, no2, no3, n2o5
      logical ldark
c
c-----Entry point
c
c-----Convert umol/m3 to ppm
c
      no   = conpig(1)/convfac
      no2  = conpig(2)/convfac
      hno3 = conpig(3)/convfac
      o3   = conpig(4)/convfac
      o3old = o3
      o3dif = 0.
c
c-----NO-NO self reaction
c
      tmpno = no/(1. + pk(9)*dt*no/3600.)
      no2 = no2 + (no - tmpno)
      no = tmpno 
c
c-----NO-O3 titration
c
      if (no.le.bdnlno .or. o3.le.bdnlo3) goto 10
      if (o3.le.no) then
        tmp = o3 - bdnlo3
        o3 = bdnlo3
        no2 = no2 + tmp
        no = amax1((no-tmp),bdnlno)
      else
        tmp = no - bdnlno
        no = bdnlno
        no2 = no2 + tmp
        o3 = amax1((o3-tmp),bdnlo3)
      endif
 10   o3dif = o3 - o3old
c
c-----HNO3 production
c
      if (o3.gt.bdnlo3) then
        if (ldark) then
c
c     Nighttime: Generate nitric acid via the N2O5+H2O channel;
c                after calculating a steady-state N2O5 concentration,
c                assume the process is limited by ozone availability 
c                and use this rate to calculate the decay of ozone;
c                note that loss of one N2O5 removes two ozone molecules 
c
          flux1  = pk(1)*no2*o3  
          flux5 = pk(5)*no2  
          flux6 = pk(6)*no2  
          flux7 = pk(7)*h2o  
          flux8 = pk(8) 
c 
          ssratio = flux6/(flux7 + flux8)  
          no3 = flux1/(flux5 + flux6 - ssratio*flux8 + 1.e-20)  
c  
          n2o5 = no3*ssratio  
          flux7 = n2o5*flux7 
c
          do3 = o3 * (1.0 - exp(-flux7*dt*2.0/(o3*3600.))) 
          do3 = amin1(o3,do3,no2)
          o3 = o3 - do3 
          no2 = no2 - do3 
          hno3 = hno3 + do3 
          o3dif = o3 - o3old
        else
c
c     Daytime: ozone photolysis generates OH radicals, which oxidize NO2
c              to nitric acid
c
          alpha = (dt/3600.)*pk(2)*pk(4)*h2o/(pk(3) + pk(4)*h2o)
          do3 = amin1(o3,alpha*o3,no2/2.)
          o3 = o3 - do3
          no2 = no2 - 2.*do3
          hno3 = hno3 + 2.*do3
          o3dif = o3 - o3old
        endif
      endif
c
c-----Convert from ppm to umol/m3 
c
      conpig(1) = no*convfac
      conpig(2) = no2*convfac
      conpig(3) = hno3*convfac
      conpig(4) = o3*convfac
      o3dif = o3dif*convfac
c
      return
      end
