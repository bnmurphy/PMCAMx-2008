      subroutine exptbl(rxntyp,rxnord,rxnpar)
c 
c-----CAMx v4.02 030709
c 
c     EXPTBL calculates a lookup table for temperature/pressure-dependent
c     rate constants using parameters from the chemistry parameters file.
c     The pressure and temperature dependance of the conversion factor from
c     cm3 molec-1 units to ppm units is accounted for.
c     The expression types supported are:
c       Type 1 Temperature independent
c       Type 2 UAM/OZIPM format Arrhenius
c       Type 3 Generalized temperature dependent
c       Type 4 Troe teperature/pressure dependent
c       Type 5 Equilibrium constant ratio to another reaction
c       Type 6 Lindemann-Hinshelwood used for OH + HNO3
c       Type 7 k = k1 + k2[M] used for OH + CO
c     During the model run rate constants are always interpolated from 
c     the lookup table.
c 
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c           
c     Modifications: 
c        none 
c 
c     Input arguments: 
c        rxntyp              rate constant expression type
c        rxnord              order of reaction
c        rxnpar              rate constant expression parameters
c 
c     Output arguments: 
c        none
c             
c     Routines Called:
c        none
c             
c     Called by: 
c        CHMPREP
c
      include "camx.prm"
      include "chmstry.com"
c
      integer rxntyp(MXRXN), rxnord(MXRXN)
      real rxnpar(MXRXN,12)
c
c-----Entry point
c
      dtemp = (TEMPHI-TEMPLO)/(NTEMPR - 1)
      dpres = (PRESHI-PRESLO)/(NPRESR - 1)
      do j = 1,NTEMPR
        tempr(j) = TEMPLO + j*dtemp
      enddo
      do k = 1,NPRESR
        presr(k) = PRESLO + k*dpres
      enddo
c
c-----Fill lookup table
c
      do j = 1,NTEMPR
        temp = tempr(j)
        do k = 1,NPRESR
          ppmpres = 1.e6*presr(k)/1013.
          factor = (298.0/temp)*(presr(k)/1013.)
          do i = 1,nreact
            if (rxntyp(i).eq.1) then
              rktbl(i,j,k) = rxnpar(i,1)*60.0
            elseif (rxntyp(i).eq.2) then
              exptmp = rxnpar(i,2)*(1./298. - 1./temp)
              rktbl(i,j,k) = rxnpar(i,1)*exp(exptmp)*60.0
            elseif (rxntyp(i).eq.3) then
              rktbl(i,j,k) = arren(rxnpar(i,1),rxnpar(i,2),
     &                       rxnpar(i,3),rxnpar(i,4),temp)*60.0
            elseif (rxntyp(i).eq.4) then
              rk0 = arren(rxnpar(i,1),rxnpar(i,2),
     &              rxnpar(i,3),rxnpar(i,4),temp)*ppmpres
              ratio = rk0/arren(rxnpar(i,5),rxnpar(i,6),
     &              rxnpar(i,7),rxnpar(i,8),temp)
              rktbl(i,j,k) = 
     &          (rk0/(1.0+ratio))*rxnpar(i,9)**(1.0/(1.0+
     &          (0.43429*alog(ratio)/rxnpar(i,10))**2))*60.0
            elseif (rxntyp(i).eq.5) then
              iref = anint(rxnpar(i,1))
              rktmp = arren(rxnpar(i,2),rxnpar(i,3),
     &                rxnpar(i,4),rxnpar(i,5),temp)
              rktbl(i,j,k) = rktbl(iref,j,k) / rktmp 
            elseif (rxntyp(i).eq.6) then
              rk1 = arren(rxnpar(i,1),rxnpar(i,2),
     &              rxnpar(i,3),rxnpar(i,4),temp)
              rk2 = arren(rxnpar(i,5),rxnpar(i,6),
     &              rxnpar(i,7),rxnpar(i,8),temp)
              rk3 = arren(rxnpar(i,9),rxnpar(i,10),
     &              rxnpar(i,11),rxnpar(i,12),temp)
              rktbl(i,j,k) = ( rk1 + (rk3*ppmpres /
     &              (1.0 + (rk3*ppmpres / rk2) ) ) )*60.0
            elseif (rxntyp(i).eq.7) then
              rk1 = arren(rxnpar(i,1),rxnpar(i,2),
     &              rxnpar(i,3),rxnpar(i,4),temp)        
              rk2 = arren(rxnpar(i,5),rxnpar(i,6),
     &              rxnpar(i,7),rxnpar(i,8),temp)
              rktbl(i,j,k) = ( rk1 + rk2*ppmpres )*60.0
            endif
            rktbl(i,j,k) = rktbl(i,j,k)*factor**(rxnord(i)-1)
          enddo
        enddo
      enddo
c
      return
      end
c
      function arren(a,ea,b,tref,temp)
c
c --- Calculate a rate constant using k = A*(T/Tref)^B*exp(Ea/T)
c
      real arren
c
      arren = a*((temp/tref)**b)*exp(ea/temp)
c
      end
