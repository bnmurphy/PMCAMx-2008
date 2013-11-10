      subroutine spfcn (n,ct,cs,ca,mw,cpx,tom,fval,oaro)
      implicit none
c
c SPFCN calculates the objective function for the bi-section solver in SOAP
c  - revised for new objective function by bkoo (05/27/05)
c     Total Organics in Mole (TOM) = SUM_i(C_i(aer)/MW_i) + C_pre/MW_pre
c     C_i(aer) = C_i(tot) - x_i * Cstar_i
c              = C_i(tot) - (C_i(aer)/MW_i/TOM) * Cstar_i
c  => C_i(aer) = C_i(tot) * TOM / (TOM + Cstar_i/MW_i)
c  => SUM_i(C_i(tot) * TOM / (TOM*MW_i + Cstar_i)) + C_pre/MW_pre - TOM = 0

c BNM - 2-18-09. Added density correction. From now on, use unit density yields
c  		 in rateslo5 file.
c	oaro(n) - density of component. g cm^-3

c
c Called by SOAP
c
      integer     n,i
      real        ct(n),cs(n),ca(n),mw(n),cpx,tom,fval,oaro(n)
c
      fval = 0.0
      do i = 1, n
        ca(i) = ct(i) * oaro(i) * tom / ( tom + cs(i) / mw(i) )
        fval  = fval + ca(i) / mw(i)
      enddo
      fval = fval + cpx - tom
c
      return
      end
