      subroutine soap(ntot,caer,cgas,tempk,convfac,
     &                iout,igrdchm,ichm,jchm,kchm,lppm,
     &                cpre,mwpre,csatT)
c BNM-122308      implicit none
c
c**********************************************************************
c                                                                     * 
c                                                                     * 
c  SOAP - A SUBROUTINE TO CALCULATE GAS-AEROSOL PARTITIONING OF       *
c                     SECONDARY ORGANIC SPECIES                       *
c                                                                     *
c                           DEVELOPED BY:                             *
c                                                                     * 
c                           ROSS STRADER                              *
c            DEPARTMENT OF CIVIL/ENVIRONMENTAL ENGINEERING            *
c             CARNEGIE MELLON UNIVERSITY, 5000 FORBES AVE             *
c                   PITTSBURGH, PENNSYLVANIA, 15213                   *
c                                                                     * 
c                               AND                                   *
c                                                                     * 
c                         SPYROS N. PANDIS                            *
c       DEPARTMENTS OF CHEMICAL ENGINEERING AND ENGINEERING AND       *
c                          PUBLIC POLICY                              *
c             CARNEGIE MELLON UNIVERSITY, 5000 FORBES AVE             *
c                   PITTSBURGH, PENNSYLVANIA, 15213                   *
c                                                                     * 
c                                                                     * 
c                                                                     * 
c                                                                     * 
c      ASSOCIATED SUBROUTINES:                                        *
c      CALC       - CALCULATES THE MOLE FRACTION OF EACH SOLUTION-    *
c                   FORMING COMPOUND IN THE PSEUDO-IDEAL SOLUTION     * 
c                   (AEROSOL PHASE)                                   *
c      FCN        - NONLINEAR SYSTEM TO BE SOLVED BY HYBRD            *
c      RNDNUM     - GENERATES RANDOM NUMBERS                          *
c      HYBRD      - SOLVES A SYSTEM OF NONLINEAR EQUATIONS            *
c                                                                     * 
c                                                                     * 
c**********************************************************************
c                                                                     * 
c                                                                     * 
c     THIS MODEL USES A PSEUDO-IDEAL SOLUTION THEORY TO PARTITION     *
c     SECONDARY ORGANIC SPECIES BETWEEN THE AEROSOL AND GAS PHASE.    *
c     THE THEORY IS OUTLINED IN:                                      * 
c                                                                     * 
c     STRADER, R. (1998), 'EVALUATION OF SECONDARY ORGANIC AEROSOL    *
c          FORMATION IN WINTER,' M.S. THESIS, DEPT. OF CIVIL/         *
c          ENVIRONMENTAL ENGINEERING, CARNEGIE MELLON UNIVERSITY      *
c                                                                     * 
c                                                                     * 
c**********************************************************************
c
c
c  MODIFICATIONS
c
c     Revised for new bi-section solver by bkoo (05/27/05)
c     Revised for new SOAP driver by bkoo (03/09/03)
c     Modified by bkoo (01/20/03)
c     Revised for CAMx by Greg Yarwood, January 2003
c
c DEFINITION OF VARIABLES:
c
c  INPUTS
c
c     ntot    - total number of CG/SOA species pairs
c     caer    - aerosol-phase concentrations of SOA species (ug/m3)
c     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
c     tempk   - temperature (K)
c     convfac - conversion factor: umol/m3 = ppm * convfac 
c     iout    - standard output file unit
c     igrdchm - index for grid containing the grid cell
c     ichm    - i index of grid cell
c     jchm    - j index of grid cell
c     kchm    - k index of grid cell
c     lppm    - gases in ppm if true, ug/m3 if false
c     cpre    - concentration of pre-existing organic aerosol (ug/m3)
c     mwpre   - molecular weight of pre-existing organic aerosol (g/mol)
c
c   OUTPUTS
c
c     caer    - aerosol-phase concentrations of SOA species (ug/m3)
c     cgas    - gas-phase concentrations of CG species (ppm or ug/m3)
c     csatT   - saturation concentrations of CG/SOA species at current T
c               (ug/m3)
c
c   VARIABLES USED WITHIN SOAP
c
c     i       - counter
c     icont   - counter
c     sum     - counter
c     nsol    - total number of solution-forming SOA species
c     cstemp  - temperatures corresponding to saturation concentrations
c               of CG/SOA species (K)
c     csat    - saturation concentrations of CG/SOA species (ug/m3)
c     deltah  - enthalpy of vaporization of CG/SOA species (J/mol)
c     flagsoap- set to 1 if CG/SOA species forms solutions; 0 if not
c     mwsoap  - molecular weights of CG/SOA species (g/mol)
c     scaer   - aerosol-phase concentrations of solution-forming 
c               SOA species (ug/m3)
c     scgas   - gas-phase concentrations of solution-forming 
c               SOA species (ug/m3)
c     scsat   - saturation concentrations of solution-forming 
c               SOA species (ug/m3)
c     sctot   - total concentrations of solution-forming SOA species
c               (ug/m3)
c     smw     - molecular weights of solution-forming SOA species 
c               (g/mol)
c     znum    - counter for number of iterations
c     conmin  - use simple solution for species below this level (ug/m3)
c     cpremin - no pre-existing organics if cpre < cpremin (ug/m3) 
c     xtol    - error tolerance for bi-section method
c
c***********************************************************************
c
c VARIABLE DECLARATION
c
      include 'soap.com'
      real        conmin, cpremin, xtol
      parameter ( conmin  = 1.e-6 )
      parameter ( cpremin = 1.01e-9 )
      parameter ( xtol    = 5.0e-5 )
c
      integer     ntot
      real        caer(ntot), cgas(ntot), ctot(NSOAP), csatT(NSOAP)
      real        smw(NSOAP), scsat(NSOAP)
      real        sctot(NSOAP), scaer(NSOAP), scgas(NSOAP)
      integer     idx(NSOAP)
c
      real        mwpre, cpre, tempk, sum, convfac
      integer     iout, igrdchm, ichm, jchm, kchm
      integer     i, icont, nsol, znum
      logical     lppm

      real        cpx, bb, cc, xend, fend, xmid, fmid, dx
c BNM - declaring variables
      dimension tarray2(2)

c      Declare Variables for Hvap look-up table
      integer tempcnt, icstar, ctemp, iter
      real csat1, csat2

c BNM - done declaring variables

c
c***********************************************************************
c
c Entry point
c
      if(ntot.ne.NSOAP) then 
        write(iout,'(//,A)') 'ERROR in SOAP:'
        write(iout,*) 'ERROR: Parameter mismatch in subroutine SOAP.'
        write(iout,'(1X,A,I2)') 'Parameter in include file: NSOAP = ',
     &                                                           NSOAP
        write(iout,'(1X,A,I2)') 'Variable in argument list: ntot = ',
     &                                                            ntot
        write(iout,*) 'Please set the parameter to be consistent.'
        call camxerr()
      endif
      do i=1,ntot
        if (lppm) cgas(i) = cgas(i)*convfac*mwsoap(i) ! convert to ug/m3
        ctot(i) = caer(i) + cgas(i)
      enddo
c
c Check PFLAG to set pre-existing organic mass
c
      if(PFLAG.eq.0) cpre = 0.0
      cpx = cpre/mwpre
c
c CHANGE SATURATION CONCENTRATIONS ACCORDING TO CURRENT TEMPERATURE
c

c BNM - timing
c	tcpu = dtime(tarray2)
c	print *,'Begining time for Hvap = ',tarray2(1)
c BNM

ccc Traditional Hvap calculation, Clausius Clapeyron ccc
cBNM      do i=1,ntot
cBNM         csatT(i)=csat(i)*(cstemp(i)/tempk)*exp((deltah(i)/8.314)
cBNM     &                *(1/cstemp(i)-1/tempk))
cBNM      enddo     

ccccccc  ------- New Hvap Look-Up Table -------- cccccccccc
c      print *,'Makes it to Hvap calc', igrdchm

ccc Using Look-up Table Compiled by Scott Epstein  ccc
c    Set Variables
c	   ntemp = 231
c          ncstar = 109


c   Search for T (seach through all T's)
c	Identify index that is upper-bound on values
c        eg. if tempk=298.3, return index for temp=298.5
      if (tempk.lt.200) then
	ctemp = 1
      elseif (tempk.gt.315) then
        ctemp = ntemp
      else
	do tempcnt = 1,ntemp
	    if (dhtemp(tempcnt).ge.tempk) then
		ctemp = tempcnt
		goto 300
	    endif
	enddo
      endif
 300  continue

c    Iterate on deltah and Cstar to get Cstemp
c	Search for Cstar index
	    do i = 1,ntot
	    	csat1 = csat(i)
		csat2 = 0
		iter = 0
		do while (abs(log10(csat1)-log10(csat2)).ge.0.1.and.iter.le.20)
		    iter = iter + 1
		    csat2 = csat1
		    do icstar = 1,ncstar
			if (dhcstar(icstar).ge.log10(csat1)) then
			    if (i.le.20) then
				deltah(i) = poadhvap(ctemp,icstar)*1000
			    else
				deltah(i) = soadhvap(ctemp,icstar)*1000
c				if (deltah(i).lt.8) deltah(i) = 8
			    endif
			    csat1=csat(i)*(cstemp(i)/tempk)*exp((deltah(i)/8.314)
     &						*(1/cstemp(i)-1/tempk))
			    goto 301
			endif
		    enddo
 301		    continue
		enddo

		csatT(i) = csat1
	    enddo

ccccc ----------- End New Hvap Look-Up Table --------------


c
c CALCULATE AEROSOL-PHASE CONCENTRATION (CAER) AND GAS-PHASE 
c CONCENTRATION (CGAS) FOR NON-SOLUTION-FORMING COMPOUNDS
c COMPOUNDS THAT HAVE A CONCENTRATION OF LESS THAN conmin ARE IGN0RED
c MAP COMPOUNDS THAT FORM SOLUTIONS ONTO ARRAYS
c
      icont=0
      do i=1,ntot
         if (flagsoap(i).eq.0) then
            cgas(i) = amin1(ctot(i), csatT(i))
            caer(i) = ctot(i) - cgas(i)
         elseif (ctot(i).lt.conmin) then
            cgas(i) = ctot(i)
            caer(i) = 0.0
         else
            icont=icont+1
            idx(icont) = i
            smw(icont)=mwsoap(i)
            scsat(icont)=csatT(i)
            sctot(icont)=ctot(i)
            scaer(icont)=caer(i)
         endif
      enddo
      nsol=icont
c
c Check for a trivial solution
c
      if (nsol.eq.0) goto 1000
      if (nsol.eq.1) then
         if (cpre.lt.cpremin) then
           scgas(1) = amin1(sctot(1), scsat(1))
           scaer(1) = sctot(1) - scgas(1)
         else ! This case has an analytical solution
           bb = scsat(1)-sctot(1)+cpx*smw(1)
           cc = -sctot(1)*cpx*smw(1)
           scaer(1) = amin1( sctot(1), .5*(-bb+SQRT(bb*bb-4.*cc)) )
           scgas(1) = sctot(1) - scaer(1)
         endif
         goto 900
      endif
      sum=0.0
      do i=1,nsol
         sum = sum + sctot(i)/scsat(i)
      enddo
      if (cpre.lt.cpremin .and. sum.le.1.0) then
         do i=1,nsol
            scgas(i)=sctot(i)
            scaer(i)=0.0
         enddo
         goto 900
      endif
c
c Find the solution using a bi-section method (approach from max)
c
      xend = 0.0
      do i = 1, nsol
        xend = xend + oaro(i)*sctot(i)/smw(i)  !BNM density correction
      enddo
      xend = xend + cpx
      call spfcn (nsol,sctot,scsat,scaer,smw,cpx,xend,fend,oaro)
      if (abs(fend).le.xtol*xend) goto 99
      if (fend.gt.0.0) then
        write (iout,'(//,a)') ' ERROR in SOAP:'
        write (iout,'(/,a)') ' ERROR: positive end point. '
        write (iout,'(/,a,e15.6)') ' fend = ',fend
        write (iout,'(/,a,e15.6)') ' xtol*xend = ',xtol*xend
        goto 50
      endif
      dx = xend - cpx
      do znum = 1, 200
        dx = 0.5 * dx
        xmid = xend - dx
        call spfcn (nsol,sctot,scsat,scaer,smw,cpx,xmid,fmid,oaro)
        if (abs(fmid).le.xtol*xmid .or. dx.le.xtol*xmid) goto 100
        if (fmid.lt.0.0) xend = xmid
      enddo
      write (iout,'(//,a)') ' ERROR in SOAP:'
      write (iout,'(/,a)') ' ERROR: max number of iterations reached'
 50   write (iout,'(a,i3,a,3i4)') ' Grid = ', igrdchm,
     &                 ' cell(i,j,k) = ', ichm, jchm, kchm
      write (iout,'(a5,2a15)') ' spec','total [ug/m3]','c* [ug/m3]'
      write (iout,'(i5,1p2e15.6)') (idx(i),sctot(i),scsat(i),i=1,nsol)
      write (iout,'(a5,e15.6)') ' cpre',cpre
      call camxerr()
c
c Converged
c
  99  xmid = xend
 100  continue
      do i=1,nsol
         scaer(i) = amin1( sctot(i), scaer(i) )
         scgas(i) = sctot(i) - scaer(i)
      enddo
c
c REMAP COMPOUNDS THAT FORM SOLUTIONS BACK ONTO ORIGINAL ARRAYS
c      
 900  continue
      do i=1,nsol
         caer(idx(i))=scaer(i)
         cgas(idx(i))=scgas(i)
      enddo
c
c Convert to ppm if inputs in ppm
c
 1000 continue
      if (lppm) then
         do i=1,ntot
            cgas(i) = cgas(i)/(convfac*mwsoap(i))
         enddo
      endif
c
      return
      end
