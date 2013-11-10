      subroutine readchm
c
c-----CAMx v4.02 030709
c
c     READCHM reads the CAMx chemistry parameter file, which defines
c     the chemical system to be simulated
c
c     Copyright 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2003
c     ENVIRON International Corporation
c  
c     Modifications:
c        10/26/99  Added check for zero Henry's law constants to avoid a
c                  divide by 0 later in the code
c        10/20/00  Added CAMx version as first record on chemistry parameters
c                  file
c        1/9/02    Aerosol size cut points and density now defined on
c                  chemistry parameters file; removed conversion of aerosol
c                  BDNL values from ug/m3 to umol/m3
c        1/18/02   Added EC, PFIN, and PCRS to mechanism 4 species list
c        12/12/02  Expanded species list for Mechanism 4
c        3/26/03   Added surface resistance scaling factor to gas params
c
c	 12/23/08  Added Enthalpy of Vaporization table look-up for SOAP module. (BNM)
c
c
c     Input arguments: 
c        none 
c 
c     Output arguments: 
c            
c     Routines called: 
c        EXPTBL
c        KTHERM
c            
c     Called by:
c        CHMPREP
c
      include 'camx.prm'
      include 'chmstry.com'
      include 'filunit.com'
      include 'flags.com'
      include 'ddmchm.com'
      include 'iehchem.com'
      include 'soap.com'
c
      parameter(ncrsspc = 2)
      character*180 record
      character*10 nametmp, splist(NSPNAM), radlist(NRADNM) 
      character*10 nmrad1(14), nmrad2(12), nmrad5(18)
      character*10 crsspc(ncrsspc),spcsz
      character*10 blank,camxv,camxvin
      character*10 tmpnam, tmpnam1
      integer mchgas(10), mchaero(10), mchrad(10), mchrxn(10)
      integer mchphot(10), mchfast(10), mchiessr(10)
      integer ipigcb4(9), ipigsap(9)
      integer rxntyp(MXRXN), rxnord(MXRXN), npar(7)
      real rxnpar(MXRXN,12)
      real kdum(MXRXN,3), tdum(3), pdum(3)
      real*8  dsec_i(MXSECT+1)
      integer omp_get_num_procs

c BNM - variables for Hvap look-up
      integer icstar, tempcnt
c

c 
c-----Data that define the mechanism/solver options
c     The fast state species must come first in SPLIST
c
      data camxv  /'VERSION4.0'/
      data blank  /'BLANK     '/
c
      data splist /'NO        ','NO2       ','O3        ',
     &             'PAN       ','CRES      ','PAN2      ',
     &             'MPAN      ','PBZN      ','NPHE      ',
     &             'RNO3      ','DCB2      ','DCB3      ',
     &             'HNO4      ','ACET      ','ALD2      ',
     &             'ALK1      ','ALK2      ','ALK3      ',
     &             'ALK4      ','ALK5      ','ARO1      ',
     &             'ARO2      ','BACL      ','BALD      ',
     &             'BCL1      ','BCL2      ','BUTA      ',
     &             'CCHO      ','CCRS      ','CPO1      ',
     &             'CPO2      ','CPO3      ','CPO4      ',
     &             'CPO5      ','CPO6      ',
     &             'CPO7      ','CPO8      ',
     &             'CPO9      ','CPO10     ','COO1      ',
     &             'COO2      ','COO3      ','COO4      ',
     &             'COO5      ','COO6      ','COO7      ',
     &             'COO8      ','COO9      ',
     &             'COO10     ','CBS1      ',
     &             'CBS2      ','CBS3      ','CBS4      ',
     &             'CAS1      ','CAS2      ',
     &             'CAS3      ','CAS4      ',
c
     &             'CL2       ','CO        ','CO2H      ',
     &             'CO3H      ','COOH      ','CPRM      ',
     &             'DCB1      ','ETH       ','ETHE      ',
     &             'ETOH      ','FCRS      ','FMCL      ',
     &             'FORM      ','FPRM      ','GLY       ',
     &             'H2O2      ','HC2H      ','HCHO      ',
     &             'HCL       ','HONO      ','HNO3      ',
     &             'HO2H      ','HOCL      ','ICL1      ',
     &             'ICL2      ','ISOP      ','ISPD      ',
     &             'MEK       ','MEOH      ','METH      ',
     &             'MGLY      ','MVK       ','NA        ',
     &             'NH3       ','NTR       ','NXOY      ',
     &             'OLE       ','OLE1      ','OLE2      ',
     &             'BPIN      ','LIMO      ','MONO      ',
     &             'SESQ      ',
     &             'OPEN      ','PAR       ','PCL       ',
     &             'PEC       ','PHEN      ','PNA       ',
     &             'PNH4      ','PNO3      ','POA       ',
     &             'PROD      ','PSO4      ','RC2H      ',
     &             'RC3H      ','RCHO      ','ROOH      ',
     &             'SO2       ',
     &             'APO1      ',
     &             'APO2      ','APO3      ','APO4      ',
     &             'APO5      ','APO6      ',
     &             'APO7      ','APO8      ',
     &             'APO9      ','APO10     ','AOO1      ',
     &             'AOO2      ','AOO3      ','AOO4      ',
     &             'AOO5      ','AOO6      ','AOO7      ',
     &             'AOO8      ','AOO9      ',
     &             'AOO10     ','ABS1      ',
     &             'ABS2      ','ABS3      ','ABS4      ',
     &             'AAS1      ','AAS2      ',
     &             'AAS3      ','AAS4      ',
     &             'SULF      ',
     &             'TERP      ','TOL       ','XN        ',
     &             'XYL       ',
     &             'APO1_1    ',
     &             'APO1_2    ','APO1_3    ','APO1_4    ',
     &             'APO1_5    ','APO1_6    ','APO1_7    ',
     &             'APO1_8    ','APO1_9    ','APO1_10   ',
     &             'APO2_1    ',
     &             'APO2_2    ','APO2_3    ','APO2_4    ',
     &             'APO2_5    ','APO2_6    ','APO2_7    ',
     &             'APO2_8    ','APO2_9    ','APO2_10   ',
     &             'APO3_1    ',
     &             'APO3_2    ','APO3_3    ','APO3_4    ',
     &             'APO3_5    ','APO3_6    ','APO3_7    ',
     &             'APO3_8    ','APO3_9    ','APO3_10   ',
     &             'APO4_1    ',
     &             'APO4_2    ','APO4_3    ','APO4_4    ',
     &             'APO4_5    ','APO4_6    ','APO4_7    ',
     &             'APO4_8    ','APO4_9    ','APO4_10   ',
     &             'APO5_1    ',
     &             'APO5_2    ','APO5_3    ','APO5_4    ',
     &             'APO5_5    ','APO5_6    ','APO5_7    ',
     &             'APO5_8    ','APO5_9    ','APO5_10   ',
     &             'APO6_1    ',
     &             'APO6_2    ','APO6_3    ','APO6_4    ',
     &             'APO6_5    ','APO6_6    ','APO6_7    ',
     &             'APO6_8    ','APO6_9    ','APO6_10   ',
     &             'APO7_1    ',
     &             'APO7_2    ','APO7_3    ','APO7_4    ',
     &             'APO7_5    ','APO7_6    ','APO7_7    ',
     &             'APO7_8    ','APO7_9    ','APO7_10   ',
     &             'APO8_1    ',
     &             'APO8_2    ','APO8_3    ','APO8_4    ',
     &             'APO8_5    ','APO8_6    ','APO8_7    ',
     &             'APO8_8    ','APO8_9    ','APO8_10   ',
     &             'APO9_1    ',
     &             'APO9_2    ','APO9_3    ','APO9_4    ',
     &             'APO9_5    ','APO9_6    ','APO9_7    ',
     &             'APO9_8    ','APO9_9    ','APO9_10   ',
     &             'APO10_1   ',
     &             'APO10_2   ','APO10_3   ','APO10_4   ',
     &             'APO10_5   ','APO10_6   ','APO10_7   ',
     &             'APO10_8   ','APO10_9   ','APO10_10  ',
     &             'AOO1_1    ',
     &             'AOO1_2    ','AOO1_3    ','AOO1_4    ',
     &             'AOO1_5    ','AOO1_6    ','AOO1_7    ',
     &             'AOO1_8    ','AOO1_9    ','AOO1_10   ',
     &             'AOO2_1    ',
     &             'AOO2_2    ','AOO2_3    ','AOO2_4    ',
     &             'AOO2_5    ','AOO2_6    ','AOO2_7    ',
     &             'AOO2_8    ','AOO2_9    ','AOO2_10   ',
     &             'AOO3_1    ',
     &             'AOO3_2    ','AOO3_3    ','AOO3_4    ',
     &             'AOO3_5    ','AOO3_6    ','AOO3_7    ',
     &             'AOO3_8    ','AOO3_9    ','AOO3_10   ',
     &             'AOO4_1    ',
     &             'AOO4_2    ','AOO4_3    ','AOO4_4    ',
     &             'AOO4_5    ','AOO4_6    ','AOO4_7    ',
     &             'AOO4_8    ','AOO4_9    ','AOO4_10   ',
     &             'AOO5_1    ',
     &             'AOO5_2    ','AOO5_3    ','AOO5_4    ',
     &             'AOO5_5    ','AOO5_6    ','AOO5_7    ',
     &             'AOO5_8    ','AOO5_9    ','AOO5_10   ',
     &             'AOO6_1    ',
     &             'AOO6_2    ','AOO6_3    ','AOO6_4    ',
     &             'AOO6_5    ','AOO6_6    ','AOO6_7    ',
     &             'AOO6_8    ','AOO6_9    ','AOO6_10   ',
     &             'AOO7_1    ',
     &             'AOO7_2    ','AOO7_3    ','AOO7_4    ',
     &             'AOO7_5    ','AOO7_6    ','AOO7_7    ',
     &             'AOO7_8    ','AOO7_9    ','AOO7_10   ',
     &             'AOO8_1    ',
     &             'AOO8_2    ','AOO8_3    ','AOO8_4    ',
     &             'AOO8_5    ','AOO8_6    ','AOO8_7    ',
     &             'AOO8_8    ','AOO8_9    ','AOO8_10   ',
     &             'AOO9_1    ',
     &             'AOO9_2    ','AOO9_3    ','AOO9_4    ',
     &             'AOO9_5    ','AOO9_6    ','AOO9_7    ',
     &             'AOO9_8    ','AOO9_9    ','AOO9_10   ',
     &             'AOO10_1   ',
     &             'AOO10_2   ','AOO10_3   ','AOO10_4   ',
     &             'AOO10_5   ','AOO10_6   ','AOO10_7   ',
     &             'AOO10_8   ','AOO10_9   ','AOO10_10  ',
     &             'ABS1_1    ',
     &             'ABS1_2    ','ABS1_3    ','ABS1_4    ',
     &             'ABS1_5    ','ABS1_6    ','ABS1_7    ',
     &             'ABS1_8    ','ABS1_9    ','ABS1_10   ',
     &             'ABS2_1    ',
     &             'ABS2_2    ','ABS2_3    ','ABS2_4    ',
     &             'ABS2_5    ','ABS2_6    ','ABS2_7    ',
     &             'ABS2_8    ','ABS2_9    ','ABS2_10   ',
     &             'ABS3_1    ',
     &             'ABS3_2    ','ABS3_3    ','ABS3_4    ',
     &             'ABS3_5    ','ABS3_6    ','ABS3_7    ',
     &             'ABS3_8    ','ABS3_9    ','ABS3_10   ',
     &             'ABS4_1    ',
     &             'ABS4_2    ','ABS4_3    ','ABS4_4    ',
     &             'ABS4_5    ','ABS4_6    ','ABS4_7    ',
     &             'ABS4_8    ','ABS4_9    ','ABS4_10   ',
     &             'AAS1_1    ',
     &             'AAS1_2    ','AAS1_3    ','AAS1_4    ',
     &             'AAS1_5    ','AAS1_6    ','AAS1_7    ',
     &             'AAS1_8    ','AAS1_9    ','AAS1_10   ',
     &             'AAS2_1    ',
     &             'AAS2_2    ','AAS2_3    ','AAS2_4    ',
     &             'AAS2_5    ','AAS2_6    ','AAS2_7    ',
     &             'AAS2_8    ','AAS2_9    ','AAS2_10   ',
     &             'AAS3_1    ',
     &             'AAS3_2    ','AAS3_3    ','AAS3_4    ',
     &             'AAS3_5    ','AAS3_6    ','AAS3_7    ',
     &             'AAS3_8    ','AAS3_9    ','AAS3_10   ',
     &             'AAS4_1    ',
     &             'AAS4_2    ','AAS4_3    ','AAS4_4    ',
     &             'AAS4_5    ','AAS4_6    ','AAS4_7    ',
     &             'AAS4_8    ','AAS4_9    ','AAS4_10   ',
     &             'POC_1     ',
     &             'POC_2     ','POC_3     ','POC_4     ',
     &             'POC_5     ','POC_6     ','POC_7     ',
     &             'POC_8     ','POC_9     ','POC_10    ',
     &             'PEC_1     ','PEC_2     ','PEC_3     ',
     &             'PEC_4     ','PEC_5     ','PEC_6     ',
     &             'PEC_7     ','PEC_8     ','PEC_9     ',
     &             'PEC_10    ','CRST_1    ','CRST_2    ',
     &             'CRST_3    ','CRST_4    ','CRST_5    ',
     &             'CRST_6    ','CRST_7    ','CRST_8    ',
     &             'CRST_9    ','CRST_10   ','PH2O_1    ',
     &             'PH2O_2    ','PH2O_3    ','PH2O_4    ',
     &             'PH2O_5    ','PH2O_6    ','PH2O_7    ',
     &             'PH2O_8    ','PH2O_9    ','PH2O_10   ',
     &             'PCL_1     ','PCL_2     ','PCL_3     ',
     &             'PCL_4     ','PCL_5     ','PCL_6     ',
     &             'PCL_7     ','PCL_8     ','PCL_9     ',
     &             'PCL_10    ','NA_1      ','NA_2      ',
     &             'NA_3      ','NA_4      ','NA_5      ',
     &             'NA_6      ','NA_7      ','NA_8      ',
     &             'NA_9      ','NA_10     ','PNH4_1    ',
     &             'PNH4_2    ','PNH4_3    ','PNH4_4    ',
     &             'PNH4_5    ','PNH4_6    ','PNH4_7    ',
     &             'PNH4_8    ','PNH4_9    ','PNH4_10   ',
     &             'PNO3_1    ','PNO3_2    ','PNO3_3    ',
     &             'PNO3_4    ','PNO3_5    ','PNO3_6    ',
     &             'PNO3_7    ','PNO3_8    ','PNO3_9    ',
     &             'PNO3_10   ','PSO4_1    ','PSO4_2    ',
     &             'PSO4_3    ','PSO4_4    ','PSO4_5    ',
     &             'PSO4_6    ','PSO4_7    ','PSO4_8    ',
     &             'PSO4_9    ','PSO4_10   ','PH2O      '/
c
      data crsspc /'CCRS      ','CPRM      '/
c
      data radlist/'O1D       ','O         ','CLO       ',
     &             'CL        ','N2O5      ','NO3       ',
     &             'OH        ','HO2       ','C2O3      ',
     &             'XO2       ','XO2N      ','TO2       ',
     &             'ROR       ','CRO       ','RO2R      ',
     &             'R2O2      ','RO2N      ','CCO3      ',
     &             'RCO3      ','MCO3      ','BZCO      ',
     &             'CXO2      ','HCO3      ','TBUO      ',
     &             'BZO       ','BZNO      '/
c
      data nmrad1 /'O1D       ','CL        ','CLO       ','O         ',
     &             'N2O5      ','NO3       ','OH        ','HO2       ',
     &             'C2O3      ','XO2       ','XO2N      ','TO2       ',
     &             'ROR       ','CRO       '/
      data nmrad2 /'O1D       ','O         ',
     &             'N2O5      ','NO3       ','OH        ','HO2       ',
     &             'C2O3      ','XO2       ','XO2N      ','TO2       ',
     &             'ROR       ','CRO       '/
      data nmrad5 /'O1D       ','O         ','N2O5      ',
     &             'NO3       ','OH        ','HO2       ',
     &             'RO2R      ','R2O2      ','RO2N      ',
     &             'CCO3      ','RCO3      ','MCO3      ','BZCO      ',
     &             'CXO2      ','HCO3      ','TBUO      ',
     &             'BZO       ','BZNO      '/
c
      data mchgas   / 34, 24, 25, 34, 91, 78,  0, 0, 0, 0 /
      data mchaero  /  0,  0,  0, 16, 37, 53,  0, 0, 0, 0 /
      data mchrad   / 14, 12, 12, 12, 18, 12,  0, 0, 0, 0 /
      data mchiessr /  4,  2,  2,  2,  2,  2,  0, 0, 0, 0 /
      data mchrxn   /110, 91, 96,100,237,100,  0, 0, 0, 0 /
      data mchphot  / 14, 11, 12, 12, 30, 12,  0, 0, 0, 0 /
      data mchfast  /  4,  4,  4,  4, 13,  4,  0, 0, 0, 0 /
      data mchidmin / 1 /
      data mchidmax / 6 /
      data ipigcb4  /  7,  9, 10, 11, 16, 17, 18, 19, 20 /
      data ipigsap  /  8, 18, 20, 19, 14, 11, 13, 12, 10 /
      data npar     /  1,  2,  4, 10,  5, 12,  8 /
      data tdum     / 298.,  273.,  298. /
      data pdum     / 1013., 1013., 491. /
c
c     Many rules about names, number and order of species are 
c     enforced here unless LCHEM is false.
c
c     Arrays MCHGAS, MCHAERO, MCHRXN and MCHPHOT allow for up to 10
c     mechanisms to be called by RADDRIVR and CHEMDRIV. 
c     Currently:
c      1 = CB4 (mech 3) + Chlorine
c      2 = CB4 + Updated radical-radical reactions
c      3 = CB4 + Carter one product isoprene mechanism
c      4 = CB4 + "1-atmosphere" aerosol mechanism
c      5 = SAPRC99 (56 specs, 221 rxns, 30 phot rxns)
c
c     The state species solver (TRAP) expects the fast species to be
c     ordered first in SPLIST (NO, NO2, O3, PAN etc) so don't change these.
c     Only state species that are named in SPLIST will be allowed.
c     Update CHMSTRY.COM if new species are added to SPLIST
c
c     The radical solvers (RADINIT, RADSLVR) expect radical species
c     in specific orders as defined in NMRAD lists, RADLIST is the
c     union of these lists in order matching CHMSTRY.COM equivalence
c
c
c-----Entry point
c 
c
c-----Mechanism ID
c
      idmech = 0
      read(ichem,'(a)') record
      camxvin = record(21:30)
      call jstlft(camxvin)
      call toupper(camxvin)
      if (camxvin.ne.camxv) then
        write(iout,'(/,a)') ' CAMx version in CHEMPARAM file INVALID'
        write(iout,'(a,a)') ' Expecting: ',camxv
        write(iout,'(a,a)') '     Found: ',camxvin
        goto 910
      endif
      read(ichem,'(a)') record
      read(record(21:80),*) idmech
      write(idiag,'(a,i4)') 'Using CHEMPARAM mechanism id ',idmech
      if (lchem) then
        if (idmech.lt.mchidmin .or. idmech.gt.mchidmax) then
          write(iout,'(/,a)') ' CHEMPARAM mechanism id INVALID'
          write(iout,'(a,i3,a,i3)')' use an id between',mchidmin,
     &                             'and',mchidmax
          goto 910
        endif
c
c --- Mechanism 4 is not allowed with OMP,
c     For compilers with no OMP support you should link in the 
c     routine in the dummy.f file.
c
        ncpus = omp_get_num_procs()
        if( (idmech .EQ. 4 .OR. idmech .EQ. 6) .AND. ncpus .GT. 1 ) then
           write(iout,'(//,a)') 'ERROR in READCHM:'
           write(iout,'(1X,2A)') 'OMP (multiprocessing) is not ',
     &                           'supported for Mechanism 4 or 6.'
           write(iout,'(1X,2A)') 'You must either compile without ',
     &                            'OMP or choose another mechanism.'
           call camxerr()
        endif
      endif
      read(ichem,'(a)') record
      write(idiag,'(a)') record(:istrln(record))
c
c-----Set radical species order. Most mechanisms use same
c     order as 2.  Save final name order.
c
      if (lchem) then
        nrad = mchrad(idmech)
        iessrad = mchiessr(idmech)
        if (nrad.gt.MXRADCL) then
          write(iout,'(/,a,i5,a,i5,a)') ' The number of radicals ',
     &         nrad, ' exceeds the value of MXRADCL ', MXRADCL,
     &         ' set in camx.prm'
          goto 910
        endif
c
        do l = 1,NRADNM
          krad(l) = nrad + 1
        enddo
c
        if (idmech.eq.1) then
          do i = 1,nrad
            do j = 1,nradnm
              if (nmrad1(i).eq.radlist(j)) then
                 krad(j) = i
                 nmrad(i) = radlist(j)
              endif
            enddo
          enddo
        endif
c
        if (idmech.eq.2.or.idmech.eq.3.or.idmech.eq.4
     &                 .or.idmech.eq.6) then
          do i = 1,nrad
            do j = 1,nradnm
              if (nmrad2(i).eq.radlist(j)) then
                 krad(j) = i
                 nmrad(i) = radlist(j)
              endif
            enddo
          enddo
        endif
c
        if (idmech.eq.5) then
          do i = 1,nrad
            do j = 1,nradnm
              if (nmrad5(i).eq.radlist(j)) then
                 krad(j) = i
                 nmrad(i) = radlist(j)
              endif
            enddo
          enddo
        endif
c
        nmrad(nrad+1) = blank
        write(idiag,'(a)') 'The radicals are'
        write(idiag,'(6a10)') (nmrad(i),i=1,nrad)
      endif
c
c-----Species number and reaction number
c
      read(ichem,'(a)') record
      read(record(21:80),*) ngas
      read(ichem,'(a)') record
      read(record(21:80),*) naero
      nspec = ngas + naero 
      if (idmech.EQ.6 .AND. naero.GT.0) then
        read(record(21:80),*) naero,nsec_c,dt_aero
        read(ichem,'(a)') record
        read(record(21:),*) (dsec_i(i),i=1,nsec_c+1)
        nspec = ngas + naero * nsec_c
      endif
      if (idmech.EQ.5 .AND. naero.GT.0) then
        read(record(21:80),*) naero,nsec_c,dt_aero
        read(ichem,'(a)') record
        read(record(21:),*) (dsec_i(i),i=1,nsec_c+1)
        nspec = ngas + naero * nsec_c
      endif
      read(ichem,'(a)') record
      read(record(21:80),*) nreact
      if (nspec.lt.1) then 
        write(iout,'(/,a,i5,a)') 
     &    ' number of GAS + AERO species on CHEMPARAM file =', nspec, 
     &    ' is less than 1' 
          goto 910 
      endif
c 
c-----Check number of species and number of reactions against chosen
c     mechanism
c 
      if (lchem) then
        if (ngas.gt.mchgas(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' number of GAS species on CHEMPARAM file =', ngas,
     &     ' greater than ngas for this mechanism =', mchgas(idmech)
          goto 910
        endif
c
        if (naero.gt.mchaero(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' number of AERO species on CHEMPARAM file =', naero,
     &     ' greater than naero for this mechanism =', mchaero(idmech)
          goto 910
        endif
c
        if (nreact.ne.mchrxn(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' number of reactions on CHEMPARAM file =',nreact,
     &     ' not equal to reactions for this mechanism =',mchrxn(idmech)
          goto 910
        endif
      endif
c
c-----Check dimensional limit of species and reaction numbers
c
      if (nspec.gt.MXSPEC) then
        write(iout,'(/,a,i5,a,i5,a)') ' The number of species ',
     &       nspec, ' exceeds the value of MXSPEC ', MXSPEC,
     &       ' set in camx.prm'
        goto 910
      endif
c
      if (nreact.gt.MXRXN) then
        write(iout,'(/,a,i5,a,i5,a)') ' The number of reactions ',
     &       nreact, ' exceeds the value of MXRXN ', MXRXN,
     &       ' set in camx.prm'
        goto 910
      endif
c
c-----Read primary photolysis ID record
c
      read(ichem,'(a)') record
      read(record(21:80),*) nphot1,(idphot1(n),n=1,nphot1) 
      if (nphot1.gt.0) then
        write(idiag,'(a)') 'The primary photolysis reactions are'
        write(idiag,'(i6)') (idphot1(n),n=1,nphot1) 
      endif
c
c-----Read secondary (scaled) photolysis ID records
c
      read(ichem,'(a)') record
      read(record(21:80),*) nphot2
      if (nphot2.gt.0) then
        if (nphot1.eq.0) then
          write(iout,'(/,a)') 
     &     'Need at least one primary photolysis reaction'
          goto 910
        endif
        do n = 1,nphot2
          read(ichem,'(a)') record
          read(record(21:80),*) idphot2(n),idphot3(n),phtscl(n)
        enddo
        write(idiag,'(a)') 'The secondary photolysis reactions are' 
        write(idiag,'(i6,a,i6,a,1pe10.3)')  
     &    (idphot2(n),' =',idphot3(n),' *',phtscl(n),n=1,nphot2)
      endif
      nphot = nphot1 + nphot2
c
c-----Check photolysis ID records
c
      if (lchem) then
c
        if (nphot1.lt.1) then 
          write(iout,'(/,a)') 
     &      ' Need at least one primary photolysis reaction' 
          goto 910 
        endif 
c
        if (nphot1.gt.MXPHT1) then
          write(iout,'(/,a,i5)') 
     &     ' Number of primary photolysis reactions exceeds max of ',
     &       MXPHT1
          goto 910
        endif
c
        if (nphot2.gt.MXPHT2) then
          write(iout,'(/,a,i5)') 
     &     ' Number of secondary photolysis reactions exceeds max of ',
     &       MXPHT2
          goto 910
        endif
c
        if (nphot.ne.mchphot(idmech)) then
          write(iout,'(/,a,i5,/,a,i5)')
     &     ' Chemistry mechanism requires', mchphot(idmech),
     &     ' photolysis reactions, but CHEMPARAM file has' , nphot
          goto 910
        endif
c
        nerr = 0
        do n = 1,nphot1
          if (idphot1(n).gt.nreact) nerr = 1
        enddo
c
        do n = 1,nphot2
          if (idphot2(n).gt.nreact) nerr = 1
          if (idphot3(n).gt.nreact) nerr = 1
          nhit = 0
          do nn = 1,nphot1
            if (idphot3(n).eq.idphot1(nn)) nhit = 1
          enddo
          if (nhit.eq.0) nerr = 1
        enddo
c
        if (nerr.eq.1) then
           write(iout,'(/,a)') 'ERROR in the CHEMPARAM file:'
           write(iout,'(a)') 
     &     ' Bad reaction number in one of the photolysis reaction IDs.'
          goto 910
        endif
      endif
c
c-----Species records, gases come first
c
      read(ichem,'(a)') record
      read(ichem,'(a)') record
      write(idiag,'(a)') 'The state species are'
      write(idiag,'(a)') record(:istrln(record))
      do l=1,ngas
        read(ichem,'(5x,a10,2e10.0,4f10.0)')
     &             spname(l),bdnl(l),henry0(l),tfact(l),
     &             diffrat(l),f0(l),rscale(l)
        rscale(l) = amin1(1.,rscale(l))
        rscale(l) = amax1(0.,rscale(l))
        write(idiag,'(i3,2x,a10,2e10.2,4f10.2)')
     &             l,spname(l),bdnl(l),henry0(l),tfact(l),
     &             diffrat(l),f0(l),rscale(l)
      enddo
c
c-----Check over the gas phase species
c
      if (lchem) then
        do l=1,ngas
          if (henry0(l).eq.0.) then
            write(iout,'(/,a,i5)')
     &      'The Henry0 value must be non-zero for species ',l
            goto 910
          endif
        enddo
c
c-----Reorder fast species (NO, NO2, O3, and PAN etc) to come first
c
        nspfst=mchfast(idmech)
        do i=1,nspfst
          do j=1,ngas
            if (splist(i).eq.spname(j)) then
              nametmp = spname(i)
              spname(i) = spname(j)
              spname(j) = nametmp
              tmp = bdnl(i)
              bdnl(i) = bdnl(j)
              bdnl(j) = tmp
              tmp = henry0(i)
              henry0(i) = henry0(j)
              henry0(j) = tmp
              tmp = tfact(i)
              tfact(i) = tfact(j)
              tfact(j) = tmp
              tmp = diffrat(i)
              diffrat(i) = diffrat(j)
              diffrat(j) = tmp
              tmp = f0(i)
              f0(i) = f0(j)
              f0(j) = tmp
              tmp = rscale(i)
              rscale(i) = rscale(j)
              rscale(j) = tmp
            endif
          enddo
        enddo
      endif
c
c-----Read the aero species, if any
c
      if (naero.ne.0) then
        read(ichem,'(a)') record
        write(idiag,'(a)') record(:istrln(record))
        if (idmech.EQ.6.or.idmech.EQ.5) then
          do iaero = 1, naero
            read(ichem,'(5x,a10,e10.0,f10.0)')
     &           tmpnam,bdnl_tmp,roprt_tmp
            nl=istrln(tmpnam)
            do isec = 1, nsec_c
              l = ngas + (iaero-1) * nsec_c + isec
              if ( isec .ge. 10 ) then
                write(tmpnam1,'(a1,i2)') '_',isec
                spname(l)=tmpnam(1:nl)//tmpnam1(1:3)
              else
                write(tmpnam1,'(a1,i1)') '_',isec
                spname(l)=tmpnam(1:nl)//tmpnam1(1:2)
              endif
              bdnl(l) = bdnl_tmp
              roprt(l) = roprt_tmp * 1.e6
              dcut(l,1) = dsec_i(isec)
              dcut(l,2) = dsec_i(isec+1)
              write(idiag,'(i3,2x,a10,e10.2,3f10.2)')
     &           l,spname(l),bdnl(l),roprt_tmp,(dcut(l,m),m=1,2)
            enddo
          enddo
        else
          do l = ngas+1,nspec
            read(ichem,'(5x,a10,e10.0,f10.0)')
     &         spname(l),bdnl(l),roprt(l)
            dcut(l,1) = 0.04
            dcut(l,2) = 2.50
            spcsz = 'FINE      '
            do j = 1, ncrsspc
              if (spname(l).eq.crsspc(j)) then
                dcut(l,1) = 4.30
                dcut(l,2) = 10.0
                spcsz = 'COARSE    '
              endif
            enddo
            write(idiag,'(i3,2x,a10,e10.2,f10.2,2x,a10)')
     &         l,spname(l),bdnl(l),roprt(l),spcsz
            roprt(l) = roprt(l)*1.e6
          enddo
        endif
      endif
c
c-----Map species names to the internal name list.  The default setting
c     to nspec+1 allows species that are in the chem solvers to be
c     omitted from the species list for this run.  Testing if a named
c     pointer is set to nspec+1 is used to identify species that are not
c     in the run.  Named pointers (e.g., kno) are set by equivalence
c     in CHMSTRY.COM
c
      if (lchem) then
        do j = 1,NSPNAM
          kmap(j) = nspec+1
        enddo
c      
        do 10 j = 1,nspec
          do i = 1,NSPNAM
            if (splist(i).eq.spname(j)) then
              kmap(i) = j
              goto 10
            endif
          enddo
          write(iout,'(/,3a)') 'species ', spname(j),
     &                        ' is not in the internal list'
          goto 910
  10    continue
        spname(nspec+1) = blank
        write(idiag,'(a)')'Internal species order is:'
        write(idiag,'(6a10)') (spname(i),i=1,nspec) 
c
c-----Check for a consistent set of ammonium/nitrate/sulfate 
c     species
c
        if (idmech.EQ.6.) then ! MECH 6
          if (kso2.eq.nspec+1    .or. kh2o2.eq.nspec+1
     &   .or. kform.eq.nspec+1   .or. khono.eq.nspec+1
     &   .or. ko3.eq.nspec+1     .or. koh.eq.nspec+1
     &   .or. kho2.eq.nspec+1    .or. kno3.eq.nspec+1
     &   .or. kno.eq.nspec+1     .or. kno2.eq.nspec+1
     &   .or. kpan.eq.nspec+1 
     &   .or. khno3.eq.nspec+1
     &   .or. knh3.eq.nspec+1    .or. ksulf.eq.nspec+1
     &   .or. khcl.eq.nspec+1  
     &   .or. kcrust_1.eq.nspec+1
     &   .or. kpoc_1.eq.nspec+1  .or. kpec_1.eq.nspec+1
     &   .or. kph2o_1.eq.nspec+1 .or. kpcl_1.eq.nspec+1
     &   .or. kna_1.eq.nspec+1   .or. kpnh4_1.eq.nspec+1
     &   .or. kpno3_1.eq.nspec+1 .or. kpso4_1.eq.nspec+1) then
            write(iout,'(/,a)') ' You must have all of the species'
            write(iout,'(a)')   ' SO2,H2O2,FORM,HONO,O3,OH,HO2,NO3,'
            write(iout,'(a)')   ' NO,NO2,PAN,CG1,CG2,CG3,CG4,HNO3,'
            write(iout,'(a)')   ' NH3,SULF,HCL,SOA1,SOA2,SOA3,SOA4,'
            write(iout,'(a)')   ' CRST,POC,PEC,PH2O,PCL,NA,PNH4,PNO3,'
            write(iout,'(a)')   ' PSO4 to use Chemistry Mechanism 6.'
            goto 910
          endif
        elseif(idmech.EQ.5 ) then  ! MECH 5
          if (knh3.lt.nspec+1 .and. kpnh4_1.lt.nspec+1 .and.
     &       kso2.lt.nspec+1 .and. ksulf.lt.nspec+1 .and.
     &        kpso4_1.lt.nspec+1 .and. khno3.lt.nspec+1 .and.
     &         kpno3_1.lt.nspec+1 ) then
            continue
          elseif (kpnh4_1.eq.nspec+1 .and. kpso4_1.eq.nspec+1 .and.
     &         kpno3_1.eq.nspec+1 ) then
            continue
          else
            write(iout,'(/,a)') ' You must have all of the '
            write(iout,'(a)')   ' ammonium/sulfate/nitrate species '
            write(iout,'(a)')   ' NH3,NH4,SO2,SULF,PSO4,HNO3,PNO3. '
            write(iout,'(a)')   ' Or: '
            write(iout,'(a)')   ' Any combination of the gas-phase'
            write(iout,'(a)')   ' species NH3,SO2,HNO3. '
            goto 910
          endif
        else                  ! OTHER THAN MECH 6 or 5
          if (knh3.lt.nspec+1 .and. kpnh4.lt.nspec+1 .and.
     &       kso2.lt.nspec+1 .and. ksulf.lt.nspec+1 .and.
     &        kpso4.lt.nspec+1 .and. khno3.lt.nspec+1 .and.
     &         kpno3.lt.nspec+1 ) then
            continue
          elseif (kpnh4.eq.nspec+1 .and. kpso4.eq.nspec+1 .and.
     &         kpno3.eq.nspec+1 ) then
            continue
          else
            write(iout,'(/,a)') ' You must have all of the '
            write(iout,'(a)')   ' ammonium/sulfate/nitrate species '
            write(iout,'(a)')   ' NH3,NH4,SO2,SULF,PSO4,HNO3,PNO3. '
            write(iout,'(a)')   ' Or: '
            write(iout,'(a)')   ' Any combination of the gas-phase'
            write(iout,'(a)')   ' species NH3,SO2,HNO3. '
            goto 910
          endif
c
c-----Check for a consistent set of sea salt species, or none
c
          if (kna.eq.nspec+1 .and. kpcl.eq.nspec+1 .and.
     &                               khcl.eq.nspec+1) then 
            continue
          elseif (kna.lt.nspec+1 .and. kpcl.lt.nspec+1 .and.
     &                               khcl.lt.nspec+1) then 
            continue
          elseif (kna.eq.nspec+1 .and. kpcl.eq.nspec+1 .and.
     &                               khcl.lt.nspec+1) then 
            continue
          else
            write(iout,'(/,a)') ' You must have all or none of the '
            write(iout,'(a)')   ' sea salt aerosol species '
            write(iout,'(a)')   ' NA, PCL, HCL.'
            write(iout,'(a)')   ' Or, just HCL.'
            goto 910
          endif
c
c-----Check for a consistent set of CG/SOA species, or none
c
c          if (kcg1.eq.nspec+1 .and. kcg2.eq.nspec+1 .and.
c     &       kcg3.eq.nspec+1 .and. kcg4.eq.nspec+1 .and.
c     &        ksoa1.eq.nspec+1 .and. ksoa2.eq.nspec+1 .and.
c     &         ksoa3.eq.nspec+1 .and. ksoa4.eq.nspec+1 ) then
c            continue
c          elseif (kcg1.lt.nspec+1 .and. kcg2.lt.nspec+1 .and.
c     &       kcg3.lt.nspec+1 .and. kcg4.lt.nspec+1 .and.
c     &        ksoa1.lt.nspec+1 .and. ksoa2.lt.nspec+1 .and.
c     &         ksoa3.lt.nspec+1 .and. ksoa4.lt.nspec+1 ) then
c            continue
c          else
c            write(iout,'(/,a)') ' You must have all or none of the '
c            write(iout,'(a)')   ' secondary organic aerosol species '
c            write(iout,'(a)')   ' CG1,CG2,CG3,CG4,SOA1,SOA2,SOA3,SOA4 '
c            goto 910
c          endif
        endif                 ! MECH 6 ?
      endif
c
c-----Set up section diameters and check parameters for AERO routines
c
      if ( lchem .AND. idmech.EQ.6 .AND. naero.GT.0 ) then
        ierr = 0 
        call aeroset(nsec_c,dsec_i,ierr)
        if ( ierr .ne. 0 ) goto 910
      endif
      if ( lchem .AND. idmech.EQ.5 .AND. naero.GT.0 ) then
        ierr = 0
        call aeroset(nsec_c,dsec_i,ierr)
        if ( ierr .ne. 0 ) goto 910
      endif
c
c-----Reaction records
c
      if (nreact.gt.0) then
        read(ichem,'(a)') record
        read(ichem,'(a)') record
        write(idiag,'(a)') 'The reaction rate parameters are'
        do i = 1,nreact
          read(ichem,'(a)') record
          read(record,*) num,rxntyp(i)
          if (rxntyp(i).lt.1 .or. rxntyp(i).gt.7) goto 902
          read(record,*,err=900) num,rxntyp(i),rxnord(i),
     &                           (rxnpar(i,j),j=1,npar(rxntyp(i)))
c
          if (rxntyp(i).eq.5) then
            iref = anint(rxnpar(i,1))
            if (i.le.iref) goto 901
          endif
c
          write(idiag,'(3i3,1p12e12.4)')
     &       num,rxntyp(i),rxnord(i),(rxnpar(i,j),j=1,npar(rxntyp(i)))
        enddo
c
c-----Populate rate constant lookup table
c
        call exptbl(rxntyp,rxnord,rxnpar)
c
c-----Provide diagnostic info for checking rate expressions
c
        write(idiag,'(/,a,/,/,a)') 
     &        'Diagnostic info for checking rate expressions',
     &        'Rates at three temps and pressures in ppm-n min-1'
        do i=1,3
          call ktherm(tdum(i), pdum(i))
          do j=1,mxrxn
            kdum(j,i)=rk(j)/60.
          enddo
        enddo
        write(idiag,'(a,3F12.1)') 'Temp= ', tdum
        write(idiag,'(a,3F12.1)') 'Pres= ', pdum
        write(idiag,'(a)') 'Rxn Type'
        write(idiag,'(2i3,1p3e12.4)')
     &       (j,rxntyp(j),kdum(j,1),
     &        kdum(j,2),kdum(j,3),j=1,nreact)
c
c-----Set reaction pointers for pig chemistry
c
        if(idmech.le.4 .OR. idmech.eq.6)then
          do i=1,9
            ipigrxn(i)=ipigcb4(i)
          enddo
        elseif(idmech.eq.5)then
          do i=1,9
            ipigrxn(i)=ipigsap(i)
          enddo
        else
          write(iout,'(/,a)') 
     &      'Dont know how to set pig rxns for mech #', idmech
          goto 910
        endif
      endif
c
c---  Set IEH solver species pointers via equivalences in iehchem.com
c
      if (lchem) then
        do i = 1,nradnm
          if (krad(i).gt.nrad) then
            irad(i)=ngas+nrad+1
          elseif (krad(i).gt.iessrad) then
            irad(i)=krad(i)-iessrad
          else
            irad(i)=krad(i)+nspfst+nrad-iessrad
          endif
        enddo
        do i=1, nspnam
          if (kmap(i).gt.nspfst) then
            imap(i)=kmap(i)+nrad
          else
            imap(i)=kmap(i)+nrad-iessrad
          endif
        enddo
      endif
c
c======================== DDM Begin =======================
c
c---  Set species pointers via equivalences in ddmchm.com
c
      if (lchem) then
        do i = 1,nradnm
          if (krad(i).le.nrad) then
            lrad(i)=krad(i)
          else
            lrad(i)=ngas+nrad+1
          endif
        enddo
        do i=1, nspnam
          lmap(i)=kmap(i)+nrad
        enddo
      endif
c
c======================== DDM End   =======================
c
      write(idiag,*)
      call flush(idiag)
c
c============ BNM - Reading Hvap tables ===================
c
c      ccc   Using Look-up Table Compiled by Scott Epstein   ccc
c       Set Variables
c            ntemp = 231
c            ncstar = 109

c       Read in all data from look-up tables
            open (96, file='/home/bnmurphy/Research/PMCAMx/PMCAMxSAPRC_delH_SEMIvol/'//
     &           'SOAP/POA_DHVAP.txt', status='OLD')
            open (97, file='/home/bnmurphy/Research/PMCAMx/PMCAMxSAPRC_delH_SEMIvol/'//
     &           'SOAP/POA_LOGCSTAR.txt', status='OLD')
            open (98, file='/home/bnmurphy/Research/PMCAMx/PMCAMxSAPRC_delH_SEMIvol/'//
     &           'SOAP/POA_T.txt', status='OLD')
            open (99, file='/home/bnmurphy/Research/PMCAMx/PMCAMxSAPRC_delH_SEMIvol/'//
     &          'SOAP/SOA_DHVAP.txt', status='OLD')
            do tempcnt = 1,ntemp
                read(98, *), dhtemp(tempcnt)
                read(96, *), (poadhvap(tempcnt,icstar),icstar=1,ncstar)
                read(99, *), (soadhvap(tempcnt,icstar),icstar=1,ncstar)
            end do
            read(97, *), (dhcstar(icstar),icstar=1,ncstar)
            close(96)
            close(97)
            close(98)
            close(99)

c=============== End Hvap Table Open and Read =============


c
      return
c
 900  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(4(a,i4))')
     &  'reaction number', num, ' of type', rxntyp(i),
     &  ' should have', npar(rxntyp(i)), ' parameters'
      goto 910
c
 901  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(a)')  'For reaction type 5'
      write(iout,'(a,i4)') 'reference reaction # ', iref
      write(iout,'(a,i4)') 'must come before this reaction', i
      goto 910
c
 902  write(iout,'(//,a)') 'ERROR: Reading CHEMPARAM file record:'
      write(iout,'(a)') record
      write(iout,'(a)') 'Reaction type out of bounds (1-5)'
      write(iout,'(a,i4)') 'for reaction ', num
      write(iout,'(a,i4)') 'value was', rxntyp(i)
      write(iout,'(a)') 'Check this is CAMx3 chemparam file'
c
 910  write(*,'(/,a)') 'ERROR in READCHM - see message in .out file'
      write(idiag,'(/,a)') 'ERROR in READCHM - see message in .out file'
      write(iout,'(//,a)') 'ERROR in READCHM:'
      call camxerr()
c
      end
