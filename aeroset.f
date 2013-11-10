      subroutine aeroset(nsec_c,dsec_i,ierr)
c
c-----PMCAMx v3.01 020531
c
c     set-up routine for AERO routines in CAMx
c
c     Called by: 
c        READCHM
c 
      include 'camx.prm' 
      include 'chmstry.com' 
      include 'filunit.com'
      include 'section.inc'
      include 'camx_aero.inc'
      include 'camx.com'
      include 'grid.com'
      include 'section_aq.inc'
      include 'dbg.inc'
c
      real*8  dsec_i(nsecp1),dsecf_i(nsecp1)
c     
c-----Entry point 
c
      ierr=0
c
c--- Check nsec_c against nsec parameters used in AERO modules
c
	if ( nsec .ne. nsect ) then
	  write(iout,*) ' ERROR:  nsec and nsect must be equal!'
	  write(iout,*) ' NSEC =',nsec, '; NSECT =',nsect
	  ierr = 1
	  return
	endif
	if ( nsec_c .ne. nsec ) then
	  write(iout,*) ' ERROR:  nsec_c and nsec must be equal!'
	  write(iout,*) ' NSEC_C =',nsec_c, '; NSEC =',nsec
	  ierr = 1
	  return
	endif
c
c   The following specifies which AERO Modules to use 
c
      chaero='EQUI'
c      chaero='MADM'
c      chaero='HYBR'
c      chaq = 'RADM'
c      chaq = 'VSRM' ! DO NOT USE THIS !!!
      chaq = 'OVSR' ! USE THIS FOR CMU'S VSRM
      lsoap = .true.
      laq   = .true.
      laero = .true.
      ldbg_aero = .false.
      lgas  = .true.
c
c   set cwmin & tamin
c
      aqcwmin = cwmin
      aqtamin = tamin
c
c-- reset dtaero to hit output interval exactly, if necessary...
c
      dtio = amin1(60.,dtinp,dtems,dtout) ! [min]
      nsteps = INT( 0.999*dtio/dt_aero ) + 1
      dtaero = dtio/nsteps
      idate = begdate
      if ( idate .ge. 99999) then
	call juldate(idate)
      endif
      do n=1,ngrid
        date_aer(n) = idate
        grd_time(n) = begtim
        call uptime(grd_time(n),date_aer(n),60.*dtaero)
	aero_dt(n)  = 0.0
      enddo
      write(*,*) 'grd_time,date_aero: ',grd_time(1),date_aer(1)
      lfrst = .true.
      if (ldbg_aero) then
        write(iout,*) ' CHAERO :',chaero
        write(iout,*) ' laq,laero,lsoap :',laq,laero,lsoap
        write(*,*) ' CHAERO :',chaero
        write(*,*) ' laq,laero,lsoap :',laq,laero,lsoap
      endif
c
c    CALCULATE THE INITIAL DIAMETERS
c 
c
c  check user supplied section cutpoints -- monotonically increasing
c
      do i=2,nsecp1
        if ( dsec_i(i) .le. dsec_i(i-1) ) then
          write(iout,'(//,a)') 'ERROR in AEROSET:'
          write(iout,'(1X,a)')
     &'Invalid section cut-points ... must be monotonically increasing!'
          write(iout,'(1X,a,/)') 'Input cut-points are :'
          do n=1,nsecp1
            write(iout,'(4X,i2,1X,D9.3)') n,dsec_i(n)
          enddo
          call camxerr()
        endif
      enddo
c
c    For the MOVING sectional approach dsec is the SECTIONAL DIAMETER
      dmin=dsec_i(1)
      dmax=dsec_i(nsecp1)
cbk      do i=1,nsecp1
cbk        dsec_c(i)=dsec_i(i)
cbk      enddo
cbk      dfmin=dmin
cbk      dfmax=dmax
      do i=1,nsecp1
        dsecf_c(i)=dsec_i(i)
      enddo

      write(idiag,*) ' '
      write(idiag,*) 'Particle section cut-points:'
      do i=1,nsecp1
        write(idiag,'(1x,i2,1x,d9.3)') i,dsecf_c(i)
      enddo
      write(idiag,*) ' '

c     set moving diameters to logarithmic mean of fixed section diameters (tmg,01/25/02)
      do i=1,nsec
        dsec_c(i)=sqrt(dsecf_c(i)*dsecf_c(i+1))
      enddo
c
      kso2_c   = kso2
      if (idmech.eq.6) then
        kform_c  = kform
        kh2o2_c  = kh2o2
      elseif (idmech.eq.5) then
        kform_c  = khcho
        kh2o2_c  = kho2h
      endif
      khono_c  = khono
      ko3_c    = ko3
      koh_c    = koh
      kho2_c   = kho2
      kno3_c   = kno3
      kno_c    = kno
      kno2_c   = kno2
      kpan_c   = kpan
      kcpo1_c  = kcpo1
      kcpo2_c  = kcpo2
      kcpo3_c  = kcpo3
      kcpo4_c  = kcpo4
      kcpo5_c  = kcpo5
      kcpo6_c  = kcpo6
      kcpo7_c  = kcpo7
      kcpo8_c  = kcpo8
      kcoo1_c  = kcoo1
      kcoo2_c  = kcoo2
      kcoo3_c  = kcoo3
      kcoo4_c  = kcoo4
      kcoo5_c  = kcoo5
      kcoo6_c  = kcoo6
      kcoo7_c  = kcoo7
      kcoo8_c  = kcoo8
      kcbs1_c  = kcbs1
      kcbs2_c  = kcbs2
      kcbs3_c  = kcbs3
      kcbs4_c  = kcbs4
      kcbs5_c  = kcbs5
      kcas1_c  = kcas1
      kcas2_c  = kcas2
      kcas3_c  = kcas3
      kcas4_c  = kcas4
      kcas5_c  = kcas5
      kcns1_c  = kcns1
      kcns2_c  = kcns2
      kcns3_c  = kcns3
      kcns4_c  = kcns4
      kcns5_c  = kcns5
      kcns6_c  = kcns6
      kcns7_c  = kcns7
      kcns8_c  = kcns8
      khno3_c  = khno3
      knh3_c   = knh3
      kh2so4_c = ksulf
      khcl_c   = khcl
      kapo1_c  = kapo1_1
      kapo2_c  = kapo2_1
      kapo3_c  = kapo3_1
      kapo4_c  = kapo4_1
      kapo5_c  = kapo5_1
      kapo6_c  = kapo6_1
      kapo7_c  = kapo7_1
      kapo8_c  = kapo8_1
      kaoo1_c  = kaoo1_1
      kaoo2_c  = kaoo2_1
      kaoo3_c  = kaoo3_1
      kaoo4_c  = kaoo4_1
      kaoo5_c  = kaoo5_1
      kaoo6_c  = kaoo6_1
      kaoo7_c  = kaoo7_1
      kaoo8_c  = kaoo8_1
      kabs1_c  = kabs1_1
      kabs2_c  = kabs2_1
      kabs3_c  = kabs3_1
      kabs4_c  = kabs4_1
      kabs5_c  = kabs5_1
      kaas1_c  = kaas1_1
      kaas2_c  = kaas2_1
      kaas3_c  = kaas3_1
      kaas4_c  = kaas4_1
      kaas5_c  = kaas5_1
      kans1_c  = kans1_1
      kans2_c  = kans2_1
      kans3_c  = kans3_1
      kans4_c  = kans4_1
      kans5_c  = kans5_1
      kans6_c  = kans6_1
      kans7_c  = kans7_1
      kans8_c  = kans8_1
      kcrst_c  = kcrust_1
      kpoc_c   = kpoc_1
      kpec_c   = kpec_1
      kph2o_c  = kph2o_1
      kpcl_c   = kpcl_1
      kna_c    = kna_1
      kpnh4_c  = kpnh4_1
      kpno3_c  = kpno3_1
      kpso4_c  = kpso4_1
      knxoy_c  = knxoy        ! RADM

c
c  set wtfacs for interface between aerosol modules
c
cbk      do l=1,mxspec
cbk        wtfac_ae(l)=1.0
cbk        wtfac_aq(l)=1.0
cbk      enddo
c
cbk      do k=1,nsec
cbk        wtfac_ae(ksoa1_c+(k-1))=wtmol(ksoa1_c+(k-1))/150.
cbk        wtfac_ae(ksoa2_c+(k-1))=wtmol(ksoa2_c+(k-1))/150.
cbk        wtfac_ae(ksoa3_c+(k-1))=wtmol(ksoa3_c+(k-1))/150.
cbk        wtfac_ae(ksoa4_c+(k-1))=wtmol(ksoa4_c+(k-1))/180.
cbk        wtfac_ae(kpoc_c+(k-1)) =wtmol(kpoc_c+(k-1))/100.
cbk        wtfac_ae(kpec_c+(k-1)) =wtmol(kpec_c+(k-1))/100.
cbk        wtfac_ae(kcrst_c+(k-1))=wtmol(kcrst_c+(k-1))/100.
cbk        wtfac_ae(kph2o_c+(k-1))=wtmol(kph2o_c+(k-1))/18.
cbk        wtfac_ae(kpcl_c+(k-1)) =wtmol(kpcl_c+(k-1))/36.5
cbk        wtfac_ae(kna_c+(k-1))  =wtmol(kna_c+(k-1))/23.
cbk        wtfac_ae(kpnh4_c+(k-1))=wtmol(kpnh4_c+(k-1))/17.
cbk        wtfac_ae(kpno3_c+(k-1))=wtmol(kpno3_c+(k-1))/63.
cbk        wtfac_ae(kpso4_c+(k-1))=wtmol(kpso4_c+(k-1))/98.
c
cbk        wtfac_aq(ksoa1_c+(k-1))=wtmol(ksoa1_c+(k-1))/150.
cbk        wtfac_aq(ksoa2_c+(k-1))=wtmol(ksoa2_c+(k-1))/150.
cbk        wtfac_aq(ksoa3_c+(k-1))=wtmol(ksoa3_c+(k-1))/150.
cbk        wtfac_aq(ksoa4_c+(k-1))=wtmol(ksoa4_c+(k-1))/180.
cbk        wtfac_aq(kpoc_c+(k-1)) =wtmol(kpoc_c+(k-1))/100.
cbk        wtfac_aq(kpec_c+(k-1)) =wtmol(kpec_c+(k-1))/100.
cbk        wtfac_aq(kcrst_c+(k-1))=wtmol(kcrst_c+(k-1))/100.
cbk        wtfac_aq(kph2o_c+(k-1))=wtmol(kph2o_c+(k-1))/18.
cbk        wtfac_aq(kpcl_c+(k-1)) =wtmol(kpcl_c+(k-1))/35.5
cbk        wtfac_aq(kna_c+(k-1))  =wtmol(kna_c+(k-1))/23.
cbk        wtfac_aq(kpnh4_c+(k-1))=wtmol(kpnh4_c+(k-1))/18.
cbk        wtfac_aq(kpno3_c+(k-1))=wtmol(kpno3_c+(k-1))/62.
cbk        wtfac_aq(kpso4_c+(k-1))=wtmol(kpso4_c+(k-1))/96.
cbk      enddo
c
c      ----- End -----
c
      return
      end
