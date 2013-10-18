c-----PMCAMx v3.01 020531
c    include for passing aerosol module  parameters
c
      real*8 dsec_c(nsecp1),dsecf_c(nsecp1)
      real*8 dmin,dmax
c
      common /sections/ dmin,dmax,dsec_c,dsecf_c

c  Common aeroopt
c  contains option flags for the PM routines.
c
      character*4 chaero,chaq
      common /aeroopt/chaero,chaq
      real aqcwmin,aqtamin
      common /aqcond/aqcwmin,aqtamin
c   
      logical lfrst
      common /init/ lfrst
c
c  Common /ptrs/ contains pointer to CAMx species requird by PM routines
c
      common /ptrs/ kso2_c,kh2o2_c,kform_c,khono_c,ko3_c,koh_c,
     &              kho2_c,kno3_c,kno_c,kno2_c,kpan_c,
     & kcpo1_c,kcpo2_c,kcpo3_c,kcpo4_c,kcpo5_c,kcpo6_c,kcpo7_c,kcpo8_c,
     & kcpo9_c,kcpo10_c,kcoo1_c,kcoo2_c,kcoo3_c,kcoo4_c,kcoo5_c,
     & kcoo6_c,kcoo7_c,kcoo8_c,kcoo9_c,kcoo10_c,
     & kcbs1_c,kcbs2_c,kcbs3_c,kcbs4_c,
     & kcas1_c,kcas2_c,kcas3_c,kcas4_c,
     &              khno3_c,knh3_c,kh2so4_c,khcl_c,              
     & kapo1_c,kapo2_c,kapo3_c,kapo4_c,kapo5_c,kapo6_c,kapo7_c,kapo8_c,
     & kapo9_c,kapo10_c,kaoo1_c,kaoo2_c,kaoo3_c,kaoo4_c,kaoo5_c,
     & kaoo6_c,kaoo7_c,kaoo8_c,kaoo9_c,kaoo10_c,
     & kabs1_c,kabs2_c,kabs3_c,kabs4_c,
     & kaas1_c,kaas2_c,kaas3_c,kaas4_c,
     &             kcrst_c,kpoc_c,kpec_c,kph2o_c,kpcl_c,kna_c,
     &              kpnh4_c,kpno3_c,kpso4_c,knxoy_c
c
c  Common /wtfac/ contains adjustment factors for aerosol mol. wts. 
c   for interface between aerosol modules = CAMx mol.wts / Aerosol module mol.wts. 
c
c      real*4 wtfac_ae(MXSPEC), wtfac_aq(MXSPEC)
c
c      common /wtfac/ wtfac_ae,wtfac_aq
