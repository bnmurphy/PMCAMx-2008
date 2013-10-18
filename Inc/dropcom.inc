c
c COMMON GROUPS AND CORRESPONDING MATRICES FOR AQUEOUS-PHASE MODULE
c
	real*4 diameter(nbounds), dd(nsp_aq), daer(nsect)
	real*4 akeq(17),akhen(21),akre(120)
	real*4 wmol(29),amol(3),gmol(22)
	common / drop / diameter, dd, daer
	common / mw / wso2, wh2o2, whcho, whcooh, wnh3, whno3, whcl
	common /aqrate2/akeq,akhen,akre
	common /aqrate3/wmol,amol,gmol
