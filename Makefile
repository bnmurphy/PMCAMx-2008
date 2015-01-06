# Just type make PMCAMx.exe DOMAIN=string

#ifndef DOMAIN
DOMAIN = otag
#endif

INC   = ./Inc

CMC   = ./CMC

DDM   = ./DDM

PA    = ./PA

RTRAC = ./Rtrac

PIG   = ./PiG

BNRY  = ./IO_bin

OSAT  = ./Osat

AERO  = ./AERO

AER   = ./AER

SOAP  = ./SOAP

TARGT = CAMx

default:
	@echo '---------------------------------------------------------'
	@echo 'To make PMCAMx.exe type make PMCAMx.exe DOMAIN=string'
	@echo '----------------------------------------------------------'

# Pavan Nandan Racherla (pavanracherla@cmu.edu) removed the -ipo flag, which performs inter procedural optimizations.
# Date: Oct 1 2008.
ifort:
	@rm -f $(INC)/camx.prm
	@csh chktracer camx.prm.$(DOMAIN) PMCAMx.exe
	@ln -s camx.prm.$(DOMAIN) $(INC)/camx.prm
	make model FC="ifort" FLGS="-I$(INC) -O2 -fpe3 -traceback -align dcommons -extend_source -convert big_endian -mcmodel=medium -shared-intel" TARGT="PMCAMx.exe" DUM=dummy

pgf90:
	@rm -f $(INC)/camx.prm
	@csh chktracer camx.prm.$(DOMAIN) PMCAMx.exe
	@ln -s camx.prm.$(DOMAIN) $(INC)/camx.prm
	make model FC="pgf90" FLGS="-I$(INC) -O2 -tp k8-64 -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio -Wl, -mcmodel=medium" TARGT="PMCAMx.exe" DUM=dummy
#	make model FC="pgf90" FLGS="-I$(INC) -O2 -tp p6 -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio -Wl,-Bstatic" TARGT="PMCAMx.exe" DUM=dummy

clean:	
	rm -f $(OBJCTS) dummy*.o

superclean:	
	rm -f $(OBJCTS) dummy*.o *.exe *.mod

OBJCTS = \
CAMx.o \
aerochem.o \
aeroset.o \
aggdep.o \
aggr00.o \
aggreg.o \
ahoprep.o \
areag.o \
average.o \
bc1grd.o \
bcmodfy.o \
caldate.o \
camxerr.o \
chemdriv.o \
chemrxn.o \
chmdat.o \
chmprep.o \
cpivot.o \
cvtdate.o \
cvtwind.o \
dateerr.o \
depsmry.o \
diffus.o \
drydep.o \
$(DUM).o \
emassign.o \
emiss.o \
emistrns.o \
expndlay.o \
exptbl.o \
fgavrg.o \
finwind.o \
fullaero.o \
getdelt.o \
getznth.o \
grdprep.o \
pspgeo.o \
hadvbot.o \
hadvppm.o \
iassgn2d.o \
iehsolv.o \
iniptr.o \
interp2d.o \
interpv.o \
intrpcnc.o \
intrpdat.o \
istrln.o \
jstlft.o \
juldate.o \
khorz.o \
kphoto.o \
ktherm.o \
lcpgeo.o \
linpac.o \
lsode.o \
luassgn.o \
massum.o \
micromet.o \
nesting.o \
nstprep.o \
parntchd.o \
plumeris.o \
pntprep.o \
ProjUtils.o\
raddrivr.o \
rassgn3d.o \
rassgn4d.o \
readaho.o \
readchm.o \
readpht.o \
scavrat.o \
setbc1d.o \
setbc.o \
startup.o \
subdomain.o \
timestep.o \
timrates.o \
toupper.o \
trap.o \
trdiag.o \
updtmet.o \
uptime.o \
utmgeo.o \
vd_aer.o \
vd_gas.o \
vdiffimp.o \
vnmshcal.o \
vrtslv.o \
wetdep.o \
wrtmass.o \
xyadvec.o \
zadvec.o \
zeros.o \
zrates.o \
$(AERO)/hlconst.o \
$(AERO)/hlindex.o \
$(AERO)/isocom_v1.6.o \
$(AERO)/isofwd_v1.6.o \
$(AERO)/isorev_v1.6.o \
$(AERO)/raqchem.o \
$(BNRY)/areaprep.o \
$(BNRY)/bndprep.o \
$(BNRY)/clciwt.o \
$(BNRY)/cncprep.o \
$(BNRY)/depprep.o \
$(BNRY)/getdepth.o \
$(BNRY)/openfils.o \
$(BNRY)/rdfgcon.o \
$(BNRY)/rdpthdr.o \
$(BNRY)/rdsumbc.o \
$(BNRY)/readar.o \
$(BNRY)/readbnd.o \
$(BNRY)/readcnc.o \
$(BNRY)/readinp.o \
$(BNRY)/readpt.o \
$(BNRY)/readzpwt.o \
$(BNRY)/srfprep.o \
$(BNRY)/wrfgcon.o \
$(BNRY)/wrfgdep.o \
$(BNRY)/wrtcon.o \
$(BNRY)/wrtdep.o \
$(CMC)/ddmjac1.o \
$(CMC)/ddmjac2.o \
$(CMC)/ddmjac3.o \
$(CMC)/ddmjac4.o \
$(CMC)/ddmjac5.o \
$(CMC)/ddmjac6.o \
$(CMC)/iejac1.o \
$(CMC)/iejac2.o \
$(CMC)/iejac3.o \
$(CMC)/iejac4.o \
$(CMC)/iejac5.o \
$(CMC)/iejac6.o \
$(CMC)/ierate1.o \
$(CMC)/ierate2.o \
$(CMC)/ierate3.o \
$(CMC)/ierate4.o \
$(CMC)/ierate5.o \
$(CMC)/ierate6.o \
$(CMC)/ierxn1.o \
$(CMC)/ierxn2.o \
$(CMC)/ierxn3.o \
$(CMC)/ierxn4.o \
$(CMC)/ierxn5.o \
$(CMC)/ierxn6.o \
$(CMC)/ieslow1.o \
$(CMC)/ieslow2.o \
$(CMC)/ieslow3.o \
$(CMC)/ieslow4.o \
$(CMC)/ieslow5.o \
$(CMC)/ieslow6.o \
$(CMC)/radslvr1.o \
$(CMC)/radslvr2.o \
$(CMC)/radslvr3.o \
$(CMC)/radslvr4.o \
$(CMC)/radslvr5.o \
$(CMC)/radslvr6.o \
$(CMC)/ratejac1.o \
$(CMC)/ratejac2.o \
$(CMC)/ratejac3.o \
$(CMC)/ratejac4.o \
$(CMC)/ratejac5.o \
$(CMC)/ratejac6.o \
$(CMC)/rateslo1.o \
$(CMC)/rateslo2.o \
$(CMC)/rateslo3.o \
$(CMC)/rateslo4.o \
$(CMC)/rateslo5.o \
$(CMC)/rateslo6.o \
$(CMC)/rxnrate1.o \
$(CMC)/rxnrate2.o \
$(CMC)/rxnrate3.o \
$(CMC)/rxnrate4.o \
$(CMC)/rxnrate5.o \
$(CMC)/rxnrate6.o \
$(DDM)/adjddmc0.o \
$(DDM)/avgrcpddm.o \
$(DDM)/bottddm.o \
$(DDM)/clrbdyddm.o \
$(DDM)/cvticddm.o \
$(DDM)/ddmchem.o \
$(DDM)/filspddm.o \
$(DDM)/hdrrcpddm.o \
$(DDM)/loaddm.o \
$(DDM)/rdarddm.o \
$(DDM)/rdbcddm.o \
$(DDM)/rdicddm.o \
$(DDM)/rdoptddm.o \
$(DDM)/rdptddm.o \
$(DDM)/specddm.o \
$(DDM)/startddm.o \
$(OSAT)/addrcp.o \
$(OSAT)/adjstc0.o \
$(OSAT)/adjstsa.o \
$(OSAT)/apcasa.o \
$(OSAT)/avgrcp.o \
$(OSAT)/avgwal.o \
$(OSAT)/clcbwt.o \
$(OSAT)/clcewt.o \
$(OSAT)/clrbdysa.o \
$(OSAT)/emprepsa.o \
$(OSAT)/filaqsa.o \
$(OSAT)/filbdysa.o \
$(OSAT)/fillar.o \
$(OSAT)/fillpt.o \
$(OSAT)/filvdsa.o \
$(OSAT)/goatsa.o \
$(OSAT)/hdrrcp.o \
$(OSAT)/hdrwsa.o \
$(OSAT)/initsa.o \
$(OSAT)/instsa.o \
$(OSAT)/osatsa.o \
$(OSAT)/pekrcp.o \
$(OSAT)/pigsa.o \
$(OSAT)/rdargrp.o \
$(OSAT)/rdfgsa.o \
$(OSAT)/rdinstsa.o \
$(OSAT)/rdoptsa.o \
$(OSAT)/rdptgrp.o \
$(OSAT)/readarsa.o \
$(OSAT)/readptsa.o \
$(OSAT)/recalib.o \
$(OSAT)/rercp.o \
$(OSAT)/resmap.o \
$(OSAT)/specsa.o \
$(OSAT)/stabsa.o \
$(OSAT)/startsa.o \
$(OSAT)/sum1grd.o \
$(OSAT)/sum1pnt.o \
$(OSAT)/sumgrps.o \
$(OSAT)/sumicwt.o \
$(OSAT)/sumwt4.o \
$(OSAT)/timadv.o \
$(OSAT)/wrfgsa.o \
$(OSAT)/xfluxsa.o \
$(OSAT)/yfluxsa.o \
$(PA)/cpamech3.o \
$(PA)/cpamech5.o \
$(PA)/cparad.o \
$(PA)/initipr.o \
$(PA)/pagrids.o \
$(PA)/pasetup.o \
$(PA)/pazero.o \
$(PA)/rdoptpa.o \
$(PA)/wrtcgcpa.o \
$(PA)/wrtfgcpa.o \
$(PA)/wrtipr.o \
$(PA)/wrtiprhdr.o \
$(PA)/wrtirr.o \
$(PA)/wrtirrhdr.o \
$(PIG)/avepig.o \
$(PIG)/greschem.o \
$(PIG)/gresdriv.o \
$(PIG)/gresgrow.o \
$(PIG)/gresmscl.o \
$(PIG)/irondriv.o \
$(PIG)/ironmscl.o \
$(PIG)/piginit.o \
$(PIG)/pigprep.o \
$(PIG)/pigwalk.o \
$(PIG)/walk1pig.o \
$(PIG)/wrtpig.o \
$(RTRAC)/chemrt.o \
$(RTRAC)/cvticrt.o \
$(RTRAC)/drydeprt.o \
$(RTRAC)/empreprt.o \
$(RTRAC)/hdrcprt.o \
$(RTRAC)/rdarrt.o \
$(RTRAC)/rdbcrt.o \
$(RTRAC)/rdchmrt.o \
$(RTRAC)/rdicrt.o \
$(RTRAC)/rdptrt.o \
$(RTRAC)/rdrcprt.o \
$(RTRAC)/startrt.o \
$(RTRAC)/wetdeprt.o \
$(RTRAC)/wrrcprt.o \
$(SOAP)/soapdat.o \
$(SOAP)/soap.o \
$(SOAP)/spfcn.o \
$(AER)/addit.o \
$(AER)/aerchem.o \
$(AER)/aero_err.o \
$(AER)/aqchem.o \
$(AER)/aqdist.o \
$(AER)/aqfex1.o \
$(AER)/aqfex2.o \
$(AER)/aqintegr1.o \
$(AER)/aqintegr2.o \
$(AER)/aqjex.o \
$(AER)/aqoperator1.o \
$(AER)/aqoperator2.o \
$(AER)/aqrates1.o \
$(AER)/aqrates2.o \
$(AER)/baddit.o \
$(AER)/block.o \
$(AER)/bmass.o \
$(AER)/coagul.o \
$(AER)/constants.o \
$(AER)/decisions.o \
$(AER)/diameter.o \
$(AER)/diff.o \
$(AER)/differ.o \
$(AER)/diffund.o \
$(AER)/dropinit.o \
$(AER)/electro.o \
$(AER)/eqpart.o \
$(AER)/eqparto.o \
$(AER)/equaer.o \
$(AER)/fullequil.o \
$(AER)/hybrid.o \
$(AER)/linint.o \
$(AER)/madm.o \
$(AER)/mass.o \
$(AER)/negchk.o \
$(AER)/newdist.o \
$(AER)/nucl.o \
$(AER)/qsaturation.o \
$(AER)/react.o \
$(AER)/state.o \
$(AER)/steady.o \
$(AER)/step.o \
$(AER)/svode.o \
$(AER)/values.o \
$(AER)/vsrm.o \
$(AER)/slsode.o

model: 	$(OBJCTS)
	$(FC) -o $(TARGT) $(FLGS) $(OBJCTS) $(LIBS)

# Modifications by Pavan Nandan Racherla (July 2 2008)--->
# Date: July 2 2008.

# Suffixes and instructions on how to compile those files:
.SUFFIXES : .f90 .f .mod .o
.f90.o :
	$(FC) -c -o $@ $(FLGS) $<
.f.o :
	$(FC) -c -o $@ $(FLGS) $<

# Build module files using a 'touch' command:
# 'touch' updates the access time and modification time/dates to the current time and date
# if the .mod file does'nt exist, 'touch' creates it with a filesize of 0
%.mod :
	@if [ ! -s $@ ] ; then \
	echo "error, module file deleted?"; rm $^; exit 1; fi
	echo "module file $@ considered updated through dependence on $^"
	touch $@

# module-object dependencies:
projutils.mod : ProjUtils.o

# Instructions to create objects wrt modules (based on the USE statements)
pspgeo.o : projutils.mod

# End of my modifications


CAMx.o 			: CAMx.f                                               \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/grid.com $(INC)/ptemiss.com $(INC)/bndary.com   \
                        $(INC)/chmstry.com $(INC)/flags.com $(INC)/ahomap.com  \
                        $(INC)/filunit.com $(INC)/tracer.com                   \
                        $(INC)/rtracchm.com $(INC)/procan.com

aerochem.o 		: aerochem.f                                           \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/chmstry.com     \
                        $(INC)/filunit.com $(INC)/procan.com

aeroset.o		: aeroset.f                                            \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com  \
                        $(INC)/section.inc $(INC)/camx_aero.inc                \
                        $(INC)/camx.com $(INC)/grid.com $(INC)/section_aq.inc  \
                        $(INC)/dbg.inc

aggr00.o 		: aggr00.f                                             \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/grid.com $(INC)/chmstry.com $(INC)/filunit.com  \
                        $(INC)/tracer.com $(INC)/procan.com

ahoprep.o 		: ahoprep.f                                            \
                        $(INC)/camx.prm $(INC)/ahomap.com $(INC)/filunit.com   \
                        $(INC)/grid.com

areag.o			: areag.f

average.o 		: average.f                                            \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com      \
                        $(INC)/chmstry.com $(INC)/tracer.com                   \
                        $(INC)/rtracchm.com $(INC)/procan.com

bcmodfy.o 		: bcmodfy.f                                            \
                        $(INC)/camx.prm $(INC)/filunit.com

camxerr.o 		: camxerr.f                                            \
                        $(INC)/camx.prm $(INC)/filunit.com

chemdriv.o 		: chemdriv.f                                           \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com      \
                        $(INC)/chmstry.com $(INC)/filunit.com                  \
                        $(INC)/ahomap.com $(INC)/procan.com                    \
                        $(INC)/tracer.com

chemrxn.o 		: chemrxn.f                                            \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/grid.com $(INC)/flags.com $(INC)/filunit.com    \
                        $(INC)/chmstry.com $(INC)/procan.com

chmdat.o 		: chmdat.f                                             \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/camx.com     \
                        $(INC)/deposit.com

chmprep.o 		: chmprep.f                                            \
                        $(INC)/camx.prm $(INC)/flags.com $(INC)/tracer.com

dateerr.o 		: dateerr.f                                            \
                        $(INC)/camx.prm $(INC)/filunit.com

depsmry.o 		: depsmry.f                                            \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/filunit.com

diffus.o 		: diffus.f                                             \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/bndary.com   \
                        $(INC)/chmstry.com $(INC)/tracer.com $(INC)/procan.com

drydep.o 		: drydep.f                                             \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com      \
                        $(INC)/deposit.com $(INC)/chmstry.com                  \
                        $(INC)/filunit.com $(INC)/section.inc

emiss.o 		: emiss.f                                              \
                        $(INC)/camx.prm $(INC)/ptemiss.com $(INC)/bndary.com   \
                        $(INC)/flags.com $(INC)/procan.com $(INC)/tracer.com

emistrns.o 		: emistrns.f                                           \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/chmstry.com $(INC)/grid.com $(INC)/flags.com    \
                        $(INC)/filunit.com $(INC)/ptemiss.com                  \
                        $(INC)/bndary.com $(INC)/procan.com                    \
                        $(INC)/tracer.com $(INC)/rtracchm.com

exptbl.o 		: exptbl.f                                             \
                        $(INC)/camx.prm $(INC)/chmstry.com

fgavrg.o 		: fgavrg.f                                             \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/grid.com $(INC)/chmstry.com $(INC)/flags.com    \
                        $(INC)/tracer.com $(INC)/procan.com

fullaero.o		: fullaero.f                                           \
                        $(INC)/aerpar.inc $(INC)/dynamic.inc                   \
                        $(INC)/droppar.inc $(INC)/camx_aero.inc $(INC)/dbg.inc

getdelt.o 		: getdelt.f                                            \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com

# Pavan:
grdprep.o		: grdprep.f                                            \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/flags.com

hadvbot.o 		: hadvbot.f                                            \
                        $(INC)/camx.prm

hadvppm.o 		: hadvppm.f                                            \
                        $(INC)/camx.prm

iehsolv.o 		: iehsolv.f                                            \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

iniptr.o 		: iniptr.f                                             \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/filunit.com     \
                        $(INC)/tracer.com

interpv.o 		: interpv.f                                            \
                        $(INC)/camx.prm

intrpdat.o 		: intrpdat.f                                           \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/filunit.com $(INC)/flags.com $(INC)/grid.com

khorz.o 		: khorz.f                                              \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/filunit.com

kphoto.o 		: kphoto.f                                             \
                        $(INC)/camx.prm $(INC)/chmstry.com

ktherm.o 		: ktherm.f                                             \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

massum.o 		: massum.f                                             \
                        $(INC)/camx.prm $(INC)/bndary.com

nesting.o 		: nesting.f                                            \
                        $(INC)/camx.prm $(INC)/grid.com

nstprep.o 		: nstprep.f                                            \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/filunit.com     \
                        $(INC)/tracer.com

parntchd.o 		: parntchd.f                                           \
                        $(INC)/camx.prm $(INC)/filunit.com

pntprep.o 		: pntprep.f                                            \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/chmstry.com     \
                        $(INC)/ptemiss.com $(INC)/filunit.com $(INC)/flags.com \
                        $(INC)/bndary.com

raddrivr.o 		: raddrivr.f                                           \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com      \
                        $(INC)/camxfld.com $(INC)/grid.com $(INC)/chmstry.com
 
readaho.o 		: readaho.f                                            \
                        $(INC)/camx.prm $(INC)/filunit.com

readchm.o 		: readchm.f                                            \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com  \
                        $(INC)/flags.com $(INC)/ddmchm.com $(INC)/iehchem.com

readpht.o 		: readpht.f                                            \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com  \
                        $(INC)/ahomap.com

setbc1d.o 		: setbc1d.f                                            \
                        $(INC)/camx.prm

setbc.o 		: setbc.f                                              \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/camxfld.com     \
                        $(INC)/chmstry.com $(INC)/tracer.com

startup.o 		: startup.f                                            \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/filunit.com $(INC)/chmstry.com                  \
                        $(INC)/ahomap.com $(INC)/grid.com $(INC)/pigsty.com    \
                        $(INC)/flags.com $(INC)/tracer.com $(INC)/procan.com

subdomain.o		: subdomain.f

timestep.o 		: timestep.f                                           \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/bndary.com $(INC)/grid.com

timrates.o 		: timrates.f                                           \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com

trap.o 			: trap.f                                               \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com
 
updtmet.o 		: updtmet.f                                            \
                        $(INC)/camx.prm $(INC)/bndary.com

vdiffimp.o 		: vdiffimp.f                                           \
                        $(INC)/camx.prm $(INC)/procan.com

vnmshcal.o 		: vnmshcal.f                                           \
                        $(INC)/camx.prm $(INC)/filunit.com

vrtslv.o 		: vrtslv.f                                             \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/tracer.com

wetdep.o 		: wetdep.f                                             \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/tracer.com $(INC)/procan.com                    \
                        $(INC)/deposit.com $(INC)/section.inc

wrtmass.o 		: wrtmass.f                                            \
                        $(INC)/camx.prm $(INC)/camxfld.com $(INC)/chmstry.com  \
                        $(INC)/filunit.com $(INC)/flags.com

xyadvec.o 		: xyadvec.f                                            \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/filunit.com $(INC)/flags.com $(INC)/tracer.com  \
                        $(INC)/procan.com

zadvec.o 		: zadvec.f                                             \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/filunit.com $(INC)/tracer.com                   \
                        $(INC)/rtracchm.com $(INC)/procan.com

zrates.o 		: zrates.f                                             \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/flags.com

$(AERO)/isocom_v1.6.o 	: $(AERO)/isocom_v1.6.f	                               \
                        $(INC)/isrpia.inc

$(AERO)/isofwd_v1.6.o 	: $(AERO)/isofwd_v1.6.f	                               \
                        $(INC)/isrpia.inc

$(AERO)/isorev_v1.6.o 	: $(AERO)/isorev_v1.6.f	                               \
                        $(INC)/isrpia.inc

$(BNRY)/areaprep.o 	: $(BNRY)/areaprep.f                                   \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/chmstry.com     \
                        $(INC)/flags.com

$(BNRY)/bndprep.o 	: $(BNRY)/bndprep.f                                    \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/bndary.com   \
                        $(INC)/grid.com $(INC)/chmstry.com $(INC)/flags.com

$(BNRY)/clciwt.o 	: $(BNRY)/clciwt.f                                     \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/chmstry.com  \
                        $(INC)/bndary.com $(INC)/grid.com $(INC)/flags.com     \
                        $(INC)/tracer.com

$(BNRY)/cncprep.o 	: $(BNRY)/cncprep.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/filunit.com $(INC)/grid.com $(INC)/chmstry.com  \
                        $(INC)/flags.com

$(BNRY)/depprep.o 	: $(BNRY)/depprep.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/filunit.com $(INC)/grid.com $(INC)/chmstry.com  \
                        $(INC)/flags.com

$(BNRY)/getdepth.o 	: $(BNRY)/getdepth.f                                   \
                        $(INC)/camx.prm $(INC)/filunit.com

$(BNRY)/openfils.o 	: $(BNRY)/openfils.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/filunit.com $(INC)/chmstry.com                  \
                        $(INC)/ahomap.com $(INC)/grid.com $(INC)/flags.com     \
                        $(INC)/tracer.com $(INC)/procan.com

$(BNRY)/rdfgcon.o 	: $(BNRY)/rdfgcon.f                                    \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/camxfld.com  \
                        $(INC)/grid.com $(INC)/chmstry.com

$(BNRY)/rdpthdr.o 	: $(BNRY)/rdpthdr.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/chmstry.com     \
                        $(INC)/ptemiss.com $(INC)/filunit.com $(INC)/flags.com \
                        $(INC)/bndary.com $(INC)/tracer.com

$(BNRY)/rdsumbc.o 	: $(BNRY)/rdsumbc.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/chmstry.com $(INC)/bndary.com $(INC)/tracer.com

$(BNRY)/readar.o 	: $(BNRY)/readar.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/chmstry.com     \
                        $(INC)/flags.com

$(BNRY)/readbnd.o 	: $(BNRY)/readbnd.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/filunit.com $(INC)/bndary.com $(INC)/grid.com   \
                        $(INC)/chmstry.com

$(BNRY)/readcnc.o 	: $(BNRY)/readcnc.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/camxfld.com     \
                        $(INC)/filunit.com $(INC)/bndary.com $(INC)/grid.com   \
                        $(INC)/chmstry.com $(INC)/flags.com

$(BNRY)/readinp.o 	: $(BNRY)/readinp.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/grid.com $(INC)/camxfld.com

$(BNRY)/readpt.o 	: $(BNRY)/readpt.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/ptemiss.com $(INC)/chmstry.com $(INC)/flags.com \
                        $(INC)/grid.com

$(BNRY)/readzpwt.o 	: $(BNRY)/readzpwt.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/flags.com

$(BNRY)/srfprep.o 	: $(BNRY)/srfprep.f                                    \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/bndary.com

$(BNRY)/wrfgcon.o 	: $(BNRY)/wrfgcon.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/flags.com       \
                        $(INC)/filunit.com $(INC)/camxfld.com                  \
                        $(INC)/chmstry.com $(INC)/grid.com

$(BNRY)/wrfgdep.o 	: $(BNRY)/wrfgdep.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/camxfld.com $(INC)/chmstry.com $(INC)/grid.com

$(BNRY)/wrtcon.o 	: $(BNRY)/wrtcon.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/grid.com        \
                        $(INC)/chmstry.com $(INC)/flags.com

$(BNRY)/wrtdep.o 	: $(BNRY)/wrtdep.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/grid.com        \
                        $(INC)/chmstry.com

$(CMC)/ddmjac1.o 	: $(CMC)/ddmjac1.f                                     \
                        $(INC)/camx.prm $(INC)/ddmchm.com $(INC)/filunit.com

$(CMC)/ddmjac2.o 	: $(CMC)/ddmjac2.f                                     \
                        $(INC)/camx.prm $(INC)/ddmchm.com $(INC)/filunit.com

$(CMC)/ddmjac3.o 	: $(CMC)/ddmjac3.f                                     \
                        $(INC)/camx.prm $(INC)/ddmchm.com $(INC)/filunit.com

$(CMC)/ddmjac4.o 	: $(CMC)/ddmjac4.f                                     \
                        $(INC)/camx.prm $(INC)/ddmchm.com $(INC)/filunit.com

$(CMC)/ddmjac5.o 	: $(CMC)/ddmjac5.f                                     \
                        $(INC)/camx.prm $(INC)/ddmchm.com $(INC)/filunit.com

$(CMC)/ddmjac6.o 	: $(CMC)/ddmjac6.f                                     \
                        $(INC)/camx.prm $(INC)/ddmchm.com $(INC)/filunit.com

$(CMC)/iejac1.o 	: $(CMC)/iejac1.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/iejac2.o 	: $(CMC)/iejac2.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/iejac3.o 	: $(CMC)/iejac3.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/iejac4.o 	: $(CMC)/iejac4.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/iejac5.o 	: $(CMC)/iejac5.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/iejac6.o 	: $(CMC)/iejac6.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierate1.o 	: $(CMC)/ierate1.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierate2.o 	: $(CMC)/ierate2.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierate3.o 	: $(CMC)/ierate3.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierate4.o 	: $(CMC)/ierate4.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierate5.o 	: $(CMC)/ierate5.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierate6.o 	: $(CMC)/ierate6.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierxn1.o 	: $(CMC)/ierxn1.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierxn2.o 	: $(CMC)/ierxn2.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierxn3.o 	: $(CMC)/ierxn3.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierxn4.o 	: $(CMC)/ierxn4.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierxn5.o 	: $(CMC)/ierxn5.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ierxn6.o 	: $(CMC)/ierxn6.f                                      \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ieslow1.o 	: $(CMC)/ieslow1.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ieslow2.o 	: $(CMC)/ieslow2.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ieslow3.o 	: $(CMC)/ieslow3.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ieslow4.o 	: $(CMC)/ieslow4.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ieslow5.o 	: $(CMC)/ieslow5.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/ieslow6.o 	: $(CMC)/ieslow6.f                                     \
                        $(INC)/camx.prm $(INC)/iehchem.com

$(CMC)/radslvr1.o 	: $(CMC)/radslvr1.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

$(CMC)/radslvr2.o 	: $(CMC)/radslvr2.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

$(CMC)/radslvr3.o 	: $(CMC)/radslvr3.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

$(CMC)/radslvr4.o 	: $(CMC)/radslvr4.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

$(CMC)/radslvr5.o 	: $(CMC)/radslvr5.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

$(CMC)/radslvr6.o 	: $(CMC)/radslvr6.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com

$(CMC)/ratejac1.o 	: $(CMC)/ratejac1.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/ratejac2.o 	: $(CMC)/ratejac2.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/ratejac3.o 	: $(CMC)/ratejac3.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/ratejac4.o 	: $(CMC)/ratejac4.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/ratejac5.o 	: $(CMC)/ratejac5.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/ratejac6.o 	: $(CMC)/ratejac6.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rateslo1.o 	: $(CMC)/rateslo1.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rateslo2.o 	: $(CMC)/rateslo2.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rateslo3.o 	: $(CMC)/rateslo3.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rateslo4.o 	: $(CMC)/rateslo4.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rateslo5.o 	: $(CMC)/rateslo5.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rateslo6.o 	: $(CMC)/rateslo6.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rxnrate1.o 	: $(CMC)/rxnrate1.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rxnrate2.o 	: $(CMC)/rxnrate2.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rxnrate3.o 	: $(CMC)/rxnrate3.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rxnrate4.o 	: $(CMC)/rxnrate4.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rxnrate5.o 	: $(CMC)/rxnrate5.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(CMC)/rxnrate6.o 	: $(CMC)/rxnrate6.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com

$(DDM)/adjddmc0.o 	: $(DDM)/adjddmc0.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(DDM)/avgrcpddm.o 	: $(DDM)/avgrcpddm.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/tracer.com      \
                        $(INC)/filunit.com

$(DDM)/bottddm.o 	: $(DDM)/bottddm.f                                     \
                        $(INC)/camx.prm $(INC)/tracer.com

$(DDM)/clrbdyddm.o 	: $(DDM)/clrbdyddm.f                                   \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/tracer.com

$(DDM)/cvticddm.o 	: $(DDM)/cvticddm.f                                    \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/camx.com

$(DDM)/filspddm.o 	: $(DDM)/filspddm.f                                    \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/chmstry.com  \
                        $(INC)/tracer.com

$(DDM)/hdrrcpddm.o 	: $(DDM)/hdrrcpddm.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/chmstry.com     \
                        $(INC)/tracer.com $(INC)/filunit.com

$(DDM)/loaddm.o 	: $(DDM)/loaddm.f                                      \
                        $(INC)/camx.prm $(INC)/tracer.com

$(DDM)/rdarddm.o 	: $(DDM)/rdarddm.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/chmstry.com $(INC)/grid.com    \
                        $(INC)/bndary.com $(INC)/tracer.com

$(DDM)/rdbcddm.o 	: $(DDM)/rdbcddm.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/chmstry.com $(INC)/grid.com $(INC)/flags.com    \
                        $(INC)/bndary.com $(INC)/tracer.com

$(DDM)/rdicddm.o 	: $(DDM)/rdicddm.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/chmstry.com $(INC)/grid.com $(INC)/flags.com    \
                        $(INC)/bndary.com $(INC)/tracer.com

$(DDM)/rdoptddm.o 	: $(DDM)/rdoptddm.f                                    \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/tracer.com

$(DDM)/rdptddm.o 	: $(DDM)/rdptddm.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/chmstry.com $(INC)/grid.com    \
                        $(INC)/bndary.com $(INC)/tracer.com

$(DDM)/specddm.o 	: $(DDM)/specddm.f                                     \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/chmstry.com  \
                        $(INC)/tracer.com

$(DDM)/startddm.o 	: $(DDM)/startddm.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/flags.com       \
                        $(INC)/tracer.com

$(OSAT)/addrcp.o 	: $(OSAT)/addrcp.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com 

$(OSAT)/adjstc0.o 	: $(OSAT)/adjstc0.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(OSAT)/adjstsa.o 	: $(OSAT)/adjstsa.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(OSAT)/apcasa.o 	: $(OSAT)/apcasa.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(OSAT)/avgrcp.o 	: $(OSAT)/avgrcp.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/tracer.com      \
                        $(INC)/filunit.com

$(OSAT)/avgwal.o 	: $(OSAT)/avgwal.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com      \
                        $(INC)/tracer.com

$(OSAT)/clcbwt.o 	: $(OSAT)/clcbwt.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/chmstry.com $(INC)/bndary.com $(INC)/tracer.com

$(OSAT)/clcewt.o 	: $(OSAT)/clcewt.f                                     \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/bndary.com   \
                        $(INC)/grid.com $(INC)/chmstry.com $(INC)/tracer.com

$(OSAT)/clrbdysa.o 	: $(OSAT)/clrbdysa.f                                   \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/tracer.com

$(OSAT)/emprepsa.o 	: $(OSAT)/emprepsa.f                                   \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/grid.com     \
                        $(INC)/flags.com $(INC)/chmstry.com $(INC)/tracer.com

$(OSAT)/filaqsa.o 	: $(OSAT)/filaqsa.f                                    \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/tracer.com

$(OSAT)/filbdysa.o 	: $(OSAT)/filbdysa.f                                   \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/tracer.com

$(OSAT)/filvdsa.o 	: $(OSAT)/filvdsa.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/tracer.com

$(OSAT)/goatsa.o 	: $(OSAT)/goatsa.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(OSAT)/hdrrcp.o 	: $(OSAT)/hdrrcp.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/tracer.com      \
                        $(INC)/filunit.com

$(OSAT)/hdrwsa.o 	: $(OSAT)/hdrwsa.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/grid.com        \
                        $(INC)/flags.com $(INC)/tracer.com $(INC)/filunit.com

$(OSAT)/initsa.o 	: $(OSAT)/initsa.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(OSAT)/instsa.o 	: $(OSAT)/instsa.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/tracer.com      \
                        $(INC)/filunit.com

$(OSAT)/osatsa.o 	: $(OSAT)/osatsa.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(OSAT)/pekrcp.o 	: $(OSAT)/pekrcp.f                                     \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/tracer.com

$(OSAT)/pigsa.o 	: $(OSAT)/pigsa.f                                      \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/pigsty.com      \
                        $(INC)/tracer.com

$(OSAT)/rdargrp.o 	: $(OSAT)/rdargrp.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/grid.com $(INC)/bndary.com     \
                        $(INC)/tracer.com

$(OSAT)/rdfgsa.o 	: $(OSAT)/rdfgsa.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com      \
                        $(INC)/filunit.com

$(OSAT)/rdinstsa.o 	: $(OSAT)/rdinstsa.f                                   \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/flags.com       \
                        $(INC)/tracer.com $(INC)/filunit.com

$(OSAT)/rdoptsa.o 	: $(OSAT)/rdoptsa.f                                    \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/grid.com     \
                        $(INC)/tracer.com

$(OSAT)/rdptgrp.o 	: $(OSAT)/rdptgrp.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/grid.com $(INC)/ptemiss.com    \
                        $(INC)/tracer.com

$(OSAT)/readarsa.o 	: $(OSAT)/readarsa.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/grid.com        \
                        $(INC)/bndary.com $(INC)/tracer.com $(INC)/filunit.com

$(OSAT)/readptsa.o 	: $(OSAT)/readptsa.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/flags.com       \
                        $(INC)/grid.com $(INC)/ptemiss.com $(INC)/tracer.com   \
                        $(INC)/filunit.com

$(OSAT)/recalib.o 	: $(OSAT)/recalib.f                                    \
                        $(INC)/camx.prm $(INC)/tracer.com

$(OSAT)/rercp.o 	: $(OSAT)/rercp.f                                      \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com      \
                        $(INC)/filunit.com

$(OSAT)/resmap.o 	: $(OSAT)/resmap.f                                     \
                        $(INC)/camx.prm $(INC)/tracer.com $(INC)/filunit.com   \
                        $(INC)/grid.com

$(OSAT)/specsa.o 	: $(OSAT)/specsa.f                                     \
                        $(INC)/camx.prm $(INC)/tracer.com $(INC)/filunit.com

$(OSAT)/stabsa.o 	: $(OSAT)/stabsa.f                                     \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/tracer.com

$(OSAT)/startsa.o 	: $(OSAT)/startsa.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/flags.com       \
                        $(INC)/tracer.com

$(OSAT)/sum1grd.o 	: $(OSAT)/sum1grd.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com

$(OSAT)/sum1pnt.o 	: $(OSAT)/sum1pnt.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/flags.com       \
                        $(INC)/tracer.com

$(OSAT)/sumgrps.o 	: $(OSAT)/sumgrps.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/bndary.com $(INC)/grid.com     \
                        $(INC)/chmstry.com $(INC)/tracer.com

$(OSAT)/sumwt4.o 	: $(OSAT)/sumwt4.f                                     \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/tracer.com

$(OSAT)/timadv.o 	: $(OSAT)/timadv.f                                     \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/filunit.com $(INC)/flags.com $(INC)/tracer.com

$(OSAT)/wrfgsa.o 	: $(OSAT)/wrfgsa.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/grid.com        \
                        $(INC)/tracer.com $(INC)/filunit.com

$(OSAT)/xfluxsa.o 	: $(OSAT)/xfluxsa.f                                    \
                        $(INC)/camx.prm $(INC)/tracer.com

$(OSAT)/yfluxsa.o 	: $(OSAT)/yfluxsa.f                                    \
                        $(INC)/camx.prm $(INC)/tracer.com

$(PA)/cpamech3.o 	: $(PA)/cpamech3.f                                     \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/tracer.com   \
                        $(INC)/procan.com

$(PA)/cpamech5.o 	: $(PA)/cpamech5.f                                     \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/tracer.com   \
                        $(INC)/procan.com

$(PA)/cparad.o 		: $(PA)/cparad.f                                       \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/tracer.com $(INC)/procan.com $(INC)/chmstry.com

$(PA)/initipr.o 	: $(PA)/initipr.f                                      \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/procan.com

$(PA)/pagrids.o 	: $(PA)/pagrids.f                                      \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/grid.com     \
                        $(INC)/procan.com

$(PA)/pasetup.o 	: $(PA)/pasetup.f                                      \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/grid.com     \
                        $(INC)/tracer.com $(INC)/filunit.com $(INC)/procan.com

$(PA)/pazero.o 		: $(PA)/pazero.f                                       \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/chmstry.com     \
                        $(INC)/tracer.com $(INC)/procan.com

$(PA)/rdoptpa.o 	: $(PA)/rdoptpa.f                                      \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/grid.com     \
                        $(INC)/procan.com

$(PA)/wrtcgcpa.o 	: $(PA)/wrtcgcpa.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/tracer.com      \
                        $(INC)/filunit.com $(INC)/flags.com

$(PA)/wrtfgcpa.o 	: $(PA)/wrtfgcpa.f                                     \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/flags.com       \
                        $(INC)/filunit.com $(INC)/camxfld.com                  \
                        $(INC)/chmstry.com $(INC)/grid.com $(INC)/tracer.com

$(PA)/wrtipr.o 		: $(PA)/wrtipr.f                                       \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com  \
                        $(INC)/procan.com

$(PA)/wrtiprhdr.o 	: $(PA)/wrtiprhdr.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/chmstry.com     \
                        $(INC)/filunit.com $(INC)/grid.com $(INC)/procan.com

$(PA)/wrtirr.o 		: $(PA)/wrtirr.f                                       \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/filunit.com  \
                        $(INC)/grid.com $(INC)/tracer.com $(INC)/procan.com

$(PA)/wrtirrhdr.o 	: $(PA)/wrtirrhdr.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/chmstry.com     \
                        $(INC)/filunit.com $(INC)/grid.com $(INC)/procan.com

$(PIG)/gresdriv.o 	: $(PIG)/gresdriv.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/pigsty.com $(INC)/chmstry.com $(INC)/ahomap.com \
                        $(INC)/procan.com $(INC)/tracer.com

$(PIG)/gresgrow.o 	: $(PIG)/gresgrow.f                                    \
                        $(INC)/camx.prm $(INC)/pigsty.com $(INC)/filunit.com

$(PIG)/gresmscl.o 	: $(PIG)/gresmscl.f                                    \
                        $(INC)/camx.prm $(INC)/pigsty.com

$(PIG)/piginit.o 	: $(PIG)/piginit.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/ptemiss.com     \
                        $(INC)/pigsty.com $(INC)/chmstry.com $(INC)/flags.com  \
                        $(INC)/filunit.com

$(PIG)/pigprep.o 	: $(PIG)/pigprep.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/flags.com       \
                        $(INC)/pigsty.com $(INC)/filunit.com                   \
                        $(INC)/camxfld.com $(INC)/ptemiss.com                  \
                        $(INC)/chmstry.com

$(PIG)/pigwalk.o 	: $(PIG)/pigwalk.f                                     \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/filunit.com     \
                        $(INC)/camxfld.com $(INC)/pigsty.com

$(PIG)/walk1pig.o 	: $(PIG)/walk1pig.f                                    \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/pigsty.com      \
                        $(INC)/bndary.com $(INC)/filunit.com

$(PIG)/wrtpig.o 	: $(PIG)/wrtpig.f                                      \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/pigsty.com   \
                        $(INC)/flags.com $(INC)/chmstry.com

$(RTRAC)/chemrt.o 	: $(RTRAC)/chemrt.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/chmstry.com     \
                        $(INC)/grid.com $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/cvticrt.o 	: $(RTRAC)/cvticrt.f                                   \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/camx.com      \
                        $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/drydeprt.o 	: $(RTRAC)/drydeprt.f                                  \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/bndary.com      \
                        $(INC)/deposit.com $(INC)/chmstry.com                  \
                        $(INC)/filunit.com $(INC)/tracer.com                   \
                        $(INC)/rtracchm.com

$(RTRAC)/empreprt.o 	: $(RTRAC)/empreprt.f                                  \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/grid.com     \
                        $(INC)/flags.com $(INC)/chmstry.com $(INC)/tracer.com  \
                        $(INC)/rtracchm.com

$(RTRAC)/hdrcprt.o 	: $(RTRAC)/hdrcprt.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/tracer.com      \
                        $(INC)/rtracchm.com $(INC)/filunit.com

$(RTRAC)/rdarrt.o 	: $(RTRAC)/rdarrt.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/chmstry.com $(INC)/grid.com    \
                        $(INC)/bndary.com $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/rdbcrt.o 	: $(RTRAC)/rdbcrt.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/chmstry.com $(INC)/grid.com $(INC)/flags.com    \
                        $(INC)/bndary.com $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/rdchmrt.o 	: $(RTRAC)/rdchmrt.f                                   \
                        $(INC)/camx.prm $(INC)/filunit.com $(INC)/grid.com     \
                        $(INC)/chmstry.com $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/rdicrt.o 	: $(RTRAC)/rdicrt.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/chmstry.com $(INC)/grid.com $(INC)/flags.com    \
                        $(INC)/bndary.com $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/rdptrt.o 	: $(RTRAC)/rdptrt.f                                    \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/filunit.com     \
                        $(INC)/flags.com $(INC)/chmstry.com $(INC)/grid.com    \
                        $(INC)/bndary.com $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/rdrcprt.o 	: $(RTRAC)/rdrcprt.f                                   \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/tracer.com      \
                        $(INC)/rtracchm.com $(INC)/procan.com $(INC)/filunit.com

$(RTRAC)/startrt.o 	: $(RTRAC)/startrt.f                                   \
                        $(INC)/camx.prm $(INC)/grid.com $(INC)/flags.com       \
                        $(INC)/tracer.com

$(RTRAC)/wetdeprt.o 	: $(RTRAC)/wetdeprt.f                                  \
                        $(INC)/camx.prm $(INC)/bndary.com $(INC)/chmstry.com   \
                        $(INC)/tracer.com $(INC)/rtracchm.com

$(RTRAC)/wrrcprt.o 	: $(RTRAC)/wrrcprt.f                                   \
                        $(INC)/camx.prm $(INC)/camx.com $(INC)/tracer.com      \
                        $(INC)/rtracchm.com $(INC)/filunit.com

$(SOAP)/soapdat.o 	: $(SOAP)/soapdat.f                                    \
                        $(INC)/soap.com

$(SOAP)/soap.o 		: $(SOAP)/soap.f                                       \
                        $(INC)/soap.com

$(AER)/addit.o		: $(AER)/addit.f

$(AER)/aerchem.o	: $(AER)/aerchem.f                                     \
                        $(INC)/dynamic.inc $(INC)/dbg.inc

$(AER)/aero_err.o	: $(AER)/aero_err.f                                    \
                        $(INC)/camx.prm $(INC)/chmstry.com $(INC)/camxfld.com  \
                        $(INC)/filunit.com $(INC)/grid.com

$(AER)/aqchem.o		: $(AER)/aqchem.f                                      \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/dropcom.inc

$(AER)/aqdist.o		: $(AER)/aqdist.f

$(AER)/aqfex1.o		: $(AER)/aqfex1.f                                      \
                        $(INC)/aerpar.inc $(INC)/droppar.inc

$(AER)/aqfex2.o		: $(AER)/aqfex2.f                                      \
                        $(INC)/aerpar.inc $(INC)/droppar.inc

$(AER)/aqintegr1.o	: $(AER)/aqintegr1.f                                   \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/math.inc

$(AER)/aqintegr2.o	: $(AER)/aqintegr2.f                                   \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/math.inc

$(AER)/aqjex.o		: $(AER)/aqjex.f

$(AER)/aqoperator1.o	: $(AER)/aqoperator1.f                                 \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/dropcom.inc

$(AER)/aqoperator2.o	: $(AER)/aqoperator2.f                                 \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/dropcom.inc

$(AER)/aqrates1.o	: $(AER)/aqrates1.f                                    \
                        $(INC)/aerpar.inc $(INC)/droppar.inc                   \
                        $(INC)/dropcom.inc $(INC)/math.inc

$(AER)/aqrates2.o	: $(AER)/aqrates2.f                                    \
                        $(INC)/aerpar.inc $(INC)/droppar.inc                   \
                        $(INC)/dropcom.inc $(INC)/math.inc

$(AER)/baddit.o		: $(AER)/baddit.f

$(AER)/block.o		: $(AER)/block.f                                       \
                        $(INC)/dynamic.inc

$(AER)/bmass.o		: $(AER)/bmass.f

$(AER)/coagul.o		: $(AER)/coagul.f                                      \
                        $(INC)/dynamic.inc

$(AER)/constants.o	: $(AER)/constants.f                                   \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/dropcom.inc

$(AER)/decisions.o	: $(AER)/decisions.f                                   \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/dropcom.inc

$(AER)/diameter.o	: $(AER)/diameter.f                                    \
                        $(INC)/dynamic.inc

$(AER)/diff.o		: $(AER)/diff.f                                        \
                        $(INC)/dynamic.inc

$(AER)/differ.o		: $(AER)/differ.f

$(AER)/diffund.o	: $(AER)/diffund.f                                     \
                        $(INC)/dynamic.inc $(INC)/equaer.inc

$(AER)/dropinit.o	: $(AER)/dropinit.f                                    \
                        $(INC)/aerpar.inc $(INC)/droppar.inc                   \
                        $(INC)/dropcom.inc $(INC)/section.inc                  \
			$(INC)/camx_aero.inc

$(AER)/electro.o	: $(AER)/electro.f

$(AER)/eqpart.o		: $(AER)/eqpart.f                                      \
                        $(INC)/dynamic.inc $(INC)/equaer.inc

$(AER)/eqparto.o	: $(AER)/eqparto.f                                     \
                        $(INC)/dynamic.inc $(INC)/soap.com

$(AER)/equaer.o		: $(AER)/equaer.f                                      \
                        $(INC)/dynamic.inc $(INC)/equaer.inc

$(AER)/fullequil.o	: $(AER)/fullequil.f

$(AER)/hybrid.o		: $(AER)/hybrid.f                                      \
                        $(INC)/dynamic.inc

$(AER)/linint.o		: $(AER)/linint.f

$(AER)/madm.o		: $(AER)/madm.f                                        \
                        $(INC)/dynamic.inc

$(AER)/mass.o		: $(AER)/mass.f

$(AER)/negchk.o		: $(AER)/negchk.f                                      \
                        $(INC)/dynamic.inc

$(AER)/newdist.o	: $(AER)/newdist.f                                     \
                        $(INC)/dynamic.inc

$(AER)/nucl.o		: $(AER)/nucl.f                                        \
                        $(INC)/dynamic.inc

$(AER)/qsaturation.o	: $(AER)/qsaturation.f

$(AER)/react.o		: $(AER)/react.f                                       \
                        $(INC)/aerpar.inc $(INC)/droppar.inc

$(AER)/state.o		: $(AER)/state.f

$(AER)/steady.o		: $(AER)/steady.f

$(AER)/step.o		: $(AER)/step.f                                        \
                        $(INC)/dynamic.inc

$(AER)/svode.o		: $(AER)/svode.f

$(AER)/values.o		: $(AER)/values.f

$(AER)/vsrm.o		: $(AER)/vsrm.f                                        \
                        $(INC)/aerpar.inc $(INC)/droppar.inc $(INC)/dropcom.inc

$(AER)/slsode.o		: $(AER)/slsode.f

