#  Makefile for After post-processing programs
#
#  Modified by Ben Murphy - 3/11/08
#
#################################################################
#
#  Portland Group Compiler Flags (pgf77)
#
#    -tp p6		specify the target processor 
#			  (p6=pentium pro/pentium II)
#    -pc 64		precision option (64=double precision)
#    -Kieee		Perform floating-point operations 
#			  in strict conformance with the IEEE 
#			  754 standard
#    -Mdalign		Align doubles in structures on 8-byte
#			  boundaries
#    -Mextend		Allow 132 column source lines
#    -Mnoframe		Do not set up a true stack frame pointer 
#			  for fxns
#    -byteswapio	Swap from big-endian to little-endian
#			  or vice-versa on input/output of
#			  unformatted FORTRAN data
#    -W1,		Pass each comma-delimeted option to 
#			  the linker
#       -Bstatic	Static linking
#    -Minfo		Emit useful info to stderr
#
#    -c			skip the link step; compile and assemble 
#			  only
#    -o			output file
#
#  Debugging Options
#
#    -g			generate symbolic debug info.
#			sets optimizations level to zero.
#
################################################################


FC = /usr/pgi/linux86-64/6.0/bin/pgf77

FLGS =  -pc 64 -Kieee -Mdalign -mcmodel medium -Mextend -Mnoframe -byteswapio -Wl,-Bstatic -Minfo -Minform=warn
#FLGS = -fast

AMB = cmp_improve.f cmp_airs.f cmp_paqs.f sub_after.f sub_after_di.f
AFT = after_saprc.f after_di_saprc.f after_di_saprc_vbins.f
AFTL = after_saprc_lay.f after_di_saprc.f after_di_saprc_vbins.f

after_saprc: after_saprc.f
	$(FC) $(FLGS) -c -g after_saprc.f
	$(FC) $(FLGS) -g after_saprc.o -o after_saprc

after_di_saprc: after_di_saprc.f
	$(FC) $(FLGS) -c -g after_di_saprc.f
	$(FC) $(FLGS) -g after_di_saprc.o -o after_di_saprc

after_di: after_di.f
	$(FC) $(FLGS) -c -g after_di.f
	$(FC) $(FLGS) -g after_di.o -o after_di

cmp_improve: cmp_improve.f sub_after.f
	$(FC) $(FLGS) -c sub_after.f cmp_improve.f
	$(FC) $(FLGS) sub_after.o cmp_improve.o -o cmp_improve

cmp_airs: cmp_airs.f sub_after.f
	$(FC) $(FLGS) -c sub_after.f cmp_airs.f
	$(FC) $(FLGS) sub_after.o cmp_airs.o -o cmp_airs

cmp_paqs: cmp_paqs.f sub_after_di.f
	$(FC) $(FLGS) -c sub_after_di.f cmp_paqs.f
	$(FC) $(FLGS) sub_after_di.o cmp_paqs.o -o cmp_paqs

afterwrap: $(AMB) $(AFT) afterwrap.f 
	$(FC) $(FLGS) -c $(AMB) $(AFT) afterwrap.f
	$(FC) $(FLGS) $(AMB:.f=.o) $(AFT:.f=.o) afterwrap.o -o afterwrap

afterwrap.lay: $(AMB) $(AFTL) afterwrap.f 
	$(FC) $(FLGS) -c $(AMB) $(AFTL) afterwrap.f
	$(FC) $(FLGS) $(AMB:.f=.o) $(AFTL:.f=.o) afterwrap.o -o afterwrap.lay

afterwrap_christos: after_saprc_2.f after_di_saprc_2.f sub_after_2.f afterwrap_christos.f after_dimap_saprc_2.f
	$(FC) $(FLGS) -c after_saprc_2.f after_di_saprc_2.f sub_after_2.f afterwrap_christos.f after_dimap_saprc_2.f
	$(FC) $(FLGS) after_saprc_2.o after_di_saprc_2.o sub_after_2.o afterwrap_christos.o after_dimap_saprc_2.o -o afterwrap_christos
	rm *.o

afterwrap_2: after_saprc_2.f after_di_saprc_2.f sub_after_2.f afterwrap_2.f after_dimap_saprc_2.f
	$(FC) $(FLGS) -c after_saprc_2.f after_di_saprc_2.f sub_after_2.f afterwrap_2.f after_dimap_saprc_2.f
	$(FC) $(FLGS) after_saprc_2.o after_di_saprc_2.o sub_after_2.o afterwrap_2.o after_dimap_saprc_2.o -o afterwrap_2
	rm *.o

afterwrap_CB4: after_saprc_2.f after_di_saprc_2.f sub_after_CB4.f afterwrap_2.f after_dimap_saprc_2.f
	$(FC) $(FLGS) -c after_saprc_2.f after_di_saprc_2.f sub_after_CB4.f afterwrap_2.f after_dimap_saprc_2.f
	$(FC) $(FLGS) after_saprc_2.o after_di_saprc_2.o sub_after_CB4.o afterwrap_2.o after_dimap_saprc_2.o -o afterwrap_CB4

avgrads: avgrads.f
	$(FC) $(FLGS) -v -c avgrads.f
	$(FC) $(FLGS) avgrads.o -o avgrads

avgOH: avgOH.f
	$(FC) $(FLGS) -c avgOH.f
	$(FC) $(FLGS) avgOH.o -o avgOH

clean:
	rm *.o avgrads
