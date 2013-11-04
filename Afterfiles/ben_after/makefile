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

FLGS =  -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio -Wl,-Bstatic -Minfo -Minform=warn
#FLGS = -fast

AMB = cmp_improve.f cmp_airs.f cmp_paqs.f sub_after.f sub_after_di.f
AFT = after_saprc.f after_di_saprc.f after_di_saprc_vbins.f
AFTL = after_saprc_lay.f after_di_saprc.f after_di_saprc_vbins.f

after_SAPRC: after_daily_map.f after_timeseries.f after_sub_SAPRC.f after_main.f after_diurnal_map.f
	$(FC) $(FLGS) -c after_daily_map.f after_timeseries.f after_sub_SAPRC.f after_main.f after_diurnal_map.f
	$(FC) $(FLGS) after_daily_map.o after_timeseries.o after_sub_SAPRC.o after_main.o after_diurnal_map.o -o after_SAPRC
	rm *.o

avgrads: avgrads.f
	$(FC) $(FLGS) -v -c avgrads.f
	$(FC) $(FLGS) avgrads.o -o avgrads

avgOH: avgOH.f
	$(FC) $(FLGS) -c avgOH.f
	$(FC) $(FLGS) avgOH.o -o avgOH

clean:
	rm *.o avgrads
