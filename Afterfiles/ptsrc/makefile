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

FC = pgf77

FLGS = -tp p6 -pc 64 -Kieee -Mdalign -Mextend -Mnoframe -byteswapio -Wl,-Bstatic -Minfo

after_SAPRC: after_SAPRC.f
	$(FC) $(FLGS) -c after_SAPRC.f
	$(FC) $(FLGS) -o after_SAPRC after_SAPRC.o
	rm after_SAPRC.o

after_SAPRC_bins: after_SAPRC_bins.f
	$(FC) $(FLGS) -c after_SAPRC_bins.f
	$(FC) $(FLGS) -o after_SAPRC_bins after_SAPRC_bins.o
	rm after_SAPRC_bins.o

after_SAPRC_neg5: after_SAPRC_neg5.f
	$(FC) $(FLGS) -c after_SAPRC_neg5.f
	$(FC) $(FLGS) -o after_SAPRC_neg5 after_SAPRC_neg5.o
	rm after_SAPRC_neg5.o

after_SAPRC_neg5_bins: after_SAPRC_neg5_bins.f
	$(FC) $(FLGS) -c after_SAPRC_neg5_bins.f
	$(FC) $(FLGS) -o after_SAPRC_neg5_bins after_SAPRC_neg5_bins.o
	rm after_SAPRC_neg5_bins.o

after_SAPRC_plus5: after_SAPRC_plus5.f
	$(FC) $(FLGS) -c after_SAPRC_plus5.f
	$(FC) $(FLGS) -o after_SAPRC_plus5 after_SAPRC_plus5.o
	rm after_SAPRC_plus5.o

after_SAPRC_plus5_bins: after_SAPRC_plus5_bins.f
	$(FC) $(FLGS) -c after_SAPRC_plus5_bins.f
	$(FC) $(FLGS) -o after_SAPRC_plus5_bins after_SAPRC_plus5_bins.o
	rm after_SAPRC_plus5_bins.o

after_SAPRC_plus25: after_SAPRC_plus25.f
	$(FC) $(FLGS) -c after_SAPRC_plus25.f
	$(FC) $(FLGS) -o after_SAPRC_plus25 after_SAPRC_plus25.o
	rm after_SAPRC_plus25.o

after_SAPRC_plus25_bins: after_SAPRC_plus25_bins.f
	$(FC) $(FLGS) -c after_SAPRC_plus25_bins.f
	$(FC) $(FLGS) -o after_SAPRC_plus25_bins after_SAPRC_plus25_bins.o
	rm after_SAPRC_plus25_bins.o

after_CAMx_AV: after_CAMx_AV.f
	$(FC) $(FLGS) -c after_CAMx_AV.f
	$(FC) $(FLGS) -o after_CAMx_AV after_CAMx_AV.o
	rm after_CAMx_AV.o

after_SAPRC_lay: after_SAPRC_lay.f
	$(FC) $(FLGS) -c after_SAPRC_lay.f
	$(FC) $(FLGS) -o after_SAPRC_lay after_SAPRC_lay.o
	rm after_SAPRC_lay.o

after_SAPRC_lay_plus5: after_SAPRC_lay_plus5.f
	$(FC) $(FLGS) -c after_SAPRC_lay_plus5.f
	$(FC) $(FLGS) -o after_SAPRC_lay_plus5 after_SAPRC_lay_plus5.o
	rm after_SAPRC_lay_plus5.o

after_SAPRC_height: after_SAPRC_height.f
	$(FC) $(FLGS) -c after_SAPRC_height.f
	$(FC) $(FLGS) -o after_SAPRC_height after_SAPRC_height.o
	rm after_SAPRC_height.o

ptsource: ptsource.f
	$(FC) $(FLGS) -c ptsource.f
	$(FC) $(FLGS) -c plumeris.f
	$(FC) $(FLGS) -o ptsource ptsource.o plumeris.o 
	rm plumeris.o
	rm ptsource.o

pntread: pntread.f
	$(FC) $(FLGS) -c pntread.f
	$(FC) $(FLGS) -o pntread pntread.o
	rm pntread.o

clean: after
	rm *.o after
