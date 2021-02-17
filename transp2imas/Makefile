# Makefile to link transp2imas with module IMAS library v. 3.8.1
#
# prerequisites:
# 1) "module add imas" (or "module add imas/3.17.2/ual/3.8.1" if the default changes)
# 2) a local TRANSP installation (at least the libraries and evironment variables used here)

F90 = gfortran
COPTS = -cpp -g -O0 -fPIE
#COPTS = -cpp -O2
INCLUDES = -I$(IMAS_PREFIX)/include/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod -I$(PSPLINE_HOME)/include
#LIBS = -L$(IMAS_PREFIX)/lib -limas-gfortran -limas \

LIBS = \
	$(LOCAL)/lib/key_access.a \
	$(LOCAL)/lib/trxplib.a \
	$(LOCAL)/lib/trdatbuf_lib.a \
	$(LOCAL)/lib/trread.a \
	$(LOCAL)/lib/tr_getnl.a \
	$(LOCAL)/lib/rp_kernel.a \
	$(LOCAL)/lib/rplot_io.a \
	$(LOCAL)/lib/xdatmgr.a \
	$(LOCAL)/lib/uflib.a \
	$(LOCAL)/lib/mds_sub.a \
	$(LOCAL)/lib/ufhdf.a \
	$(LOCAL)/lib/old_xplasma.a \
	$(LOCAL)/lib/xplasma2.a \
	$(LOCAL)/lib/geqdsk_mds.a \
	$(LOCAL)/lib/mdstransp.a \
	$(LOCAL)/lib/nscrunch.a \
	$(LOCAL)/lib/fluxav.a \
	$(LOCAL)/lib/r8bloat.a \
	$(LOCAL)/lib/mclib.a \
	$(LOCAL)/lib/r4smlib.a \
	$(LOCAL)/lib/smlib.a \
	$(LOCAL)/lib/interp_sub.a \
	$(LOCAL)/lib/comput.a \
	$(LOCAL)/lib/vaxonly.a \
	$(PSPLINE_HOME)/lib/libpspline.a \
	$(LOCAL)/lib/ezcdf.a \
	$(LOCAL)/lib/lsode.a \
	$(LOCAL)/lib/lsode_linpack.a \
	$(LOCAL)/lib/portlib.a \
	-L$(IMAS_PREFIX)/lib -limas-gfortran -limas \
	-L$(MDSPLUS_DIR)/lib -lMdsLib -lMdsShr \
	$(NETCDFHOME)/lib/libnetcdff.a $(NETCDFHOME)/lib/libnetcdf.so $(BLAS_LIBS) \
	-lxml2 -lrt -lpthread -L$(EBROOTGCC)/lib64 -lstdc++ -ldl

transp2imas : transp2imas_mod.o transp2imas.o t1mhdgeq.o
	@echo "   * $^ -> $@"
	$(F90) $(COPTS) $^ $(LIBS) -o $@

transp2imas_mod.o : transp2imas_mod.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

t1mhdgeq.o : t1mhdgeq.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

transp2imas.o : transp2imas.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

.PHONY : clean
clean :
	rm -f transp2imas eqdskin fort.* transp2imas.msgs *.o *.mod *~
