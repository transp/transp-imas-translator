# Makefile to link transp2imas with unoptimized (-g -O0) build of IMAS library v. 3.7.0
#
# prerequisites:
# 1) "module add imas" (or "module add imas/3.15.1/ual/3.7.0" if the default changes)
# 2) a local TRANSP installation (at least the libraries and evironment variables used here)

F90 = gfortran
#COPTS = -cpp -g -O0 -fsanitize=address
COPTS = -cpp -g -O0
#COPTS = -cpp -O2
INCLUDES = -I/home/ITER/carlssj/projects/imas-3.7.0/access-layer/fortraninterface/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
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
	$(LOCAL)/lib/smlib.a \
	$(LOCAL)/lib/interp_sub.a \
	$(LOCAL)/lib/comput.a \
	$(LOCAL)/lib/vaxonly.a \
	$(LOCAL)/lib/pspline.a \
	$(LOCAL)/lib/ezcdf.a \
	$(LOCAL)/lib/lsode.a \
	$(LOCAL)/lib/lsode_linpack.a \
	$(LOCAL)/lib/portlib.a \
	/home/ITER/carlssj/projects/imas-3.7.0/access-layer/fortraninterface/libimas-gfortran.a \
	/home/ITER/carlssj/projects/imas-3.7.0/access-layer/lowlevel/libimas.a \
	$(MDSPLUS_DIR)/lib/libMdsObjectsCppShr.a $(MDSPLUS_DIR)/lib/libMdsLib.a \
	$(MDSPLUS_DIR)/lib/libMdsIpShr.a $(MDSPLUS_DIR)/lib/libTdiShr.a $(MDSPLUS_DIR)/lib/libTreeShr.a \
	$(MDSPLUS_DIR)/lib/libMdsShr.a $(MDSPLUS_DIR)/lib/libMdsIpShr.a \
	$(NETCDFHOME)/lib/libnetcdff.a $(NETCDFHOME)/lib/libnetcdf.a $(EBROOTOPENBLAS)/lib/libopenblas.a \
	-lxml2 -lrt -lpthread -L$(EBROOTGCC)/lib64 -lstdc++ -ldl

transp2imas : transp2imas_mod.o transp2imas.o t1mhdgeq.o
	@echo "   * $^ -> $@"
	$(F90) $(COPTS) $^ $(LIBS) -o $@

transp2imas_mod.o : transp2imas_mod.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

t1mhdgeq.o : t1mhdgeq.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

transp2imas.o : transp2imas.f90
	git rev-parse HEAD > hash.txt && $(F90) $(COPTS) -I. $(INCLUDES) -c $<

.PHONY : clean
clean :
	rm -f transp2imas eqdskin fort.* transp2imas.msgs *.o *.mod hash.txt *~