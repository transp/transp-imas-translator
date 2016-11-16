F90 = gfortran
COPTS = -cpp -fPIC -g -O0
INCLUDES = -I$(UAL)/fortraninterface/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
LIBS = $(UAL)/fortraninterface/libimas-gfortran.a $(UAL)/lowlevel/libimas.a $(MDSPLUS_DIR)/lib/libMdsObjectsCppShr.a \
$(MDSPLUS_DIR)/lib/libMdsLib.a $(MDSPLUS_DIR)/lib/libMdsIpShr.a $(MDSPLUS_DIR)/lib/libTdiShr.a $(MDSPLUS_DIR)/lib/libTreeShr.a $(MDSPLUS_DIR)/lib/libMdsShr.a \
-L$(MDSPLUS_DIR)/lib -lMdsIpShr $(NETCDFHOME)/lib/libnetcdff.a $(NETCDFHOME)/lib/libnetcdf.a -lhdf5_hl -lhdf5 \
-L$(HOME)/local/mkl-11.3.2.181/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64_lin -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
-lxml2 -lpthread -L/home/jcarlsso/local/gcc-6.1.0/lib64 -lstdc++ -ldl

LIBS_transp = \
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
	$(LOCAL)/lib/portlib.a

OBJ = transp2imas_mod.o transp2imas.o t1mhdgeq.o

transp2imas: $(OBJ)
	@echo "   * $^ -> $@"
	$(F90) $(OBJ) $(COPTS) $(LIBS_transp) $(INCLUDES) $(LIBS) -I. -o $@

transp2imas_mod.o: transp2imas_mod.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

t1mhdgeq.o: t1mhdgeq.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

transp2imas.o: transp2imas.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<
	
clean:
	rm -f transp2imas eqdskin fort.* transp2ids.* *.o *.mod *~
