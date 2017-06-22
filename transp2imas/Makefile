F90 = gfortran
COPTS = -cpp -fPIC -O0 -g -fsanitize=address
#COPTS = -cpp -fPIC -O0 -g
INCLUDES = -I$(UAL)/fortraninterface/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
LIBS = $(UAL)/fortraninterface/libimas-gfortran.a $(UAL)/lowlevel/libimas.a $(MDSPLUS_DIR)/lib/libMdsObjectsCppShr.a \
$(MDSPLUS_DIR)/lib/libMdsLib.a $(MDSPLUS_DIR)/lib/libMdsIpShr.a $(MDSPLUS_DIR)/lib/libTdiShr.a $(MDSPLUS_DIR)/lib/libTreeShr.a $(MDSPLUS_DIR)/lib/libMdsShr.a \
-L$(MDSPLUS_DIR)/lib -lMdsIpShr $(NETCDFHOME)/lib/libnetcdff.a $(NETCDFHOME)/lib/libnetcdf.a \
-L$(HOME)/local/mkl-11.3.2.181/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64_lin -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
-lxml2 -lpthread -L/home/jcarlsso/local/$(GCC)/lib64 -lstdc++ -ldl

# SPECIFIC TO ITER HPC CLUSTER
ifneq (,$(findstring iter.org, $(shell hostname)))
  F90 = ifort -g
  TRANSP_DIR = /home/ITER/andrer/transpsvn/Linux_x86_64/local
  LOCAL=$(TRANSP_DIR)
  INCLUDES = `pkg-config imas-ifort netcdf netcdf-fortran --cflags`  -I$(LOCAL)/mod
  LIBS = `pkg-config imas-ifort imas-lowlevel netcdf netcdf-fortran --libs` 
  LIBS += $(MDSPLUS_DIR)/lib/libMdsObjectsCppShr.a $(MDSPLUS_DIR)/lib/libMdsLib.a $(MDSPLUS_DIR)/lib/libMdsIpShr.a \
	$(MDSPLUS_DIR)/lib/libTdiShr.a $(MDSPLUS_DIR)/lib/libTreeShr.a $(MDSPLUS_DIR)/lib/libMdsShr.a \
	-L$(MDSPLUS_DIR)/lib -lMdsIpShr -L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential \
	-lmkl_core -lpthread -limf -lm
endif

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

transp2imas : transp2imas_mod.o transp2imas.o t1mhdgeq.o
	@echo "   * $^ -> $@"
	$(F90) transp2imas_mod.o transp2imas.o t1mhdgeq.o $(COPTS) $(LIBS_transp) $(INCLUDES) $(LIBS) -I. -o $@

transp2imas_mod.o : transp2imas_mod.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

t1mhdgeq.o : t1mhdgeq.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

transp2imas.o : transp2imas.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

.PHONY : clean
clean :
	rm -f transp2imas eqdskin fort.* transp2imas.msgs *.o *.mod *~
