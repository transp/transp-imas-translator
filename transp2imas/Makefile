# Makefile to link transp2imas with module IMAS library v. 3.8.1
#
# prerequisites:
#  -- a local TRANSP installation that is based on the GCC compiler

F90 = gfortran
COPTS = -cpp -g -O0 -fPIE
#COPTS = -cpp -O2
INCLUDES = -I$(IMAS_PREFIX)/include/gfortran -I$(NETCDF_FORTRAN_HOME)/include -I$(LOCAL)/mod -I$(PSPLINE_HOME)/include
#LIBS = -L$(IMAS_PREFIX)/lib -limas-gfortran -limas \

LIBS = \
	$(LOCAL)/lib/libkey_access.a \
	$(LOCAL)/lib/libtrxplib.a \
	$(LOCAL)/lib/libtrdatbuf_lib.a \
	$(LOCAL)/lib/libtrread.a \
	$(LOCAL)/lib/libtr_getnl.a \
	$(LOCAL)/lib/librp_kernel.a \
	$(LOCAL)/lib/librplot_mod.a \
	$(LOCAL)/lib/librplot_io.a \
	$(UFILES_ROOT)/lib/libureadsub.a \
        $(LOCAL)/lib/libtokyr.a \
	$(LOCAL)/lib/libxdatmgr.a \
	$(UFILES_ROOT)/lib/libuflib.a \
	$(UFILES_ROOT)/lib/libmds_sub.a \
	$(UFILES_ROOT)/lib/libufhdf.a \
	$(UFILES_ROOT)/lib/libufport.a \
	$(LOCAL)/lib/libold_xplasma.a \
	$(LOCAL)/lib/libxplasma2.a \
	$(LOCAL)/lib/libgeqdsk_mds.a \
	$(UFILES_ROOT)/lib/libmdstransp.a \
	$(LOCAL)/lib/libnscrunch.a \
	$(LOCAL)/lib/libfluxav.a \
	$(LOCAL)/lib/libr8bloat.a \
	$(LOCAL)/lib/libmclib.a \
	$(UFILES_ROOT)/lib/libr4smlib.a \
	$(LOCAL)/lib/libsmlib.a \
	$(UFILES_ROOT)/lib/libinterp_sub.a \
	$(LOCAL)/lib/libcomput.a \
	$(LOCAL)/lib/libvaxonly.a \
	$(PSPLINE_HOME)/lib/libpspline.a \
	$(LOCAL)/lib/libezcdf.a \
	$(LOCAL)/lib/liblsode.a \
	$(LOCAL)/lib/liblsode_linpack.a \
	$(LOCAL)/lib/libportlib.a \
	-L$(IMAS_PREFIX)/lib -limas-gfortran -limas \
	-L$(MDSPLUS_DIR)/lib -lMdsLib -lMdsShr \
	$(NETCDF_FORTRAN_HOME)/lib/libnetcdff.a $(NETCDF_C_HOME)/lib/libnetcdf.so $(BLAS_LIBS) \
	-lxml2 -lrt -lpthread -lblas -L$(EBROOTGCC)/lib64 -lstdc++ -ldl

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
