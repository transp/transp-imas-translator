# Makefile to link imas2transp with unoptimized (-g -O0) build of IMAS library v. imas/3.17.0/ual/3.8.0
#
# prerequisites:
# 1) IMAS_INSTALL_DIR=/home/ITER/carlssj/local/imas (has unoptimized build of IMAS library)
# 2) a local TRANSP installation (at least the libraries and evironment variables used here)

F90 = gfortran
#COPTS = -cpp -g -O0 -fsanitize=address
COPTS = -cpp -g -O0 -fdefault-real-8
#COPTS = -cpp -O2
#INCLUDES = -I/home/ITER/carlssj/projects/imas-3.7.0/access-layer/fortraninterface/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
#INCLUDES = -I$(IMAS_INSTALL_DIR)/include/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
INCLUDES = -I$(IMAS_PREFIX)/include/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod -I$(PSPLINE_HOME)/include
#LIBS = -L$(IMAS_PREFIX)/lib -limas-gfortran -limas \

LIBS = \
	$(UFILES_ROOT)/lib/libuflib.a \
	$(UFILES_ROOT)/lib/libufhdf.a \
	$(UFILES_ROOT)/lib/libufport.a \
	$(LOCAL)/lib/libvaxonly.a \
	$(LOCAL)/lib/libtokyr.a \
	$(UFILES_ROOT)/lib/libmds_sub.a \
	$(UFILES_ROOT)/lib/libmdstransp.a \
	$(PSPLINE_HOME)/lib/libpspline.a \
	$(LOCAL)/lib/libportlib.a \
	$(LOCAL)/lib/libezcdf.a \
	$(LOCAL)/lib/libcomput.a \
	$(L_IMAS) \
	$(MDSPLUS_DIR)/lib/libMdsObjectsCppShr.a $(MDSPLUS_DIR)/lib/libMdsLib.a \
	$(MDSPLUS_DIR)/lib/libTdiShr.a $(MDSPLUS_DIR)/lib/libTreeShr.a \
	$(MDSPLUS_DIR)/lib/libMdsShr.a $(MDSPLUS_DIR)/lib/libMdsIpShr.so \
	$(NETCDF_FORTRAN_HOME)/lib/libnetcdff.a $(NETCDF_C_HOME)/lib/libnetcdf.so $(BLAS_LIB) \
	-lxml2 -lrt -lpthread -lstdc++ -ldl
#	-L$(IMAS_INSTALL_DIR)/lib -l:libimas-gfortran.so.3.17.0 -limas \

imas2transp : put_data_to_ufiles_module.o put_data_to_ufiles.o imas2transp.o
	@echo "   * $^ -> $@"
	$(F90) $(COPTS) $^ $(LIBS) -o $@

put_data_to_ufiles_module.o : put_data_to_ufiles_module.f90
	$(F90) $(COPTS) -c $<

put_data_to_ufiles.o : put_data_to_ufiles.f90
	$(F90) $(COPTS) -I. -c $<

imas2transp.o : imas2transp.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

put_nubeam_data_to_ufiles : put_nubeam_data_to_ufiles.o
	$(F90) $(COPTS) $^ $(LIBS) -o $@

put_nubeam_data_to_ufiles.o : put_nubeam_data_to_ufiles.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

.PHONY : clean
clean :
	rm -f imas2transp put_nubeam_data_to_ufiles *.o *.mod *~ A?????.??? gindex.dat ids_*.geq
