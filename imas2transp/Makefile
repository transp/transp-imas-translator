# Makefile to link imas2transp with unoptimized (-g -O0) build of IMAS library v. imas/3.17.0/ual/3.8.0
#
# prerequisites:
# 1) IMAS_INSTALL_DIR=/home/ITER/carlssj/local/imas (has unoptimized build of IMAS library)
# 2) a local TRANSP installation (at least the libraries and evironment variables used here)

F90 = gfortran
#COPTS = -cpp -g -O0 -fsanitize=address
COPTS = -cpp -g -O0
#COPTS = -cpp -O2
#INCLUDES = -I/home/ITER/carlssj/projects/imas-3.7.0/access-layer/fortraninterface/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
INCLUDES = -I$(IMAS_INSTALL_DIR)/include/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
#LIBS = -L$(IMAS_PREFIX)/lib -limas-gfortran -limas \

LIBS = \
	$(LOCAL)/lib/uflib.a \
	$(LOCAL)/lib/ufhdf.a \
	$(LOCAL)/lib/vaxonly.a \
	$(LOCAL)/lib/mds_sub.a \
	$(LOCAL)/lib/mdstransp.a \
	$(LOCAL)/lib/pspline.a \
	$(LOCAL)/lib/portlib.a \
	$(LOCAL)/lib/ezcdf.a \
	$(LOCAL)/lib/comput.a \
	-L$(IMAS_INSTALL_DIR)/lib -l:libimas-gfortran.so.3.17.0 -limas \
	$(MDSPLUS_DIR)/lib/libMdsObjectsCppShr.a $(MDSPLUS_DIR)/lib/libMdsLib.a \
	$(MDSPLUS_DIR)/lib/libMdsIpShr.a $(MDSPLUS_DIR)/lib/libTdiShr.a $(MDSPLUS_DIR)/lib/libTreeShr.a \
	$(MDSPLUS_DIR)/lib/libMdsShr.a $(MDSPLUS_DIR)/lib/libMdsIpShr.a \
	$(NETCDFHOME)/lib/libnetcdff.a $(NETCDFHOME)/lib/libnetcdf.a $(EBROOTOPENBLAS)/lib/libopenblas.a \
	-lxml2 -lrt -lpthread -L$(EBROOTGCC)/lib64 -lstdc++ -ldl

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