F90 = gfortran
COPTS = -cpp -fPIC -g -O0
INCLUDES = -I$(UAL)/fortraninterface/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod
LIBS = $(UAL)/fortraninterface/libimas-gfortran.a $(UAL)/lowlevel/libimas.a $(MDSPLUS_DIR)/lib/libMdsObjectsCppShr.a \
$(MDSPLUS_DIR)/lib/libMdsLib.a $(MDSPLUS_DIR)/lib/libMdsIpShr.a $(MDSPLUS_DIR)/lib/libTdiShr.a $(MDSPLUS_DIR)/lib/libTreeShr.a $(MDSPLUS_DIR)/lib/libMdsShr.a \
-L$(MDSPLUS_DIR)/lib -lMdsIpShr $(NETCDFHOME)/lib/libnetcdff.a $(NETCDFHOME)/lib/libnetcdf.a -lhdf5_hl -lhdf5 \
-L$(HOME)/local/mkl-11.3.2.181/compilers_and_libraries_2016.2.181/linux/mkl/lib/intel64_lin -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
-lxml2 -lpthread -L/home/jcarlsso/local/gcc-6.1.0/lib64 -lstdc++ -ldl

LIBS_ufile = \
	$(LOCAL)/lib/uflib.a \
	$(LOCAL)/lib/ufhdf.a \
	$(LOCAL)/lib/vaxonly.a \
	$(LOCAL)/lib/mds_sub.a \
	$(LOCAL)/lib/mdstransp.a
LIBS_ezspline = \
	$(LOCAL)/lib/pspline.a \
	$(LOCAL)/lib/portlib.a \
	$(LOCAL)/lib/ezcdf.a \
	$(LOCAL)/lib/comput.a

imas2transp : put_data_to_ufiles_module.o put_data_to_ufiles.o imas2transp.o
	$(F90) $(COPTS) -o $@ $^ $(LIBS_ufile) $(LIBS_ezspline) $(LIBS)

put_data_to_ufiles_module.o : put_data_to_ufiles_module.f90
	$(F90) $(COPTS) -c $<

put_data_to_ufiles.o : put_data_to_ufiles.f90
	$(F90) $(COPTS) -I. -c $<

imas2transp.o : imas2transp.f90
	$(F90) $(COPTS) -I. $(INCLUDES) -c $<

.PHONY : clean
clean :
	rm -f imas2transp *.o *.mod *~
