F90 = gfortran
COPTS = -cpp -O2
INCLUDES = -I$(IMAS_PREFIX)/include/gfortran -I$(NETCDFHOME)/include -I$(LOCAL)/mod

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
	/work/imas/core/installer/src/3.12.0/ual/3.6.3/fortraninterface/libimas-gfortran.a \
	/work/imas/core/installer/src/3.12.0/ual/3.6.3/lowlevel/libimas.a \
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
	rm -f imas2transp put_nubeam_data_to_ufiles *.o *.mod *~ A32???.???