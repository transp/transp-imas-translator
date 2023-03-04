# Makefile to link transp2imas with module TRANSP/IMAS conversion tools
#

all : 
	cd transp2imas; make
	cd imas2transp; make

.PHONY : clean
clean : 
	cd transp2imas; make clean
	cd imas2transp; make clean
  
