# transp-imas-translator

Source code for transp-to-imas and imas-to-transp data translators.
TRANSP and IMAS must be built before these translators.
TRANSP and IMAS environment variables must also be properly defined,
including ids_path and MDSPLUS_TREE_BASE_0.
The Makefiles should be modified for your system.

The transp-to-imas translator was recently updated to generate a
summary IDS (in addition to core_profiles, equilibrium, et cetera).
The newest version of the code has only been tested on the ITER cluster
with version 3.7.0 of IMAS (access-layer / data-dictionary). More
details on how to build the code can be found in transp2imas/README.

The input for transp2imas is a TRANSP name list and NetCDF output file,
e.g. 38601I02TR.DAT and 38601I02.CDF, which are distributed with the
source code.

The input for imas2transp are IMAS files ids_<treeid>.tree, 
ids_<treeid>.datafile and ids_<treeid>.characteristics located in the 
directory $MDSPLUS_TREE_BASE_0 .

Current lead developer: Johan Carlsson, PPPL

Original developer: Jin Chen, PPPL
