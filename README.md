# transp-imas-translator

Source code for transp-to-imas and imas-to-transp data translators.
TRANSP and IMAS must be built before these translators.
TRANSP and IMAS environment variables must also be properly defined,
including UAL, ids_path, treename and MDSPLUS_TREE_BASE_0.
The Makefiles should be modified for your system.

The input for transp2imas is a TRANSP name list and NetCDF output file,
e.g. 38601I02TR.DAT and 38601I02.CDF, which are distributed with the
source code.

Current lead developer: Johan Carlsson, PPPL

Original developer: Jin Chen, PPPL
