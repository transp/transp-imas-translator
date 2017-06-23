# --------------------------------------------------------------------------------
# Instructions to install, compile and execute the translator on ITER hpc cluster
# --------------------------------------------------------------------------------

1) Configure your environment:
module purge
module load imas/3.10.0/ual/3.5.4
module load intel/12.0.2
module load MDSplus/7.1.14-goolf-1.5.16-Java-1.7.0_79-Python-2.7.9
module load netcdf/4.3.2
module load netcdf-fortran/4.4.2-intel-12.0.2
imasdb test
export LOCAL=[transp-path] (e.g. $THOME/transpsvn/Linux_x86_64/local)

2 Compile the translator:
cd [your-path]/transp-imas-translator/transp2imas
make clean
make

3) Execution:
One needs to either copy the NetCDF and U-files to where the transp2imas executable is present or to go to 
where these files are located and execute the translator in this folder. In the present example, the folder, 
located in the official ITER scenario folder (temporarily replacing the future scenario database), contains 
the NetCDF file 79691H01.CDF and U-file (or namelist) 79691H01TR.DAT

-----> Example:
cd /work/imas/shared/scenarios/IOS/ITER_baseline/JET/transp_input/
[your-path]/transp-imas-translator/transp2imas/transp2imas 79691H01

-----> Typical logfile:
Open shot and write in IMAS !
transp runid: 79691H01
ids shot: 796
ids run: 91
[....]
transp2imas: save eq ids
transp2imas: save cp ids
transp2imas: save ct ids
transp2imas: save nbi ids
Close shot in IMAS!

4) What you get:
The IMAS converted file appears in the local database $HOME/public/imasdb/test/3/0 . It consists of three files: 
ids_7960091.tree, ids_7960091.datafile, ids_7960091.characteristics (shot 796, run 91) and contains the 
equilibrium, core_profiles core_transport and nbi IDSs.


