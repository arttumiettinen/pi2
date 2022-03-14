
This file contains notes about compiling and running pi2 in various cluster environments.

NOTE: Typically the modules mentioned below must be loaded when compiling AND
when running the compiled code!



PSI Ra cluster
--------------

The required module load commands are:
module load gcc/9.3.0
module load Python/3.6.3

Please do not load Python module when compiling as there seems to be
some incompatibilities between some gcc and Python modules.

Additionally, there is a Ra-specific Makefile that you should activate by command
ln Makefile.local.Ra Makefile.local

Compile using command
make -j16 NO_OPENCL=1

The file slurm_config_psi_ra.txt can be used as a default SLURM configuration file, i.e. replace slurm_config.txt by slurm_config_psi_ra.txt. You might need to edit it, though.



CSC Puhti cluster
-----------------

The required module load commands are:
module load gcc/9.1.0
module load fftw/3.3.8-omp

Compile using command
make -j16 NO_OPENCL=1

The file slurm_config_csc_puhti.txt can be used as a template SLURM configuration file, i.e. replace slurm_config.txt by slurm_config_csc_puhti.txt. You will need to edit it to enter your project number.



