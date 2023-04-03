Pi2 on a compute cluster
------------------------

Building pi2 to run on a compute cluster is largely the same process than building for a single workstation.
Please see build instructions for the operating system the cluster is running (typically [Linux](build_instructions_linux.md)).

This file contains further notes about compiling and running pi2 in various cluster environments.

**Note**
Typically the modules mentioned below must be loaded when compiling AND when running the compiled code!



PSI Ra cluster
--------------

The required module load commands are:
```
module load gcc/9.3.0
module load Python/3.6.3
```

Please do not load Python module when compiling as there seems to be
some incompatibilities between some gcc and Python modules.

Additionally, there is a Ra-specific Makefile that you should activate by command
```
ln Makefile.local.Ra Makefile.local
```

Compile using command
```
make -j16 NO_OPENCL=1
```

The file `slurm_config_psi_ra.txt` can be used as a default SLURM configuration file,
i.e., replace `slurm_config.txt` by `slurm_config_psi_ra.txt`. You might need to edit it, though.
Typically in Ra cluster, good performance is achieved when
```
use_nn5 = False
```



CSC Puhti cluster
-----------------

The required module load command is:
```
module load fftw/3.3.10-mpi-omp
```

Compile using command
```
make -j16 NO_OPENCL=1
```

The file slurm_config_csc_puhti.txt can be used as a template SLURM configuration file,
i.e., replace `slurm_config.txt` by `slurm_config_csc_puhti.txt`. You will need to edit it
to enter your project number.

Typically in Puhti cluster, good performance is achieved when NN5 is enabled with
```
use_nn5 = True.
```
Note that setting `use_nn5` to `False` might make the _entire cluster_ unusable and the
admins will get angry at you. This has something to do with interplay between pi2 and Lustre
filesystem. NN5 chunk size setting
```
chunk_size = [1024, 1024, 1024]
```
has worked fine, but other sizes might be faster.

CSC prefers low number of long jobs, so I suggest to use
```
max_parallel_submit_count = 10
promote_threshold = 3
```
In Puhti, jobs seem to fail into weird errors if the entire node is not reserved for the job.
Therefore, it is suggested to set `sbatch` arguments such that
```
--tasks 1 --cpus-per-task=40
```
and job init commands to
```
job_init_commands = export OMP_NUM_THREADS=40
```



ESRF clusters
-------------

Some of the required libraries are available from the modules system, e.g.,
```
module load fftw
```
The rest you will need to install locally or compile yourself.
Do not forget to add the locally installed libraries in the PATH.

Additionally, there is an ESRF-specific Makefile that you should activate by command
```
ln Makefile.local.ESRF Makefile.local
```

Compile pi2 using command
```
make -j16 NO_OPENCL=1
```

The file slurm_config_esrf.txt can be used as a starting point for SLURM
settings file. You will need to modify it, though! Please see notes in the
template file.

