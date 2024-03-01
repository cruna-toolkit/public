# Installation
## Preparatory steps

### 1.) Configuration
First create your own config file in your main CRUNA directory:  
`cp cruna.conf.example cruna.conf`  
Inside `cruna.conf` you can adjust the `USER`, and the `CASE`. 
Exemplary, [cruna.conf.example](../../cruna.conf.example) holds the user `unit_test` and the case `acoustic_gauss_pulse`. A `USER` and a `CASE` require corresponding case files located in `./usr/$USER/$CASE/*`.  
A GNU and an Intel build is supported by setting `COMP = gnu` or `intel`.
Furthermore `cruna.conf` specifies a `$GENERIC_MAKEFILE` and a `$LOCAL_MAKEFILE` (see next section).

### 2.) Makefile setup
The `makefile` is located in the main directory. This makefile always includes `cruna.conf` and then (target dependent) several sub-makefiles which all can be found in `./makefiles/*`. To build CRUNA universal compiler flags (`$GENERIC_MAKEFILE`) and machine specific options (`$LOCAL_MAKEFILE`) are included always. 
In most cases the default flags of the included generic makefile (`./makefiles/generic`) do not need to be changed. In contrast, a local makefile (`./makefiles/local/your-local-makefile`) needs to be written. We recommend to copy an existing local makefile (e.g. `./makefiles/local/with_system_libs.example`) and adopt it to your system. This local makefile should state the names of the compiler wrappers and provide all the details how to link which system library. All libraries (apart of the compiler and MPI itself) listed in the following step need to be linked here.

## Provide existing system libraries (recommended)

Instead of reinventing the wheel we recommend to base CRUNA on existing libraries (in case you have access to a system with optimized builds):

### 3.) Necessary compiler and libraries for an Intel build
- load your compiler and Intel-MKL (providing \$MKLROOT)
- load Intel-MPI (providing mpiifort)
- load HDF5 (providing \$HDF5\_ROOT)
> Tested combinations are (`Intel 2021.2` `Intel-MPI 2021.7`) or (`Intel 18` `Intel-MPI 2018.5`) and `HDF5 1.10.5` or `HDF5 1.12.0`.

### 4.) Necessary compiler and libraries for a GNU build
- load your compiler 
- load OpenMPI (providing mpifort) 
- load OpenBLAS (providing \$OPENBLAS_DIR)
- load FFTW3 (providing \$FFTW3\_ROOT)
- load HDF5 (providing \$HDF5\_ROOT)
> Tested combinations are (`GCC 9.3.*` `OpenMPI 3.1.5` `OpenBLAS 0.3.7` `FFTW 3.3.8` `HDF5 1.10.6`).

Only in case your system does not natively provide all necessary libraries the `makefile` provides addtional targets (for `OpenBLAS 0.3.14`, `FFTW 3.3.8`, `HDF5 1.12.0`). To build these please refer to the [generic library build](build_own_generic_libraries.md).

## Install CRUNA

`make clean`  
`make -j cruna`  

