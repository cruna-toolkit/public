### Build your own generic libaries

In case, you miss OpenBLAS, FFTW3, or HDF5 the `makefile` provides additional targets to install these libraries from source. The versions of the included sources are (`OpenBLAS 0.3.14` `FFTW 3.3.8` `HDF5 1.12.0`). Alternatively, please consider downloading an update from the original source distributor. To improve the performance of a library, we suggest adapting the generic build options (for example `./makefiles/makefile_ext_hdf5`) to your specific system.

Load your compiler and MPI.
> A tested combination is (`GCC 9.3.*` `OpenMPI 3.1.5`).

`make clean_all`  # cleans all libraries inside `./lib`  
`make -j folders`  # creates directories  
`make -j openBLAS`  
`make -j fftw3`  
`make -j hdf5`

All build artefacts will be stored inside `./lib`
