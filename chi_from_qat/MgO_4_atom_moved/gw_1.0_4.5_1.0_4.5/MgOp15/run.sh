gfortran -g -C -fbacktrace V2.F90 -L/opt/intel/composer_xe_2013.2.146/mkl/lib/intel64 -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -o cent_2.out
gfortran -g -C -fbacktrace cent_1.F90 -L/opt/intel/composer_xe_2013.2.146/mkl/lib/intel64 -lmkl_rt -lmkl_scalapack_lp64 -lmkl_blacs_openmpi_lp64 -o cent_1.out
./cent_1.out
./cent_2.out
echo
echo
echo
