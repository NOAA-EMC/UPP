macro (setDiscover)
  message("Setting paths for Discover")
  option(BUILD_CORELIBS "Build core libs from source" ON)
  option(FIND_HDF5 "Try to Find HDF5 libraries" OFF)
  option(FIND_HDF5_HL "Try to Find HDF5 libraries" ON)
  option(HDF5_USE_STATIC_LIBRARIES "Use only static libraries for HDF5" OFF)
  set(BUILD_NCDIAG_SERIAL FALSE)
  set(HOST_FLAG "-xHOST" CACHE INTERNAL "Host Flag")
  set(MKL_FLAG "-mkl"  CACHE INTERNAL "MKL Flag")
  set(GSI_Intel_Platform_FLAGS "-DPOUND_FOR_STRINGIFY -O3 -fp-model source -assume byterecl -convert big_endian -g -traceback -D_REAL8_ ${OpenMP_Fortran_FLAGS} ${MPI_Fortran_COMPILE_FLAGS}" CACHE INTERNAL "GSI Fortran Flags")
  set(ENKF_Platform_FLAGS "-O3 ${HOST_FLAG} -warn all -implicitnone -traceback -fp-model strict -convert big_endian -DGFS -D_REAL8_ ${MPI3FLAG} ${OpenMP_Fortran_FLAGS}" CACHE INTERNAL "ENKF Fortran Flags")
  set(host "Discover" CACHE INTERNAL "")
  
  if( ENV{BASEDIR} )
    set(BASEDIR $ENV{BASEDIR}/Linux CACHE INTERNAL "")
  endif()

  set(ENV{MPI_HOME} $ENV{MPI_ROOT} )

endmacro()

