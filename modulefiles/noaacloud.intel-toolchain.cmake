# intel-toolchain.cmake for noaacloud build environments

# Prefer absolute paths to avoid environment reliance
set(CMAKE_C_COMPILER /apps/oneapi/mpi/latest/bin/mpiicc CACHE FILEPATH "")
set(CMAKE_CXX_COMPILER /apps/oneapi/mpi/latest/bin/mpiicpc CACHE FILEPATH "")
set(CMAKE_Fortran_COMPILER /apps/oneapi/mpi/latest/bin/mpiifort CACHE FILEPATH "")

# Archiver and ranlib (must be compatible with libimf.so and Intel’s linker)
set(CMAKE_AR /apps/oneapi/compiler/latest/linux/bin-llvm/llvm-ar CACHE FILEPATH "")
set(CMAKE_RANLIB /apps/oneapi/compiler/latest/linux/bin-llvm/llvm-ranlib CACHE FILEPATH "")
set(CMAKE_LINKER /apps/oneapi/compiler/latest/linux/bin/intel64/xild CACHE FILEPATH "")
