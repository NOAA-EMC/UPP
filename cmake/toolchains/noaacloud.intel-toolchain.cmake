# intel-toolchain.cmake for noaacloud build environments

# Prefer absolute paths to avoid environment reliance
set(CMAKE_C_COMPILER /apps/oneapi/mpi/latest/bin/mpiicc CACHE FILEPATH "")
set(CMAKE_CXX_COMPILER /apps/oneapi/mpi/latest/bin/mpiicpc CACHE FILEPATH "")
set(CMAKE_Fortran_COMPILER /apps/oneapi/mpi/latest/bin/mpiifort CACHE FILEPATH "")

# Checking for llvm bins
if(EXISTS "/apps/oneapi/compiler/latest/linux/bin-llvm")
	set(LLVM_BIN "/apps/oneapi/compiler/latest/linux/bin-llvm")
	set(INTEL_COMPILER_BIN "/apps/oneapi/compiler/latest/linux/bin/intel64")
elseif(EXISTS "/apps/oneapi/compiler/latest/bin/compiler")
	set(LLVM_BIN "/apps/oneapi/compiler/latest/bin/compiler")
	set(INTEL_COMPILER_BIN "/apps/oneapi/compiler/latest/bin")
else()
	message(WARNING "LLVM bin path not found. Using system defaults for AR, RANLIB, and LINKER.")
endif()

# Archiver and ranlib (must be compatible with libimf.so and Intel’s linker)
if(DEFINED LLVM_BIN)
	set(CMAKE_AR "${LLVM_BIN}/llvm-ar" CACHE FILEPATH "")
	set(CMAKE_RANLIB "${LLVM_BIN}/llvm-ranlib" CACHE FILEPATH "")
	set(CMAKE_LINKER "${INTEL_COMPILER_BIN}/xild" CACHE FILEPATH "")
endif()
