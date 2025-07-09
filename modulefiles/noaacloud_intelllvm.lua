help([[
Load environment to build UPP on NOAA Cloud with LLVM compilers
]])

prepend_path("MODULEPATH", "/contrib/spack-stack-rocky8/spack-stack-1.9.1/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")
prepend_path("MODULEPATH", "/apps/modules/modulefiles")

local gcc_ver=os.getenv("gcc_ver") or "13.2.0"
local stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
local stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"

load(pathJoin("gnu", gcc_ver))
load(pathJoin("stack-oneapi", stack_oneapi_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("upp_common")

setenv("I_MPI_CC", "/apps/oneapi/compiler/2024.2/bin/icx")
setenv("I_MPI_CXX", "/apps/oneapi/compiler/2024.2/bin/icpx")
setenv("I_MPI_F90", "/apps/oneapi/compiler/2024.2/bin/ifort")
setenv("CC","/apps/oneapi/mpi/latest/bin/mpiicx")
setenv("CXX","/apps/oneapi/mpi/latest/bin/mpiicpx")
setenv("FC","/apps/oneapi/mpi/latest/bin/mpiifx")

whatis("Description: UPP build environment")
