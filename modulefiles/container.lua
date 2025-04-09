help([[
Load environment to build UPP on Container
]])

prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.6.0/envs/fms-2024.01/install/modulefiles/Core")
--prepend_path("MODULEPATH", "/apps/modules/modulefiles")
prepend_path("MODULEPATH", "/root/modulefiles")

local stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
local stack_impi_ver=os.getenv("stack_impi_ver") or "2021.9.0"

load("gnu")
load(pathJoin("stack-intel", stack_intel_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
unload("gnu")

load("tbb/latest")
load("compiler-rt/latest")
load("oclfpga/latest")
load("compiler/latest")
load("mpi/latest")

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("upp_common")

setenv("CC","mpiicc")
setenv("CXX","mpiicpc")
setenv("FC","mpiifort")

whatis("Description: UPP build environment")
