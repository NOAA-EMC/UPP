help([[
Build environment for UPP in container 
]])

prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/intel-oneapi-mpi/2021.9.0/intel/2021.10.0")
prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/intel/2021.10.0")

local stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
local stack_impi_ver=os.getenv("stack_impi_ver") or "2021.9.0"
local cmake_ver=os.getenv("cmake_ver") or "3.23.1"

load("gnu")
load(pathJoin("stack-intel", stack_intel_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
unload("gnu")

load(pathJoin("cmake", cmake_ver))
load("upp_common")

whatis("Description: UPP environment in container with Intel Compilers")
