-- ---------------------------------------------------------------------------
-- David Huber 06/2021, Set up config. with the hpc-stack NCEPLIBS.
-- Innocent Souopgui 11/2023, Update to use spack-stack
-- David Huber 1/24, Update to use spack-stack v1.6.0
-- ---------------------------------------------------------------------------


prepend_path("MODULEPATH", "/data/prod/jedi/spack-stack/spack-stack-1.6.0/envs/upp-addon-env/install/modulefiles/Core")

local stack_intel_ver=os.getenv("stack_intel_ver") or "2021.5.0"
local stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.0"
local cmake_ver=os.getenv("cmake_ver") or "3.23.1"

load(pathJoin("stack-intel", stack_intel_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
load(pathJoin("cmake", cmake_ver))

load("upp_common")

setenv("CC","mpiicc")
setenv("CXX","mpiicpc")
setenv("FC","mpiifort")

whatis("Description: UPP build environment")

