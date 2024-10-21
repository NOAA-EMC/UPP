help([[
Load environment to build UPP on orion
]])


prepend_path("MODULEPATH", "/work/noaa/epic/role-epic/spack-stack/orion/spack-stack-1.8.0/envs/ue-intel-2021.9.0/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2021.9.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.9.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

python_ver=os.getenv("python_ver") or "3.10.8"
load(pathJoin("python", python_ver))

load("upp_common")

setenv("CC","mpiicc")
setenv("CXX","mpiicpc")
setenv("FC","mpiifort")

whatis("Description: UPP build environment")
