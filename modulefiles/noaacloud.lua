help([[
Load environment to build UPP on NOAA Cloud
]])

prepend_path("MODULEPATH", "/contrib/spack-stack-rocky8/spack-stack-1.8.0/envs/ue-intel-2021.10.0/install/modulefiles/Core")
prepend_path("MODULEPATH", "/apps/modules/modulefiles")
prepend_path("MODULEPATH", "/apps/oneapi/modulefiles")

load("gnu/9.2.0")

stack_intel_ver = os.getenv("stack_intel_ver") or "2021.10.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_impi_ver = os.getenv("stack_impi_ver") or "2021.10.0"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("upp_common")

whatis("Description: UPP build environment")
