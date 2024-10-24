help([[
  This module loads libraries required for building and running UPP
  on the NOAA RDHPC machine Gaea C5 using Intel-2023.2.0.
]])

whatis([===[Loads libraries needed for building the UPP on Gaea C5 ]===])

prepend_path("MODULEPATH", "/autofs/ncrc-svm1_proj/epic/spack-stack/spack-stack-1.6.0/envs/upp-addon-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2023.2.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.28"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("upp_common")

unload("darshan-runtime")
unload("cray-libsci")

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

setenv("CMAKE_Platform","gaeac5.intel")
