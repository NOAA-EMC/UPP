help([[
  This module loads libraries required for building and running UPP
  on the NOAA RDHPC machine Gaea using Intel-2023.2.0.
]])

whatis([===[Loads libraries needed for building the UPP on Gaea ]===])

prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/c6/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
load(pathJoin("stack-oneapi", stack_oneapi_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.32"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("upp_common")
load("zlib/1.2.13")

unload("darshan-runtime")

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

setenv("CMAKE_Platform","gaea.intel")
