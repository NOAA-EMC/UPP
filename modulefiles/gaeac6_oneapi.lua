help([[
  Load environment to build UPP on Gaea6
]])

whatis([===[Loads libraries needed for building the UPP on Gaea6 ]===])

prepend_path("MODULEPATH", "/ncrc/proj/epic/spack-stack/c6/spack-stack-2.1.1/envs/ue-oneapi-2025.2.1/modules/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2025.2.1"
load(pathJoin("stack-intel-oneapi-compilers", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.32"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.31.8"
load(pathJoin("cmake", cmake_ver))

load("upp_common")


setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

setenv("CMAKE_Platform","gaeac6.oneapi")
