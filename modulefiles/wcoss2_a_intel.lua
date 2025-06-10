help([[
Load environment to build UPP on WCOSS2 Acorn
]])

prepend_path("MODULEPATH", "/lfs/h1/emc/nceplibs/noscrub/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core")

stack_intel_ver=os.getenv("stack_intel_ver") or "2022.0.2.262"
load(pathJoin("stack-intel", stack_intel_ver))

stack_cray_ver=os.getenv("stack_cray_ver") or "8.1.9"
load(pathJoin("stack-cray-mpich", stack_cray_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("upp_common")

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

whatis("Description: UPP build environment")
