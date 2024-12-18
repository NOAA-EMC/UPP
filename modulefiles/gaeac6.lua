help([[
Load environment to build UPP on Jet
]])


prepend_path("MODULEPATH", "/autofs/ncrc-svm1_proj/epic/spack-stack/spack-stack-1.6.0/envs/unified-env-c6/install/modulefiles/Core")

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.5.0"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))

stack_intel_ver=os.getenv("stack_intel_ver") or "2023.2.0"
load(pathJoin("stack-intel", stack_intel_ver))

stack_cray_mpich_ver=os.getenv("stack_cray_mpich_ver") or "8.1.29"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.23.1"
load(pathJoin("cmake", cmake_ver))

load("upp_common")

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

whatis("Description: UPP build environment")
