help([[
Load environment to build UPP on ursa
]])

setenv("LMOD_TMOD_FIND_FIRST","yes")
prepend_path("MODULEPATH", "/lustre/desc1/scratch/epicufsrt/contrib/modulefiles_extra")
prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")
-- prepend_path("MODULEPATH", "/glade/work/epicufsrt/contrib/spack-stack/derecho/spack-stack-1.9.2/envs/ue-oneapi-2024.2.1/install/modulefiles/cray-mpich/8.1.29-4natrhl/gcc/12.2.0")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
load(pathJoin("stack-oneapi", stack_oneapi_ver))

-- stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"
-- load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

stack_impi_ver=os.getenv("stack_cray_mpich_ver") or "8.1.29"
load(pathJoin("stack-cray-mpich", stack_cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

load("upp_common")

setenv("CC", "mpicc")
setenv("CXX", "mpic++")
setenv("FC", "mpifort")
setenv("F90", "mpifort")
--setenv("I_MPI_CC", "icx")
--setenv("I_MPI_CXX", "icpx")
--setenv("I_MPI_FC", "ifort")
--setenv("I_MPI_F77", "ifort")
--setenv("I_MPI_F90", "ifort")

whatis("Description: UPP build environment")
