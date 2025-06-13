help([[
Load environment to build UPP on ursa
]])

prepend_path("MODULEPATH", "/contrib/spack-stack/spack-stack-1.9.1/envs/ue-oneapi-2024.2.1/install/modulefiles/Core")
prepend_path("MODULEPATH", "/scratch4/NAGAPE/epic/role-epic/modulefiles")

stack_oneapi_ver=os.getenv("stack_oneapi_ver") or "2024.2.1"
load(pathJoin("stack-oneapi", stack_oneapi_ver))

stack_impi_ver=os.getenv("stack_impi_ver") or "2021.13"
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

cmake_ver=os.getenv("cmake_ver") or "3.27.9"
load(pathJoin("cmake", cmake_ver))

local ufs_modules = {
  {["jasper"]          = "2.0.32" },
  {["zlib"]            = "1.2.13" },
  {["libpng"]          = "1.6.37" },
  {["hdf5"]            = "1.14.3" },
  {["netcdf-c"]        = "4.9.2"  },
  {["netcdf-fortran"]  = "4.6.1"  },
  {["bacio"]           = "2.4.1"  },
  {["crtm"]            = "2.4.0.1"},
  {["g2"]              = "3.5.1"  },
  {["g2tmpl"]          = "1.13.0" },
  {["ip"]              = "5.1.0"  },
  {["w3emc"]           = "2.10.0" },
  {["nemsio"]          = "2.5.4"  },
  {["sigio"]           = "2.3.3"  },
  {["wrf-io"]          = "1.2.0"  },
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

setenv("CC", "mpiicx")
setenv("CXX", "mpiicpx")
setenv("FC", "mpiifort")
setenv("I_MPI_CC", "icx")
setenv("I_MPI_CXX", "icpx")
setenv("I_MPI_FC", "ifort")
setenv("I_MPI_F77", "ifort")
setenv("I_MPI_F90", "ifort")

whatis("Description: UFS build environment")
