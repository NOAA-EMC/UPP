help([[
Build environment for UPP in container 
]])

prepend_path("MODULEPATH", "/opt/spack-stack/spack-stack-1.6.0/envs/unified-env/install/modulefiles/Core")

local stack_intel_ver=os.getenv("stack_intel_ver") or "2021.10.0"
local stack_impi_ver=os.getenv("stack_impi_ver") or "2021.9.0"
local cmake_ver=os.getenv("cmake_ver") or "3.23.1"

load("gnu")
load(pathJoin("stack-intel", stack_intel_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))
unload("gnu")

load(pathJoin("cmake", cmake_ver))

local ufs_modules = {
  {["jasper"]          = "2.0.32" },
  {["libpng"]          = "1.6.37" },
  {["hdf5"]            = "1.14.0" },
  {["netcdf-c"]        = "4.9.2"  },
  {["netcdf-fortran"]  = "4.6.1"  },
  {["bacio"]           = "2.4.1"  },
  {["crtm"]            = "2.4.0.1"},
  {["g2"]              = "3.5.1"  },
  {["g2tmpl"]          = "1.13.0" },
  {["ip"]              = "4.3.0"  },
  {["w3emc"]           = "2.10.0" },
  {["nemsio"]          = "2.5.4"  },
  {["sigio"]           = "2.3.2"  },
  {["wrf-io"]          = "1.2.0"  },
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

whatis("Description: UPP environment in container with Intel Compilers")
