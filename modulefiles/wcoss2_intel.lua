help([[
Load environment to build UPP on WCOSS2
]])

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.5.0"
intel_ver=os.getenv("intel_ver") or "19.1.3.304"
craype_ver=os.getenv("craype_ver") or "2.7.17"
cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.19"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))
load(pathJoin("intel", intel_ver))
load(pathJoin("craype", craype_ver))
load(pathJoin("cray-mpich", cray_mpich_ver))

cmake_ver=os.getenv("cmake_ver") or "3.20.2"
load(pathJoin("cmake", cmake_ver))

local upp_modules = {
  {["jasper"]          = "2.0.25" },
  {["zlib"]            = "1.2.11"  },
  {["libpng"]          = "1.6.37" },
  {["hdf5-D"]          = "1.14.0" },
  {["netcdf-D"]        = "4.9.2"  },
  {["pnetcdf-D"]       = "1.12.2"  },
  {["bacio"]           = "2.4.1"  },
  {["crtm"]            = "2.4.0.1"},
  {["g2"]              = "3.5.1"  },
  {["g2tmpl"]          = "1.17.0" },
  {["ip"]              = "4.0.0"  },
  {["sp"]              = "2.3.3"  },
  {["w3emc"]           = "2.12.0" },
  {["nemsio"]          = "2.5.2"  },
  {["sigio"]           = "2.3.2"  },
  {["wrf_io"]          = "1.2.0"  },
}

for i = 1, #upp_modules do
  for name, default_version in pairs(upp_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

whatis("Description: post build environment")
