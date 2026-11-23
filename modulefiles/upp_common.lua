whatis("Description: UPP build environment common libraries")

help([[Load UPP common libraries]])

local ufs_modules = {
  {["jasper"]          = "4.2.4" },
  {["zlib"]            = "1.2.13"  },
  {["libpng"]          = "1.6.37" },
  {["hdf5"]            = "1.14.5" },
  {["netcdf-c"]        = "4.9.2"  },
  {["netcdf-fortran"]  = "4.6.1"  },
  {["bacio"]           = "2.6.0"  },
  {["crtm"]            = "3.1.3"},
  {["g2"]              = "3.5.1"  },
  {["g2tmpl"]          = "1.17.0" },
  {["ip"]              = "5.4.0"  },
  {["w3emc"]           = "2.13.0" },
  {["nemsio"]          = "2.5.5"  },
  {["sigio"]           = "2.3.3"  },
  {["wrf-io"]          = "1.3.0"  },
  {["wgrib2"]          = "3.8.0"  },
  {["prod_util"]       = "2.1.2"  },
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end
