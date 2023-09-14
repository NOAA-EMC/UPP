help([[
  This module loads libraries required for building and running UPP
  on the NOAA RDHPC machine Gaea using Intel-2022.0.2
]])

whatis([===[Loads libraries needed for building the UPP on Gaea ]===])

load_any(pathJoin("cmake", os.getenv("cmake_ver") or "3.20.1"),"cmake")

prepend_path("MODULEPATH","/lustre/f2/dev/role.epic/contrib/hpc-stack/intel-classic-2022.0.2/modulefiles/stack")
load(pathJoin("hpc", os.getenv("hpc_ver") or "1.2.0"))

load(pathJoin("intel-classic", os.getenv("intel_classic_ver") or "2022.0.2"))
load(pathJoin("cray-mpich", os.getenv("cray_mpich_ver") or "7.7.20"))
load(pathJoin("hpc-intel-classic", os.getenv("hpc_intel_classic_ver") or "2022.0.2"))
load(pathJoin("hpc-cray-mpich", os.getenv("hpc_cray_mpich_ver") or "7.7.20"))
load(pathJoin("libpng", os.getenv("libpng_ver") or "1.6.37"))

-- Needed at runtime:
load("alps")

local ufs_modules = {
  {["jasper"]      = "2.0.25"},
  {["zlib"]        = "1.2.11"},
  {["libpng"]      = "1.6.37"},
  {["hdf5"]        = "1.10.6"},
  {["netcdf"]      = "4.7.4"},
  {["pio"]         = "2.5.7"},
  {["esmf"]        = "8.3.0b09"},
  {["fms"]         = "2022.04"},
  {["bacio"]       = "2.4.1"},
  {["crtm"]        = "2.4.0"},
  {["g2"]          = "3.4.5"},
  {["g2tmpl"]      = "1.10.2"},
  {["ip"]          = "3.3.3"},
  {["sp"]          = "2.3.3"},
  {["w3emc"]       = "2.9.2"},
  {["gftl-shared"] = "v1.5.0"},
  {["mapl"]        = "2.22.0-esmf-8.3.0b09"},
  {["nemsio"]      = "2.5.4"},
  {["sigio"]       = "2.3.2"},
  {["sfcio"]       = "1.4.1"},
  {["wrf_io"]      = "1.2.0"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end

setenv("CC","cc")
setenv("FC","ftn")
setenv("CXX","CC")
setenv("CMAKE_Platform","gaea.intel")
