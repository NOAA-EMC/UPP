whatis("Description: UPP build environment common libraries")

help([[Load UFS Model common libraries]])

local ufs_modules = {
  {["jasper"]          = "2.0.32"},
  {["zlib"]            = "1.2.13"},
  {["libpng"]          = "1.6.37"},
  {["hdf5"]            = "1.14.0"},
  {["netcdf-c"]        = "4.7.4"},
  {["netcdf-fortran"]  = "4.5.3"},
  {["parallelio"]      = "2.5.10"},
  {["bacio"]           = "2.4.1"},
  {["crtm"]            = "2.3.0"},
  {["g2"]              = "3.4.5"},
  {["g2tmpl"]          = "1.10.0"},
  {["ip"]              = "3.3.3"},
  {["sp"]              = "2.3.3"},
  {["w3emc"]           = "2.9.2"},
  {["nemsio"]          = "2.5.4"},
  {["sigio"]           = "2.3.2"},
  {["sfcio"]           = "1.4.1"},
  {["wrf-io"]          = "1.2.0"},
}

for i = 1, #ufs_modules do
  for name, default_version in pairs(ufs_modules[i]) do
    local env_version_name = string.gsub(name, "-", "_") .. "_ver"
    load(pathJoin(name, os.getenv(env_version_name) or default_version))
  end
end
