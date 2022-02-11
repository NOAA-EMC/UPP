help([[
Load environment to build upp lib on HERA
]])

prepend_path("MODULEPATH", "/apps/contrib/NCEP/libs/hpc-stack/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.1.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2018.4"
load(pathJoin("hpc-intel",hpc_intel_ver))
impi_ver=os.getenv("impi_ver") or "2018.0.4"
load(pathJoin("hpc-impi", impi_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.10.6"
load(pathJoin("hdf5", hdf5_ver))
netcdf_ver=os.getenv("netcdf_ver") or "4.7.4"
load(pathJoin("netcdf", netcdf_ver))

jasper_ver=os.getenv("jasper_ver") or "2.0.25"
load(pathJoin("jasper", jasper_ver))
libpng_ver=os.getenv("libpng_ver") or "1.6.37"
load(pathJoin("libpng", libpng_ver))
zlib_ver=os.getenv("zlib_ver") or "1.2.11"
load(pathJoin("zlib", zlib_ver))

g2_ver=os.getenv("g2_ver") or "3.4.5"
load(pathJoin("g2", g2_ver))
ip_ver=os.getenv("ip_ver") or "3.3.3"
load(pathJoin("ip", ip_ver))
g2tmpl_ver=os.getenv("g2tmpl_ver") or "1.9.1"
load(pathJoin("g2tmpl", g2tmpl_ver))
w3nco_ver=os.getenv("w3nco_ver") or "2.4.1"
load(pathJoin("w3nco", w3nco_ver))
bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))
gfsio_ver=os.getenv("gfsio_ver") or "1.4.1"
load(pathJoin("gfsio", gfsio_ver))
sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))
crtm_ver=os.getenv("crtm_ver") or "2.3.0"
load(pathJoin("crtm", crtm_ver))
w3emc_ver=os.getenv("w3emc_ver") or "2.9.2"
load(pathJoin("w3emc", w3emc_ver))

whatis("Description: upp lib build environment")
