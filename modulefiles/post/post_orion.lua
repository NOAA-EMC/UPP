help([[
Load environment to build post on ORION
]])

prepend_path("MODULEPATH", "/apps/contrib/NCEP/hpc-stack/libs/hpc-stack-gfsv16/modulefiles/stack")

hpc_ver=os.getenv("hpc_ver") or "1.2.0"
load(pathJoin("hpc", hpc_ver))

hpc_intel_ver=os.getenv("hpc_intel_ver") or "2018.4"
load(pathJoin("hpc-intel",hpc_intel_ver))
impi_ver=os.getenv("impi_ver") or "2018.4"
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
g2tmpl_ver=os.getenv("g2tmpl_ver") or "1.9.1"
load(pathJoin("g2tmpl", g2tmpl_ver))
w3nco_ver=os.getenv("w3nco_ver") or "2.4.1"
load(pathJoin("w3nco", w3nco_ver))
bacio_ver=os.getenv("bacio_ver") or "2.4.1"
load(pathJoin("bacio", bacio_ver))
gfsio_ver=os.getenv("gfsio_ver") or "1.4.1"
load(pathJoin("gfsio", gfsio_ver))
ip_ver=os.getenv("ip_ver") or "3.3.3"
load(pathJoin("ip", ip_ver))
sp_ver=os.getenv("sp_ver") or "2.3.3"
load(pathJoin("sp", sp_ver))
crtm_ver=os.getenv("crtm_ver") or "2.4.0"
load(pathJoin("crtm", crtm_ver))
w3emc_ver=os.getenv("w3emc_ver") or "2.9.2"
load(pathJoin("w3emc", w3emc_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.2"
load(pathJoin("nemsio", nemsio_ver))
sigio_ver=os.getenv("sigio_ver") or "2.3.2"
load(pathJoin("sigio", sigio_ver))
sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
load(pathJoin("sfcio", sfcio_ver))
wrf_io_ver=os.getenv("wrf_io_ver") or "1.2.0"
load(pathJoin("wrf_io", wrf_io_ver))

whatis("Description: post build environment")
