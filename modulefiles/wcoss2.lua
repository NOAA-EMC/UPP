help([[
Load environment to build post on WCOSS2
]])

PrgEnv_intel_ver=os.getenv("PrgEnv_intel_ver") or "8.3.3"
intel_ver=os.getenv("intel_ver") or "19.1.3.304"
craype_ver=os.getenv("craype_ver") or "2.7.17"
cray_mpich_ver=os.getenv("cray_mpich_ver") or "8.1.19"
load(pathJoin("PrgEnv-intel", PrgEnv_intel_ver))
load(pathJoin("intel", intel_ver))
load(pathJoin("craype", craype_ver))
load(pathJoin("cray-mpich", cray_mpich_ver))

hdf5_ver=os.getenv("hdf5_ver") or "1.14.0"
netcdf_ver=os.getenv("netcdf_ver") or "4.9.2"
load(pathJoin("hdf5-B", hdf5_ver))
load(pathJoin("netcdf-B", netcdf_ver))

jasper_ver=os.getenv("jasper_ver") or "2.0.25"
libpng_ver=os.getenv("libpng_ver") or "1.6.37"
zlib_ver=os.getenv("zlib_ver") or "1.2.11"
load(pathJoin("jasper", jasper_ver))
load(pathJoin("libpng", libpng_ver))
load(pathJoin("zlib", zlib_ver))

g2_ver=os.getenv("g2_ver") or "3.4.5"
g2tmpl_ver=os.getenv("g2tmpl_ver") or "1.10.2"
bacio_ver=os.getenv("bacio_ver") or "2.4.1"
ip_ver=os.getenv("ip_ver") or "4.0.0"
sp_ver=os.getenv("sp_ver") or "2.3.3"
crtm_ver=os.getenv("crtm_ver") or "2.4.0"
w3emc_ver=os.getenv("w3emc_ver") or "2.9.2"
load(pathJoin("g2", g2_ver))
load(pathJoin("g2tmpl", g2tmpl_ver))
load(pathJoin("bacio", bacio_ver))
load(pathJoin("ip", ip_ver))
load(pathJoin("sp", sp_ver))
load(pathJoin("crtm", crtm_ver))
load(pathJoin("w3emc", w3emc_ver))

nemsio_ver=os.getenv("nemsio_ver") or "2.5.2"
sigio_ver=os.getenv("sigio_ver") or "2.3.2"
sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
wrf_io_ver=os.getenv("wrf_io_ver") or "1.2.0"
load(pathJoin("nemsio", nemsio_ver))
load(pathJoin("sigio", sigio_ver))
load(pathJoin("sfcio", sfcio_ver))
load(pathJoin("wrf_io", wrf_io_ver))

setenv("CC","cc")
setenv("CXX","CC")
setenv("FC","ftn")

whatis("Description: post build environment")
