-- ---------------------------------------------------------------------------
-- David Huber 06/2021, Set up config. with the hpc-stack NCEPLIBS.
-- Innocent Souopgui 11/2023, Update to use spack-stack
-- ---------------------------------------------------------------------------

help([[
Loads modules required for building upp
]])

prepend_path("MODULEPATH", "/data/prod/jedi/spack-stack/spack-stack-1.5.0/envs/unified-env/install/modulefiles/Core")


local wrf_io_ver=os.getenv("wrf_io_ver") or "1.2.0"
local stack_intel_ver=os.getenv("stack_intel_ver") or "2021.5.0"
local stack_impi_ver=os.getenv("stack_impi_ver") or "2021.5.0"
local cmake_ver=os.getenv("cmake_ver") or "3.23.1"
local jasper_ver=os.getenv("jasper_ver") or "2.0.32"
local png_ver=os.getenv("png_ver") or "1.6.37"
local zlib_ver=os.getenv("zlib_ver") or "1.2.13"
local netcdf_c_ver=os.getenv("netcdf_c_ver") or "4.9.2"
local netcdf_fortran_ver=os.getenv("netcdf_fortran_ver") or "4.6.0"
local bacio_ver=os.getenv("bacio_ver") or "2.4.1"
local crtm_ver=os.getenv("crtm_ver") or "2.4.0"
local g2_ver=os.getenv("g2_ver") or "3.4.5"
local g2tmpl_ver=os.getenv("g2tmpl_ver") or "1.10.2"
local ip_ver=os.getenv("ip_ver") or "4.3.0"
local nemsio_ver=os.getenv("nemsio_ver") or "2.5.4"
local sfcio_ver=os.getenv("sfcio_ver") or "1.4.1"
local sigio_ver=os.getenv("sigio_ver") or "2.3.2"
local w3emc_ver=os.getenv("w3emc_ver") or "2.10.0"
local sp_ver=os.getenv("sp_ver") or "2.3.3"

load(pathJoin("stack-intel", stack_intel_ver))
load(pathJoin("stack-intel-oneapi-mpi", stack_impi_ver))

load(pathJoin("cmake", cmake_ver))

load(pathJoin("jasper", jasper_ver))
load(pathJoin("zlib", zlib_ver))
load(pathJoin("libpng", png_ver))

load(pathJoin("netcdf-c", netcdf_c_ver))
load(pathJoin("netcdf-fortran", netcdf_fortran_ver))

load(pathJoin("bacio", bacio_ver))
load(pathJoin("crtm", crtm_ver))
load(pathJoin("g2", g2_ver))
load(pathJoin("g2tmpl", g2tmpl_ver))
load(pathJoin("ip", ip_ver))
load(pathJoin("nemsio", nemsio_ver))
load(pathJoin("sfcio", sfcio_ver))
load(pathJoin("sigio", sigio_ver))
load(pathJoin("sp", sp_ver))
load(pathJoin("w3emc", w3emc_ver))
load(pathJoin("wrf-io", wrf_io_ver))

