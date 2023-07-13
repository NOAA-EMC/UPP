
# Unified Post-Processing (UPP)

The Unified Post Processor (UPP) software package is a software
package designed to generate useful products from raw model
output.

The UPP is currently used in operations with the Global Forecast
System (GFS), GFS Ensemble Forecast System (GEFS), North American
Mesoscale (NAM), Rapid Refresh (RAP), High Resolution Rapid Refresh
(HRRR), Short Range Ensemble Forecast (SREF), and Hurricane WRF (HWRF)
applications. It is also used in the Unified Forecasting System (UFS),
including the Rapid Refresh Forecast System (RRFS), Hurricane Application
Forecasting System (HAFS), and the Medium Range Weather (MRW) and Short 
Range Weather (SRW) Applications.

The UPP provides the capability to compute a variety of diagnostic
fields and interpolate to pressure levels or other vertical
coordinates.

UPP also incorporates the Joint Center for Satellite Data Assimilation
(JCSDA) Community Radiative Transfer Model (CRTM) to compute model
derived brightness temperature (TB) for various instruments and
channels. This additional feature enables the generation of a number
of simulated satellite products including GOES products.

Output from the UPP is in National Weather Service (NWS) and World
Meteorological Organization (WMO) GRIB2 format and can be used
directly by visualization, plotting, or verification packages, or for
further downstream post-processing, e.g. statistical post-processing
techniques.

Examples of UPP products include:

- T, Z, humidity, wind, cloud water, cloud ice, rain, and snow on pressure levels
- SLP, shelter level T, humidity, and wind fields
- Precipitation-related fields
- PBL-related fields
- Severe weather products (e.g. CAPE, Vorticity, Wind shear)
- Radiative/Surface fluxes
- Cloud related fields
- Aviation products
- Radar reflectivity products
- Satellite look-alike products


## User Support
Support for the UFS UPP is provided through [GitHub Discussions](https://github.com/NOAA-EMC/UPP/discussions).

## Documentation 
User Guide for latest public release: https://upp.readthedocs.io/en/latest/.

Technical code-level documentation: https://noaa-emc.github.io/UPP/.

## Developer Information
Please see review the [wiki](https://github.com/NOAA-EMC/UPP/wiki)

## Authors

NCEP/EMC Developers

Code Managers: Wen Meng, Huiya Chuang, Kate Fossell

## Prerequisites

The UPP requires certain NCEPLIBS packages to be installed via the 
HPC-Stack project. For instructions on installing these packages as a 
bundle via HPC-Stack, see: https://hpc-stack.readthedocs.io/en/latest/.
Users may instead install packages via spack-stack. For instructions,
see: https://spack-stack.readthedocs.io/en/latest/.
The `UPP/modulefiles` directory indicates which package versions are 
used and supported on Level 1 systems. 

Required NCEPLIBS packages:

- [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2)
- [NCEPLIBS-g2tmpl](https://github.com/NOAA-EMC/NCEPLIBS-g2tmpl)
- [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)
- [NCEPLIBS-ip](https://github.com/NOAA-EMC/NCEPLIBS-ip)
- [NCEPLIBS-bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)
- [NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc)
- [CRTM](https://github.com/noaa-emc/emc_crtm)

Also required to build NCEPpost executable (cmake option
BUILD_POSTEXEC):

- [NCEPLIBS-sigio](https://github.com/NOAA-EMC/NCEPLIBS-sigio)
- [NCEPLIBS-sfcio](https://github.com/NOAA-EMC/NCEPLIBS-sfcio)
- [NCEPLIBS-nemsio](https://github.com/NOAA-EMC/NCEPLIBS-nemsio)
- [NCEPLIBS-gfsio](https://github.com/NOAA-EMC/NCEPLIBS-gfsio)

The [NCEPLIBS-wrf_io](https://github.com/NOAA-EMC/NCEPLIBS-wrf_io)
library is required to build with NCEPpost with WRF-IO library (cmake
option BUILD_WITH_WRFIO).

The following third-party libraries are required:

- [netcdf](https://github.com/Unidata/netcdf)
- [netcdf-c](https://github.com/Unidata/netcdf-c)
- [netcdf-fortran](https://github.com/Unidata/netcdf-fortran)
- [Jasper](https://github.com/jasper-software/jasper)
- [libpng](http://www.libpng.org/pub/png/libpng.html)
- [zlib](https://zlib.net/)
- [hdf5](https://github.com/HDFGroup/hdf5)

## Building

Builds include:

- Inline post (UPP library): Currently only supported for the GFS, RRFS,
  HAFS, and the UFS-MRW Application.

- Offline post (UPP executable): Supported for Regional applications
  including SRW, RRFS, HAFS, and standalone applications of UPP.


CMake is used to manage all builds of the UPP. 
The script `UPP/tests/compile_upp.sh` can be used to automatically
build UPP on fully supported platforms where HPC-stack is supported.
Details in this script can be used to build on new platforms


## Disclaimer

The United States Department of Commerce (DOC) GitHub project code is
provided on an "as is" basis and the user assumes responsibility for
its use. DOC has relinquished control of the information and no longer
has responsibility to protect the integrity, confidentiality, or
availability of the information. Any claims against the Department of
Commerce stemming from the use of its GitHub project will be governed
by all applicable Federal law. Any reference to specific commercial
products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of
Commerce. The Department of Commerce seal and logo, or the seal and
logo of a DOC bureau, shall not be used in any manner to imply
endorsement of any commercial product or activity by DOC or the United
States Government.

## UPP Terms of Use Notice

The UPP Terms of Use Notice is available at: https://github.com/NOAA-EMC/UPP/wiki/UPP-Terms-of-Use-Notice
