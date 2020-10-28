
# Unified Post-Processing (UPP)

The Unified Post Processing (UPP) package provides tools for post
processing for NCEP models. For historical reasons, the UPP repository
is called EMC_post.

This is part of the [NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS)
project. It also serves as a standalone package that is distributed to
the community. It is built outside of nceplibs for some UFS
applications (e.g. SRW v1.0.0).

## Authors

NCEP/EMC Developers

Code Manager: Wen Meng, Huiya Chuang, Kate Fossell

## Prerequisites

This package requires the following NCEPLIBS packages:

- [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2)
- [NCEPLIBS-g2tmpl](https://github.com/NOAA-EMC/NCEPLIBS-g2tmpl)
- [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)
- [NCEPLIBS-ip](https://github.com/NOAA-EMC/NCEPLIBS-ip)
- [NCEPLIBS-bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)
- [NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc)
- [NCEPLIBS-w3nco](https://github.com/NOAA-EMC/NCEPLIBS-w3nco)
- [CRTM](https://github.com/noaa-emc/emc_crtm)

Also required to build NCEPpost executable (cmake option
BUILD_POSTEXEC):

- [NCEPLIBS-sigio](https://github.com/NOAA-EMC/NCEPLIBS-sigio) -
- [NCEPLIBS-sfcio](https://github.com/NOAA-EMC/NCEPLIBS-sfcio) -
- [NCEPLIBS-nemsio](https://github.com/NOAA-EMC/NCEPLIBS-nemsio) -
- [NCEPLIBS-gfsio](https://github.com/NOAA-EMC/NCEPLIBS-gfsio)

The [NCEPLIBS-wrf_io](https://github.com/NOAA-EMC/NCEPLIBS-wrf_io)
library is required to build with NCEPpost with WRF-IO library (cmake
option BUILD_WITH_WRFIO).

The following third-party libraries are required:

- [netcdf-c](https://github.com/Unidata/netcdf-c)
- [netcdf-fortran](https://github.com/Unidata/netcdf-fortran)
- [Jasper](https://github.com/jasper-software/jasper)
- [libpng](http://www.libpng.org/pub/png/libpng.html)
- [libz](https://zlib.net/)

## Building

Builds include:

- Operational use GNC build as Wen described for both library and
  executable (library used for GFS only at this time)

- MRW App uses UPP packaged with nceplibs and cmake to build/run with
  executable (via release/public-v1 branch).

- SRW App uses UPP repo branch/tag directly and uses cmake to
  build/run with executable (via release/public-v2 branch).

- Community standalone uses UPP repo branch/tag directly and uses
  cmake to build/run with executable (via release/public-v2
  branch). For these procedures, we add a
  -DCMAKE_PREFIX_PATH=${INSTALL_PREFIX} where INSTALL_PREFIX is the
  location of the nceplibs installation as a dependency requirement.

```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install 
make 
make test
make install
```

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
