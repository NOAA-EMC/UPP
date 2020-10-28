
# EMC_post

The EMC_post package provides tools for post processing for NCEP models.

This is part of the [NCEPLIBS](https://github.com/NOAA-EMC/NCEPLIBS)
project.

## Authors

NCEP/EMC Developers

Code Manager: WenMeng

## Prerequisites

This package requires the following NCEPLIBS packages:

- [NCEPLIBS-g2](https://github.com/NOAA-EMC/NCEPLIBS-g2)
- [NCEPLIBS-g2tmpl](https://github.com/NOAA-EMC/NCEPLIBS-g2tmpl)
- [NCEPLIBS-sp](https://github.com/NOAA-EMC/NCEPLIBS-sp)
- [NCEPLIBS-ip](https://github.com/NOAA-EMC/NCEPLIBS-ip)
- [NCEPLIBS-bacio](https://github.com/NOAA-EMC/NCEPLIBS-bacio)
- [NCEPLIBS-w3emc](https://github.com/NOAA-EMC/NCEPLIBS-w3emc)
- CRTM

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

## Building

```
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=/path/to/install 
make -j2
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
