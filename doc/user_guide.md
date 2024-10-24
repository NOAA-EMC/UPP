@mainpage

# Unified Post-Processing (UPP)

The UPP code is publicly available on GitHub at https://github.com/NOAA-EMC/UPP.

## Documentation for Previous Versions of UPP

* [UPP Version UPP-SRW-v2.2.0](upp-srw-v2.2.0/index.html)
* [UPP Version 11.0.0](upp_v11.0.0/index.html)

## Background Information

The Unified Post Processor (UPP) software package is a software
package designed to generate useful products from raw model
output. The UPP is currently used in operations with the Global
Forecast System (GFS), GFS Ensemble Forecast System (GEFS), North
American Mesoscale (NAM), Rapid Refresh (RAP), High Resolution Rapid
Refresh (HRRR), Short Range Ensemble Forecast (SREF), Hurricane WRF
(HWRF) applications, and is also used in Unified Forecasting System
(UFS) applications. The UPP provides the capability to compute a
variety of diagnostic fields and interpolate to pressure levels or
other vertical coordinates. UPP also incorporates the Joint Center for
Satellite Data Assimilation (JCSDA) Community Radiative Transfer Model
(CRTM) to compute model derived brightness temperature (TB) for
various instruments and channels. This additional feature enables the
generation of a number of simulated satellite products including GOES
products. Output from the UPP is in National Weather Service (NWS) and
World Meteorological Organization (WMO) GRIB2 format and can be used
directly by visualization, plotting, or verification packages, or for
further downstream post-processing, e.g. statistical post-processing
techniques. Examples of UPP products include:

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

## Prerequisites
The UPP requires certain NCEPLIBS packages to be installed via the spack-stack project. For instructions on installing these packages as a bundle via spack-stack, see: https://spack-stack.readthedocs.io/en/latest/. The UPP/modulefiles directory indicates which package versions are used and supported on Level 1 systems.

## Community Support
Community support for the Unified Forecast System (UFS) UPP in FV3-based applications is provided by the
Earth Prediction Innovation Center (EPIC). Community support for the UPP with WRF is no longer available. 

* Support for the UFS UPP is provided through [GitHub Discussions](https://github.com/NOAA-EMC/UPP/discussions).
* The UPP User's Guide for the latest standalone public release is [UPP v11.0.0](https://upp.readthedocs.io/en/upp_v11.0.0/).
* The UPP User's Guide for develop branch is [UPP develop](https://upp.readthedocs.io/en/develop/).
* The [UPP wiki](https://github.com/NOAA-EMC/UPP/wiki) includes relevant information and links for users and developers. 
* Instructions on technical code documentation are available in a set of [Doxygen Documentation Slides](https://github.com/NOAA-EMC/UPP/wiki/DoxygenDocumentation.pdf).

Code Managers: Wen Meng (EMC), Huiya Chuang (EMC), Fernando Andrade-Maldonado (EPIC)

