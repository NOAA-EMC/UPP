************
Introduction
************

The Unified Post Processor (UPP) software package is a software package designed to generate useful products from raw model output. The UPP is currently used in operations with the Global Forecast System (GFS), GFS Ensemble Forecast System (GEFS), North American Mesoscale (NAM), Rapid Refresh (RAP), High Resolution Rapid Refresh (HRRR), Short Range Ensemble Forecast (SREF), Hurricane WRF (HWRF) applications, and is also used in Unified Forecasting System (UFS) applications. The UPP provides the capability to compute a variety of diagnostic fields and interpolate to pressure levels or other vertical coordinates. UPP also incorporates the Joint Center for Satellite Data Assimilation (JCSDA) Community Radiative Transfer Model (CRTM) to compute model derived brightness temperature (TB) for various instruments and channels. This additional feature enables the generation of a number of simulated satellite products including GOES products. Output from the UPP is in National Weather Service (NWS) and World Meteorological Organization (WMO) GRIB2 format and can be used directly by visualization, plotting, or verification packages, or for further downstream post-processing, e.g. statistical post-processing techniques.
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

Support for the UFS UPP is provided through the UFS Forum by the Developmental Testbed Center (DTC) for FV3-based applications.
