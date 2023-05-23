**********
Regridding
**********

Users that wish to interpolate their unipost output to a different grid may do so with the *wgrib2*
utility. The general format for re-gridding to various common projections are outlined in the following
examples.

*Wgrib2* is a versatile program that has the ability to convert grib2 files from one grid to another
for various user-defined grids as well as pre-defined NCEP grids. Complete documentation with examples
of re-gridding for all available grid definitions can be found at:

http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/new_grid.html

.. _Examples-of-wgrib2:

==================
Examples of wgrib2
==================

**Example 1: Latitude-Longitude Grid**

*-new_grid latlon lon0:nlon:dlon lat0:nlat:dlat outfile*

+----------+------------------------------------------+
| Variable | Description                              |
+==========+==========================================+
| lon0     | Longitude of first grid point in degrees |
+----------+------------------------------------------+
| nlon     | Number of longitudes                     |
+----------+------------------------------------------+
| dlon     | Grid resolution in degrees of longitude  |
+----------+------------------------------------------+
| lat0     | Latitude of first grid point in degrees  |
+----------+------------------------------------------+
| nlat     | Number of latitudes                      |
+----------+------------------------------------------+
| dlat     | Grid resolution in degrees of latitude   |
+----------+------------------------------------------+

**Example 2: Lambert Conic Conformal Grid**

*-new_grid lambert:lov:latin1:latin2 lon0:nx:dx lat0:ny:dy outfile*

+----------+-----------------------------------------------------------------+
| Variable | Description                                                     |
+==========+=================================================================+
| lov      | Longitude where y axis is parallel to meridian in degrees       |
+----------+-----------------------------------------------------------------+
| latin1   | First latitude from pole which cuts the secant cone in degrees  |
+----------+-----------------------------------------------------------------+
| latin2   | Second latitude from pole which cuts the secant cone in degrees |
+----------+-----------------------------------------------------------------+
| lon0     | Longitude of the first grid point in degrees                    |
+----------+-----------------------------------------------------------------+
| nx       | Total number of grid points along x                             |
+----------+-----------------------------------------------------------------+
| dx       | Grid cell size in meters in x direction                         |
+----------+-----------------------------------------------------------------+
| lat0     | Latitude of the first grid point in degrees                     |
+----------+-----------------------------------------------------------------+
| ny       | Total number of grid points along y                             | 
+----------+-----------------------------------------------------------------+
| dy       | Grid cell size in meters in y direction                         |
+----------+-----------------------------------------------------------------+

**Example 3: Polar Stereographic Grid**

*-new_grid nps(or SPS):lov:lad lon0:nx:dx lat0:ny:dy outfile*

+----------+-----------------------------------------------------------+
| Variable | Description                                               |
+==========+===========================================================+
| nps/sps  | North/south polar stereographic                           |
+----------+-----------------------------------------------------------+
| lov      | Longitude where y axis is parallel to meridian in degrees |
+----------+-----------------------------------------------------------+
| lad      | Latitude where dx and dy are specified                    |
+----------+-----------------------------------------------------------+
| lon0     | Longitude of the first grid point in degrees              |
+----------+-----------------------------------------------------------+
| nx       | Total number of grid points along x                       |
+----------+-----------------------------------------------------------+
| dx       | Grid cell distance in meters in x direction at lad        |
+----------+-----------------------------------------------------------+
| lat0     | Latitude of the first grid point in degrees               |
+----------+-----------------------------------------------------------+
| ny       | Total number of grid points along y                       |
+----------+-----------------------------------------------------------+
| dy       | Grid cell distance in meters in y direction at lad        |
+----------+-----------------------------------------------------------+

**Winds**

*-new_grid_winds grid(or earth)*

+----------+----------------------------------------------+
| Variable | Description                                  |
+==========+==============================================+
| grid     | U-wind goes from grid (i,J) to (i+1,j)       |
+----------+----------------------------------------------+
| earth    | U-wind goes eastward, V-wind goes northward  |
+----------+----------------------------------------------+

**Interpolation**

The default interpolation type is bilinear, but it can be set to another type (e.g. neighbor, budget).
 
*-new_grid_interpolation type*

**Operational Example**

Interpolates to a 0.25 degree latitude-longitude grid using various interpolation types depending on
the variable.

.. code-block:: console

    wgrib2 infile -set_grib_type same -new_grid_winds earth |
    -new_grid_interpolation bilinear |
    -if ":(CRAIN|CICEP|CFRZR|CSNOW|ICSEV):" -new_grid_interpolation neighbor -fi |
    -set_bitmap 1 -set_grib_max_bits 16 |
    -if ":(APCP|ACPCP|PRATE|CPRAT):" -set_grib_max_bits 25 -fi |
    -if ":(APCP|ACPCP|PRATE|CPRAT|DZDT):" -new_grid_interpolation budget -if |
    -new_grid "latlon 0:1440:0.25 90:721:-0.25" outfile

**Note:** *wgrib2 is not distributed within the UFS weather application. Users may download and install
from http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/.*
