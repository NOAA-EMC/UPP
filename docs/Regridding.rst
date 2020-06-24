**********
Regridding
**********

Users that wish to interpolate their unipost output to a different grid
may do so with the *wgrib2* utility. The general format for re-gridding
to a lat-lon grid is given in the example.

==================
Examples of wgrib2
==================

*Wgrib2* is a versatile program that has the ability to convert
grib2 files from one grid to another for various user-defined grids as
well as pre-defined NCEP grids. Complete documentation with examples of
re-gridding for all available grid definitions can be found at:

http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/new_grid.html.

Sample command line usage for calling wgrib2:

   *wgrib2 -new\_grid\_winds W -new\_grid A B C outfile*

Where,

  **W** = earth or grid

          earth: winds oriented to the earths north and south directions

          grid: winds are rotated so that north is relative to the grid

  **A**, **B**, and **C** represent the output grid description

  Sample lat-lon grid description:

  **A** = latlon

  **B** = lon0:nlon:dlon

          lon0 is longitude of first grid point in degrees

          nlon is number of longitudes

          dlon is grid resolution in degrees of longitude

  **C** = lat0:nlat:dlat

          lat0 is latitude of first grid point

          nlat is number of latitudes

          dlat is grid resolution in degrees of latitude

**Note:** *wgrib2 is not distributed within the UFS weather
application. Users may download and install from
http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/.*
