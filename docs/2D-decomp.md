# 2-D Decomposition Overview

**Author:** George Vandenberghe

**Date:** June 2022

## Comparison of 1D vs. 2D Decomposition
The 1D decomposition can read state from a model forecast file, either by reading on rank 0 and scattering, or by doing MPI_IO on the model history file using either nemsio, sigio, or netcdf serial or parallel I/O. Very old post tags also implement the more primitive full state broadcast or (a performance bug rectified 10/17) read the entire state on all tasks. This is mentioned in case a very old tag is encountered.  

The 2D decomposition only supports MPI_IO, namely NetCDF Parallel I/O. But the code is backwards compatible and all I/O methods remain supported for the 1D decomposition cases and works for all cases currently supported by older 1D tags and branches.

## 2D Decomposition Design 

The 2D decomposition operates on subdomains with some latitudes and some longitudes.  The subdomains are lon-lat rectangles rather than strips. This means state must be chopped into pieces in any scatter operation and the pieces reassembled in any gather operation that requires a continuous in memory state. I/O and halo exchanges both require significantly more bookkeeping.

The structural changes needed for the 2D decomposition are implemented in MPI_FIRST.f and CTLBLK.f. The CTLBLK.f routine contains numerous additional variables describing left and right domain boundaries. Many additional changes are also implemented in EXCH.f to support 2D halos.  Many additional routines required addition of the longitude subdomain limits but changes to the layouts are handled in CTLBLK.f and the "many additional routines" do not require additional changes when subdomain shapes are changed and have not been a trouble point.

Both MPI_FIRST.f and EXCH.f contain significant additional test code to exchange arrays containing grid coordinates and ensure EXACT matches for all exchanges before the domain exchanges are performed. This is intended to trap errors in the larger variety of 2D decomposition layouts that are possible and most of it can eventually be removed or made conditional at build and run time.

Indices and variables to facilitate the 2D decomposition are found in CTLBLK.f and shared in the rest of UPP through use of CTLBLK.mod.

