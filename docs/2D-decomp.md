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

The following is found in CTLBLK.f and shared in the rest of UPP through use of CTLBLK.mod:

| Variable | Type | Description |
|---------|------|-------------|
| im | integer |  full longitude domain|
| jm | integer |  full latitude domain|
| | | | 
| jsta | integer |  start latitude on a task subdomain|
| jend | integer |   end  latitude on a task subdomain|
| ista | integer  | start longitude on a task subdomain|
| iend  | integer |  end   longitude on a task subdomain|
| | | | 
| ista_2l | integer |start longitude -2 of the subdomain|
| iend_2u | integer |end   longitude +2 of the subdomain|
| jsta_2l | integer |start latitude  -2 of the subdomain|
| jend_2u | integer |end   latitude  +2 of the subdomain|

The shape of the subdomain is ista_2l:iend_2u,jsta_2l:jend_2u so it includes the halos although the halos are not populated until exchange is done in EXCH.f.  Because of halos we need more bounds defined:

| Variable | Type | Description |
|---------|------|-------------|
| jsta_m |  integer |  Beginning latitude loop index in subdomain for halo depth 1 |
| jend_m |   integer |   ending latitude loop index in subdomain for halo depth 1      |
| jsta_m2 |  integer |   second latitude below begin latitude of subdomain for halo depth 2 (in NGMFLD.f) |
| jend_m2 |   integer |  second latitude above end latitude of subdomain for halo depth 2 ( in NGMFLD.f) |

Note:<ul><li>In interior subdomains these are the same as jsta and jend.</li><li>In boundary subdomains these loop indices define a smaller subset of the subdomain since halos are not defined on the full domain boundaries and stencils must be restricted to valid full domain points. </li></ul>

| Variable | Type | Description |
|---------|------|-------------| 
|  ista_m |  integer  |  begining longitude loop index in subdomain for halo depth 1|
|  iend_m |  integer  |  end longitude loop index in subdomain  for halo depth 1 |
|  ista_m2 |  integer  |  second longitude before begin longitude for halo depth 2  (not used as of 6/22)|
|  iend_m2 |  integer  | second longitude after end  longitude    for halo depth 2  (not used as of 6/22) |

Note:<ul><li>In interior subdomains these are the same as ista and iend.</li><li>In boundary subdomains these loop indices define a smaller subset of the subdomain since halos are not defined on the full domain boundaries and stencils must be restricted to valid full domain points.</li></ul>


| Variable | Type | Description |
|---------|------|-------------| 
|  ileft |  integer |  MPI rank containing the last longitude before ista 
|  iright |  integer |  MPI rank containing the first longitude after iend
|  iup   |  integer |   MPI rank containing the first latitude after jend
|  idn |  integer |     MPI rank containing the  last latitude before  jsta
| | | |
|  ileftb  | integer |MPI rank containing the last longitude before ista but for cyclic boundary conditions where "last" at the beginning is the other end of the domain (apparently unused and replaced with local calculation) |
|  irightb | integer |  MPI rank containing the first longitude after iend but for cyclic boundary conditions where "first" at the beginning is the other end of the domain (apparently unused and replaced with local calculation) |

