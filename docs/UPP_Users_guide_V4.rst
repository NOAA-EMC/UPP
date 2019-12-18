.. role:: math(raw)
   :format: html latex
.. role:: bolditalic
   :class: bolditalic

======================================================
User’s Guide for the NCEP Unified Post Processor (UPP)
======================================================

**Version 4**

*Acknowledgments:*

*The adaptation of the original WRF Post Processor package and Users
Guide (by Mike Baldwin of NSSL/CIMMS and Hui-Ya Chuang of NCEP/EMC) was
done by Lígia Bernardet (NOAA/ESRL/DTC) in collaboration with Dusan
Jovic (NCEP/EMC), Robert Rozumalski (COMET), Wesley Ebisuzaki
(NWS/HQTR), and Louisa Nance (NCAR/RAL/DTC). Upgrades to WRF Post
Processor versions 2.2 and higher were performed by Hui-Ya Chuang, Dusan
Jovic and Mathew Pyle (NCEP/EMC). Transitioning of the documentation
from the WRF Post Processor to the Unified Post Processor was performed
by Nicole McKee (NCEP/EMC), Hui-ya Chuang (NCEP/EMC), and Jamie Wolff
(NCAR/RAL/DTC). Implementation of the Community Unified Post Processor
was performed by Tricia Slovacek and Kate Fossell (NCAR/RAL/DTC).*

NCEP Unified Post Processor (UPP)
=================================

UPP Introduction
================

The Unified Post Processor (UPP) software package is a software package 
designed to generate useful products from raw model output. The UPP is 
currently used in operations with the Global Forecast System (GFS), GFS 
Ensemble Forecast System (GEFS), North American Mesoscale (NAM), Rapid 
Refresh (RAP), High Resolution Rapid Refresh (HRRR), Short Range Ensemble 
Forecast (SREF),  Hurricane WRF (HWRF) applications, and is also used in 
Unified Forecasting System (UFS) applications.  The UPP provides the capability 
to compute a variety of diagnostic fields and interpolate to pressure levels 
or other vertical coordinates. UPP also incorporates the Joint Center for 
Satellite Data Assimilation (JCSDA) Community Radiative Transfer Model (CRTM) 
to compute model derived brightness temperature (TB) for various instruments 
and channels. This additional feature enables the generation of a number of 
simulated satellite products including GOES products.  Output from the UPP is 
in National Weather Service (NWS) and World Meteorological Organization (WMO) 
GRIB2 format and can be used directly by visualization, plotting, or verification 
packages, or for further downstream post-processing, e.g. statistical 
post-processing techniques.

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

Community user support is currently provided for FV3-based applications, and on 
a limited basis, WRF-ARW applications.

http://www.dtcenter.org/upp/users/docs/user_guide/crtm_ug/CRTM_User_Guide.pdf.

Software Installation/Getting Started
=====================================

System Requirements, Libraries, and Compilers
---------------------------------------------
UPP has been tested on LINUX platforms (with PGI, Intel and GFORTRAN
compilers).

Before installing the UPP code, it is necessary to ensure that you have
the required libraries available on your system. These libraries
include:

-  Unidata’s NetCDF library

   https://www.unidata.ucar.edu/software/netcdf/docs/getting_and_building_netcdf.html
   
-  The NCEP libraries for the UPP application

   https://github.com/NCAR/NCEPlibs
   
-  The JasPer, PNG, and Zlib libraries, required for GRIB2 capabilities. NCEP
   provides these necessary codes for download:

   http://www.nco.ncep.noaa.gov/pmb/codes/GRIB2/
   
-  GNU M4, used to build the WRF I/O libraries included in UPPV4.0.
   This is usually bundled with an operating system or can be downloaded:

   https://www.gnu.org/software/m4/
   
The UPP has some sample visualization scripts included to create graphics using:

-  GrADS (http://cola.gmu.edu/grads/gadoc/gadoc.php)

-  GEMPAK (http://www.unidata.ucar.edu/software/gempak/index.html)

**Note:** *These are not part of the UPP installation and need to be
installed separately if one would like to use either plotting package.*

The UPP and all libraries required by the UPP need to be built with the same compiler.

Obtaining the UPP Code
----------------------

The most recent release UPPV4.0.1 is available on github. To obtain the code, move to the diretory where you would like the repository to reside and enter on the command line:

   git clone -b DTC_post --recurse-submodules https://github.com/NOAA-EMC/EMC_post

This will clone the DTC_post branch of EMC_post, including the CRTM submodule.

UPP Directory Structure
-----------------------

Under the main directory of **UPPV4.0** reside seven subdirectories
(\* indicates directories that are created after the configuration
step):

**docs**:

**exec**\*: Contains the *unipost* executable after successful compilation

**include**\*: Source include modules built/used during compilation of UPP

**lib**\*: Libraries built/used by UPP not in NCEPlibs

**parm**: Contains the parameter files, which can be modified by the
user to control how the post processing is performed.

**scripts**: contains sample running scripts to process model history files.

-  **run\_unipost**: runs *unipost*.

**sorc**: Contains source codes for:

-  **arch**: Machine dependent configuration build scripts used to construct *configure.upp*

-  **ncep_post.fd**: Source code for *unipost*

-  **comlibs**: Contains source code subdirectories for the UPP libraries not included in NCEPlibs

   -  **crtm2**: Community Radiative Transfer Model library

   -  **wrfmpi\_stubs**: Contains some *C* and *FORTRAN* codes
      to generate *libmpi.a* library used to replace MPI calls for
      serial compilation.

   -  **xml**: XML support for the GRIB2 parameter file

Installing the UPP Code
-----------------------

Environment variables to compatible versions of the required libraries must be set
before beginning the installation.

To reference the netCDF libraries, the configure script checks for an
environment variable (**NETCDF**) first, then the system default
(**/user/local/netcdf**), and then a user supplied link
(*./netcdf\_links*). If none of these resolve a path, the user will
be prompted by the configure script to supply a path.

To reference the NCEP libraries, the configure script will also check that you set the environment variable **NCEPLIBS_DIR** to the location where the NCEPlibs were built.

Type configure, and provide the required info. For example:

   *./configure*

You will be given a list of choices for your computer.

1. Linux x86\_64, PGI compiler (serial)

2. Linux x86\_64, PGI compiler (dmpar)

3. Linux x86\_64, Intel compiler (serial)

4. Linux x86\_64, Intel compiler (dmpar)

5. Linux x86\_64, Intel compiler, SGI MPT (serial)

6. Linux x86\_64, Intel compiler, SGI MPT (dmpar)

7. Linux x86\_64, gfortran compiler (serial)

8. Linux x86\_64, gfortran compiler (dmpar)

Choose one of the configure options listed. Check the
*configure.upp* file created and edit for compile options/paths, if
necessary. For debug flag settings, the configure script can be run with
a *-d* switch or flag.

To compile UPP, enter the following command:

   *./compile >& compile\_upp.log &*

When compiling UPP with distributed memory (serial) this command should
create 2 (3) libraries in **EMC_post/lib/** (*libCRTM.a*, (*libmpi.a*), *libxmlparse.a*) and the UPP executable in **exec/** (*unipost.exe*).

To remove all built files, as well as the *configure.upp*, type:

   *./clean*

This action is recommended if a mistake is made during the installation
process or a change is made to the configuration or build environment.
There is also a *clean -a* option which will revert back to a
pre-install configuration.

UPP Functionalities
-------------------

The UPP,

-  is compatible with WRF v3.7 and higher for Ferrier physics.

-  can be used to post-process WRF-ARW, WRF-NMM, NMMB, GFS, CFS, and FV3
   forecasts (community support only provided for WRF-ARW and global FV3
   forecasts).

-  can ingest WRF history files (wrfout\*) in netCDF format.

-  can ingest FV3 history files (dyn\*/phy\* or gfs\*) in binarynemsiompioo.
   format.

**Unipost**

   -  Interpolates the forecasts from the models native vertical
      coordinate to NWS standard output levels (e.g., pressure, height)
      and computes mean sea level pressure. If the requested parameter
      is on a models native level, then no vertical interpolation is
      performed.

   -  Computes diagnostic output quantities (e.g., convective available
      potential energy, helicity, relative humidity). A full list of
      fields that can be generated by *unipost* is provided in
      https://dtcenter.org/upp/users/docs/tables/UPP_GRIB2_Table.pdf.

   -  Outputs the results in NWS and WMO standard GRIB2 format
      (for GRIB documentation, see
      http://www.nco.ncep.noaa.gov/pmb/docs/).

   -  Except for new capabilities of post processing GFS/CFS and
      additions of many new variables, UPP uses the same algorithms to
      derive most existing variables as were used in WPP. The only three
      exceptions/changes from the WPP are:

      -  Computes RH w.r.t. ice for GFS, but w.r.t. water for all other
         supported models. WPP computed RH w.r.t. water only.

      -  The height and wind speed at the maximum wind level is computed
         by assuming the wind speed varies quadratically in height in
         the location of the maximum wind level. The WPP defined maximum
         wind level at the level with the maximum wind speed among all
         model levels. The static tropopause level is obtained by
         finding the lowest level that has a temperature lapse rate of
         less than 2 K/km over a 2 km depth above it. The WPP defined
         the tropopause by finding the lowest level that has a mean UPP
         V3: Users Guide 8 temperature lapse rate of 2 K/km over three
         model layers.

Setting up the model to interface with UPP
-----------------------------------------------------

The unipost program is currently set up to read a large number of fields
from the model history files. This configuration stems from
NCEP’s need to generate all of its required operational products. When
using the netCDF or NEMS binary read, this program is configured such
that it will run successfully even if an expected input field is missing
from the history file as long as this field is not required
to produce a requested output field. If the pre-requisites for a
requested output field are missing from the history file,
unipost will abort at run time.

Take care not to remove fields from the model history files, which may
be needed for diagnostic purposes by the UPP package. For example, if
isobaric state fields are requested, but the pressure fields on model
interfaces are not available in the history
file, unipost will abort at run time. In general, the default fields
available in the model history files are sufficient to run UPP.

UPP is written to process a single forecast hour, therefore, having a
single forecast per output file is optimal. However, for WRF based
forecasts, UPP can be run across multiple forecast times in a single
output file to extract a specified forecast hour.

UPP Control File Overview
=========================

GRIB2 control file
------------------

-  For outputting GRIB2 format using version 4.0, a preprocessing step
   is required by the user to convert the xml file
   *parm/postcntrl.xml* to a flat text file
   *parm/postxconfig-NT.txt*. The flat file is quicker to process
   than the xml file. The user will first need to edit the
   *postcntrl.xml* file to declare which fields are to be output
   from UPP.

-  In order to ensure that the user-edited xml files are error free, XML
   stylesheets (*parm/EMC\_POST\_CTRL\_Schema.xsd* and
   *EMC\_POST\_Avblflds\_Schema.xsd*) are used to validate both the
   *postcntrl.xml* and *post\_avblflds.xml* files, respectively.
   Confirmation of validation will be given (e.g. postcntrl.xml
   validates) or otherwise return errors if it does not match the
   schema. To run the validation:

   *xmllint --noout --schema EMC\_POST\_CTRL\_Schema.xsd
   postcntrl.xml*

   *xmllint --noout --schema EMC\_POST\_Avblflds\_Schema.xsd
   post\_avblflds.xml*

-  Once the xmls are validated, the user will need to generate the flat
   file. Edit the *parm/makefile* if necessary to point to the
   desired flat file directory and xmls. The makefile will call the perl
   program *parm/POSTXMLPreprocessor.pl* to generate the post flat
   file *postxconfig-NT.txt*. Generate the flat file:

   *make*

Controlling which variables unipost outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To output a field, the body of the control file needs to contain an
entry for the appropriate variable. For variables found on isobaric or
height levels, the desired levels to be output must be listed (see next
section: *Controlling which levels unipost outputs*). If an entry for
a particular field is not yet available in the control file, it  may be
added to the control file with the appropriate entries for that field.

Controlling which levels unipost outputs
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The <level> tag in the postcntrl.xml is used to list the desired levels
for output. The following levels are currently available for output:

-  For isobaric output, 46 levels are possible, from 2 to 1000 hPa (*2,
   5, 7, 10, 20, 30, 50, 70 mb and then every 25 mb from 75 to 1000
   mb*). The complete list of levels is specified in
   *src/unipost/CTLBLK.f*.

   -  Modify specification of variable LSMDEF to change the number of
      pressure levels: LSMDEF=47

   -  Modify specification of SPLDEF array to change the values of
      pressure levels: (/200.,500.,700.,1000.,2000.,3000.
      &,5000.,7000.,7500.,10000.,12500.,15000.,17500.,20000., …/)

-  For model-level output, all model levels are possible, from the
   highest to the lowest.

-  When using the Noah LSM, the soil layers are 0-10 cm, 10-40 cm,
   40-100 cm, and 100-200 cm.

-  When using the RUC LSM, the soil levels are 0 cm, 1 cm, 4 cm, 10 cm,
   30 cm, 60 cm, 100 cm, 160 cm, and 300 cm. (For the old RUC LSM,
   there are only 6 layers and if using this, you will need to change
   “RUC LSM” from 9 to 6 in the WRFPOST.f routine.)

-  When using Pliem-Xiu LSM, there are two layers: 0-1 cm, 1-100 cm

-  For low, mid, and high cloud layers, the layers are :math:`\geq`\ 642
   hPa, :math:`\geq`\ 350 hPa, and <350 hPa, respectively.

-  For PBL layer averages, the levels correspond to 6 layers with a
   thickness of 30 hPa each.

-  For flight level, the levels are 30 m, 50 m, 80 m, 100 m, 305 m, 457
   m, 610 m, 914 m, 1524 m, 1829 m, 2134 m, 2743 m, 3658 m, 4572 m, 6000
   m, 7010 m.

-  For AGL radar reflectivity, the levels are 4000 and 1000 m (see
   Appendix A for details).

-  For surface or shelter-level output, the <level> is not necessary.

Running UPP
===========

Scripts for running the UPP package are included in the tar file:

   *run\_unipost*

Before running any of the above listed scripts, perform the following instructions:

#. *cd* to your *DOMAINPATH* directory.

#. Make a directory to put the UPP results.

   *mkdir postprd*

#. Make a directory to put a copy of *postxconfig-NT.txt*.

   *mkdir parm*

#. Copy over the relevant control file to your working directory to
   customize *unipost*

   For WRF, copy the *UPPV4.0/parm/postxconfig-NT-WRF.txt* or for
   FV3, copy the *UPPV4.0/parm/postxconfig-NT-GFS.txt*.

   **Note:** *If you modified postcntrl.xml to reflect desired fields
   and levels, you will need to be sure that you generated the new flat file
   (please reference the section 'GRIB2 Control File' on how to do this)*

#. Copy over the (*UPPV4.0/scripts/run\_unipost*) script to **postprd/**.

#. Edit the run script as outlined below. Once these directories are set
   up and the edits outlined below are complete, the script can be run
   interactively from the **postprd** directory by simply typing the
   script name on the command line.

Overview of the scripts to run the UPP
--------------------------------------

**Note** *: It is recommended that the user refer to the
run\_unipost scripts in the script/ while reading this overview.*

User modified variables are all contained at the top of
the *run\_unipost* script in one user-edit section, along with a
brief description. Descriptions below follow the *run\_unipost*
script.

#. Set up basic path variables:

   **TOP\_DIR** : Top level directory for source code (*UPPV4.0*)

   **DOMAINPATH** : Working directory for this run

   **UNIPOST\_HOME** : Location of the *UPPV4.0* build directory

   **POSTEXEC** : Location of the *UPPV4.0* executables

   **modelDataPath** : Location of the model output data files to be
   processed (e.g. **wrfprd/** for WRF-based runs).

   **txtCntrlFile** : Name and location of *postxconfig-NT.txt*
   file that lists desired fields for GRIB2 format.
   This file is generated by the user following the steps listed above
   in the *'GRIB2 Control File'* section.

   **Note:** *For WRF, the scripts are configured such that
   unipost expects the WRF history files (wrfout
   files) to be in wrfprd/, the postxconfig-NT.txt file to be
   in parm/ and the postprocessor working directory to be
   called postprd/, all under DOMAINPATH*

   This set up is for user convenience to have a script ready to run,
   paths may be modified but be sure to check the run script to make
   sure settings are correct.

#. Specify dynamic core being run (ARW or FV3)

   **dyncore**: What model core is used (ARW or FV3)

#. Specify the format for the input model files and output UPP files.

   **inFormat**: Format of the model data

                 arw - "netcdf"
                 fv3 - "binarynemsiompiio"

   **outFormat**: Format of output from UPP

                  "grib2"

#. Specify the forecast cycles to be post-processed

   **startdate** : Forecast start date (YYYYMMDDHH)

   **fhr** : First forecast hour to be post-processed

   **lastfhr** : Last forecast hour to be post-processed

   **incrementhr** : Increment (in hours) between forecast files

                    \*Do not set to 0 or the script will loop continuously

#. Set up how many domains will be post-processed

   **domain\_list** : List of domains for run (e.g. d01 d02)

#. Set/uncomment the run command for your system. (i.e. serial, mpirun,
   etc).

   **RUN\_COMMAND** : System run command for serial or parallel runs

   -  The default execution command in the distributed scripts is for a
      single processor:

      *./unipost.exe > unipost\_${domain}.${fhr}.out 2>&1*

   -  To run unipost using mpi (dmpar compilation), the command line
      should be:

      >> LINUX-MPI systems: *mpirun -np N unipost.exe > outpost
      2>&1*

      (**Note:** *on some systems a host file also needs to be
      specified: -machinefile "host"*)

      >> IBM: *mpirun.lsf unipost.exe < itag > outpost*

      >> SGI MPT: *mpiexec\_mpt unipost.exe < itag > outpost*

#. Set naming convention for prefix and extension of output file name

   -  **comsp** is the initial string of the output file name (by default
      it is not setand the prefix of the output file will be the string set
      in the <datset> tag of the *postcntrl.xml*).

   -  **tmmark** is used for the file extension (in
      *run\_unipost*, *tmmark=tm00*, if not set, it is set to .GrbF)

The itag that will be read in by *unipost.exe* from
stdin (unit 5) is generated automatically in the *run\_unipost*
script based on the user-defined options above. It should not be
necessary to edit this. For description purposes, the namelist
(*itag*) contains 6 lines for FV3:

#. Name of the FV3 (pressure level) output file to be posted.

#. Format of FV3 model output (netcdf, binarynemsio).

#. Format of UPP output (GRIB2)

#. Forecast valid time (not model start time) in FV3 format (the
   forecast time desired to be post-processed).

#. Dynamic core used (GFS).

#. Name of the FV3 (surface) output file to be post-processed.

**Note:** *If the third line (i.e., UPP output type) is not set, UPP will default the output file format to "grib1".
UPP output for FV3 only supports GRIB2.*

If scripts *run\_unipostandgrads* or *run\_unipostandgempak* are
used, additional steps are taken to create image files (see
Visualization section below).

Upon a successful run, *unipost* will generate output files *GFSPRS.hh* in the postprocessor working directory, where *nn* refers to
the domain id and *hh* denotes the forecast hour. In addition, the script
*run\_unipostandgrads* will produce a suite of png images named
*variablehh\_GrADS.png*, and the script *run\_unipostandgempak*
will produce a suite of gif images named *variablehh.gif*.

If the run did not complete successfully, a log file in the
post-processor working directory called *unipost.hh.out*, where *nn* is the domain id and *hh* is the
forecast hour, may be consulted for further information.

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

**Note:** *wgrib2 is not distributed within the UPP tar
file. Users may download and install from
http://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/.*

Visualization with UPP
======================

GEMPAK 
-------

The GEMPAK utility *nagrib* is able to decode GRIB files whose
navigation is on any non-staggered grid. Hence, GEMPAK is able to decode
GRIB files generated by the UPP package and plot horizontal fields or
vertical cross sections.

A sample script named *run\_unipostandgempak*, which is included in
the scripts directory of the tar file, can be used to run *unipost*
and plot the following fields using GEMPAK:

-  **Sfcmap\_dnn\_hh.gif**: mean SLP and 6 hourly precipitation

-  **PrecipType\_dnn\_hh.gif**: precipitation type (just snow and
   rain)

-  **850mbRH\_dnn\_hh.gif**: 850 mb relative humidity

-  **850mbTempandWind\_dnn\_hh.gif**: 850 mb temperature and wind
   vectors

-  **500mbHandVort\_dnn\_hh.gif**: 500 mb geopotential height and
   vorticity

-  **250mbWindandH\_dnn\_hh.gif**: 250 mb wind speed isotacs and
   geopotential height

This script can be modified to customize fields for output. GEMPAK has
an online users guide at:

http://www.unidata.ucar.edu/software/gempak/help_and_documentation/manual/.

In order to use the script *run\_unipostandgempak*, it is necessary
to set the environment variable *GEMEXEC* to the path of the GEMPAK
executables. For example,

   *setenv GEMEXEC /usr/local/gempak/bin*

GrADS
-----

The GrADS utilities *g2ctl.pl* and *gribmap*
are able to decode GRIB2 files whose navigation is on any
non-staggered grid. These utilities and instructions on how to use them
to generate GrADS control files are available from:

http://www.cpc.ncep.noaa.gov/products/wesley/g2ctl.html (GRIB2).

The GrADS package is available from:
http://cola.gmu.edu/grads/gadoc/gadoc.php.

GrADS has an online Users Guide at:
http://cola.gmu.edu/grads/gadoc/users.html

A list of basic commands for GrADS can be found at:
http://cola.gmu.edu/grads/gadoc/reference_card.pdf.

A sample script named *run\_unipostandgrads*, which is included in
the scripts directory of the Unified Post Processing package, can be
used to run *unipost* and plot the following fields using GrADS:

-  **Sfcmaphh\_dnn\_GRADS.png**: mean SLP and 6-hour accumulated
   precipitation.

-  **850mbRHhh\_dnn\_GRADS.png**: 850 mb relative humidity

-  **850mbTempandWindhh\_dnn\_GRADS.png**: 850 mb temperature and wind
   vectors

-  **500mbHandVorthh\_dnn\_GRADS.png**: 500 mb geopotential heights
   and absolute vorticity

-  **250mbWindandHhh\_dnn\_GRADS.png**: 250 mb wind speed isotacs and
   geopotential heights

In order to use the script *run\_unipostandgrads*, it is necessary
to:

#. Set the environmental variable GADDIR to the path of the GrADS fonts
   and auxiliary files. For example:

   *setenv GADDIR /usr/local/grads/data*

#. Add the location of the GrADS executables to the PATH. For example:

   *setenv PATH /usr/local/grads/bin:$PATH*

#. Link script *cbar.gs* to the post-processor working directory.
   (This scripts is provided in UPP package, and the
   *run\_unipostandgrads* script makes a link from **scripts/** to
   **postprd/**.) To generate the plots above, GrADS script
   *cbar.gs* is invoked. This script can also be obtained from the
   GrADS library of scripts at
   http://cola.gmu.edu/grads/gadoc/library.html.

Fields produced by unipost
==========================

The 2 tables described below contain documentation regarding the fields
that are available for output by UPP for GRIB2.

Grib2 Table:
------------

https://dtcenter.org/upp/users/docs/tables/UPP_GRIB2_Table.pdf

This table lists basic and derived fields currently produced by unipost
for grib2. The abbreviated names listed in the second column of each
table describe how the fields should be entered in the
*post\_cntrl.xml*.

Appendix A: UPPV3.1+ Reflectivity field descriptions
====================================================

Reflectivities are filled/computed depending on the model core and
microphysics options.

UPP uses model derived reflectivity (REFL\_10CM from WRF) for model runs
using the Thompson microphysics option (mp=8). Other combinations use
algorithms within UPP code.

Work is underway to provide more user flexibility when selecting
reflectivity computations. For more information on model computed
reflectivity, e.g. REFL\_10CM, please see model documentation.

Relevant routines for reflectivity. Some or all of these may need to be
modified if the user desires to change where/how reflectivity is
processed. It is recommended that the user have knowledge of the model
output and reflectivity computations before trying to modify the UPP
code. Email `upp-help@ucar.edu <mailto:>`__ for further questions.

**INITPOST**\* - Separate routines for each different model core (e.g.
ARW, FV3, etc.). Reads model fields, e.g. REFL\_10CM, REFD\_MAX

**MDLFLD.f** - Computes DBZ or fills DBZ arrays with model computed
Reflectivity. - Fills 3-D model level reflectivity array (UPP ID: 250) -
Fills 2-D composite reflectivity array (UPP ID: 252)

**MDL2AGL.f** - Interpolates relevant DBZ array to AGL reflectivity
(UPP ID: 253) - Outputs model computed maximum hourly reflectivity
(REFD\_MAX; UPP ID: 421)

**MDL2P.f** - Interpolates relevant DBZ array to pressure levels (UPP
ID: 251)


Acknowledgement

If significant help was provided via the UPP helpdesk for work resulting
in a publication, please acknowledge the Developmental Testbed Center
UPP Team.

For referencing this document please use:

UPP Users Guide V3.0, 34 pp. [available online at
http://www.dtcenter.org/upp/users/docs/user\_guide/V3/upp\_users\_guide.pdf]
