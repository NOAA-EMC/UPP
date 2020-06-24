******************
Inputs and Outputs
******************

This section describes the input files used when running the UPP and the resulting output files

===========
Input files
===========

The UPP requires the following input files:
 - The itag namelist
 - The GRIB2 control file
 - Additional data files (e.g. lookup tables, coefficient files for satellite)

----
ITAG
----

The *itag* namelist that is read in by *unipost.exe* from
stdin (unit 5) is generated automatically in the *run\_unipost*
script based on user-defined options. It should not be
necessary to edit this. For description purposes, the namelist
(*itag*) contains 7 lines for FV3:

#. Name of the FV3 (pressure level) output file to be posted.

#. Format of FV3 model output (netcdf, binarynemsio).

#. Format of UPP output (GRIB2)

#. Forecast valid time (not model start time) in FV3 format (the
   forecast time desired to be post-processed).

#. Dynamic core used (GFS).

#. Name of the FV3 (surface) output file to be post-processed.

#. Name of configuration file (postxconfig-NT.txt)

------------
Control File
------------

The user interacts with unipost through the control file to define what fields and levels to output. It is composed of a header and a body. The header specifies the output file information. The body includes which fields and levels to process.

A default control file, *postxconfig-NT.txt*, is provided and read by unipost. For users wishing to customize the control file to add or remove fields and/or levels, they may do so by modyfying the postcntrl.xml and then remaking the text file required by unipost.

The `GRIB2 Output Table <https://dtcenter.org/sites/default/files/community-code/upp-grib2-table\_0.pdf>`_ lists basic and derived fields currently produced by unipost
for grib2.

Controlling which variables unipost outputs
-------------------------------------------

To output a field, the body of the control file needs to contain an
entry for the appropriate variable. For variables found on isobaric or
height levels, the desired levels to be output must be listed (see next
section: *Controlling which levels unipost outputs*). If an entry for
a particular field is not yet available in the control file, it  may be
added to the control file with the appropriate entries for that field.

Controlling which levels unipost outputs
----------------------------------------

The <level> tag in the postcntrl.xml is used to list the desired levels
for output. The following levels are currently available for output:

- For isobaric output, 46 levels are possible, from 2 to 1000 hPa (*2,
  5, 7, 10, 20, 30, 50, 70 mb and then every 25 mb from 75 to 1000
  mb*). The complete list of levels is specified in *src/unipost/CTLBLK.f*.
  
   - Modify specification of variable LSMDEF to change the number of
     pressure levels: LSMDEF=47
   - Modify specification of SPLDEF array to change the values of
     pressure levels: (/200.,500.,700.,1000.,2000.,3000.
     &,5000.,7000.,7500.,10000.,12500.,15000.,17500.,20000., …/)
      
- For model-level output, all model levels are possible, from the
  highest to the lowest.
- When using the Noah LSM, the soil layers are 0-10 cm, 10-40 cm,
  40-100 cm, and 100-200 cm.
- When using the RUC LSM, the soil levels are 0 cm, 1 cm, 4 cm, 10 cm,
  30 cm, 60 cm, 100 cm, 160 cm, and 300 cm. (For the old RUC LSM,
  there are only 6 layers and if using this, you will need to change
  “RUC LSM” from 9 to 6 in the WRFPOST.f routine.)
- When using Pliem-Xiu LSM, there are two layers: 0-1 cm, 1-100 cm
- For low, mid, and high cloud layers, the layers are :math:`\geq`\ 642
  hPa, :math:`\geq`\ 350 hPa, and <350 hPa, respectively.
- For PBL layer averages, the levels correspond to 6 layers with a
  thickness of 30 hPa each.
- For flight level, the levels are 30 m, 50 m, 80 m, 100 m, 305 m, 457
  m, 610 m, 914 m, 1524 m, 1829 m, 2134 m, 2743 m, 3658 m, 4572 m, 6000
  m, 7010 m.
- For AGL radar reflectivity, the levels are 4000 and 1000 m (see
  Appendix A for details).
- For surface or shelter-level output, the <level> is not necessary.

Creating the Flat Text File
---------------------------

For outputting GRIB2 format using version 4.0, a preprocessing step
is required by the user to convert the xml file
*parm/postcntrl.xml* to a flat text file
*parm/postxconfig-NT.txt*. The flat file is quicker to process
than the xml file. The user will first need to edit the
*postcntrl.xml* file to declare which fields are to be output
from UPP.

In order to ensure that the user-edited xml files are error free, XML
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

Once the xmls are validated, the user will need to generate the flat
file. Edit the *parm/makefile* if necessary to point to the
desired flat file directory and xmls. The makefile will call the perl
program *parm/POSTXMLPreprocessor.pl* to generate the post flat
file *postxconfig-NT.txt*. Generate the flat file:

  *make*

============
Output Files
============

Upon a successful run, *unipost* will generate GRIB2 output files
*GFSPRS.hh* in the postprocessor working directory, where *hh* denotes
the forecast hour. These files will include all fields that were
requested in the control file.

If the run did not complete successfully, a log file in the
post-processor working directory called *unipost.hh.out*, where *hh*
is the forecast hour, may be consulted for further information.
