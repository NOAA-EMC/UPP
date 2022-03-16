.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

******************
Inputs and Outputs
******************

This section describes the input files used when running the UPP and the resulting output files.

===========
Input files
===========

The UPP requires the following input files:
 - The itag namelist file
 - The GRIB2 control file
 - Additional data files (e.g. lookup tables, coefficient files for satellite)

----
ITAG
----

The file called :bolditalic:`itag` is a text file that contains the fortran namelist &model_inputs. It is
read in by the :bolditalic:`upp.x` executable from stdin (unit 5) and is
generated automatically within the UFS application workflow or stand-alone run script based on
user-defined options. It should not be necessary to edit this. For description purposes, the namelist
&model_inputs (:bolditalic:`itag` file) contains 7 lines for FV3:

#. Name of the FV3 (pressure level) output file to be posted.

#. Format of FV3 model output (netcdf, binarynemsiompiio).

#. Format of UPP output (GRIB2)

#. Forecast valid time (not model start time) in YYYY-MM-DD_HH:00:00 format (the forecast time desired
   to be post-processed).

#. Model used (GFS, FV3R - regional FV3 also known as the LAM - Limited Area Model).

#. Name of the FV3 (surface) output file to be post-processed.

#. Name of configuration file (postxconfig-NT.txt)

------------
Control File
------------

The user interacts with the UPP through the control file to define what fields and levels to output. It
is composed of a header and a body. The header specifies the output file information. The body includes
which fields and levels to process.

A default control file, :bolditalic:`postxconfig-NT.txt`, is provided and read by the UPP. For users
wishing to customize the control file to add or remove fields and/or levels, they may do so by
modifying the :bolditalic:`postcntrl.xml` and then remaking the text file as described in the later section
:ref:`Creating the Flat Text File`.

.. Note::
   The control file names :bolditalic:`postxconfig-NT.txt` and :bolditalic:`postcntrl.xml` are generic
   names and are different depending on the application used.

The tables below list all fields that are included in the control files for the various UFS
applications. All fields in the tables may not be present in your output depending on whether the field
dependencies are available in your model output.

UFS MRW Table (GFS model)
 - :doc:`MRW_GFSPRS_table`

UFS SRW Tables (LAM - Limited Area Model)
 - :doc:`SRW_BGDAWP_table`
 - :doc:`SRW_BGRD3D_table`

Controlling which variables the UPP outputs
-------------------------------------------

To output a field, the body of the control file needs to contain an entry for the appropriate variable.
If an entry for a particular field is not yet available in the control file, it may be added to the
control file with the appropriate entries for that field. For variables found on vertical levels (e.g.
isobaric or height levels), the desired levels to be output must be listed (see next section:
:ref:`control_levels`). A list of available Grib2 fields that can be output by UPP can be found in the 
table :doc:`UPP_GRIB2_Table`. Please note that some fields are dependent on model, physics, and other fields.

.. _control_levels:

Controlling which levels the UPP outputs
----------------------------------------

The <level> tag in the postcntrl.xml file is used to list the desired levels for output. The following
levels are currently available for output:

- For isobaric output, 46 levels are possible, from 2 to 1000 hPa (*2, 5, 7, 10, 20, 30, 50, 70 mb and
  then every 25 mb from 75 to 1000 mb*). The complete list of levels is specified in
  :bolditalic:`sorc/ncep_post.fd/CTLBLK.f`.
  
   - Modify specification of variable LSMDEF to change the number of pressure levels: LSMDEF=47
   - Modify specification of SPLDEF array to change the values of pressure levels:
     (/200.,500.,700.,1000.,2000.,3000.,5000.,7000.,7500.,10000.,12500.,15000.,17500.,20000., …/)
      
- For model-level output, all model levels are possible, from the highest to the lowest.
- When using the Noah LSM, the soil layers are 0-10 cm, 10-40 cm, 40-100 cm, and 100-200 cm.
- When using the RUC LSM, the soil levels are 0 cm, 1 cm, 4 cm, 10 cm, 30 cm, 60 cm, 100 cm, 160 cm,
  and 300 cm. (For the old RUC LSM, there are only 6 layers and if using this, you will need to change
  “RUC LSM” from 9 to 6 in the :bolditalic:`sorc/ncep_post.fd/WRFPOST.f` routine.)
- When using Pliem-Xiu LSM, there are two layers: 0-1 cm, 1-100 cm
- For low, mid, and high cloud layers, the layers are :math:`\geq`\ 642 hPa, :math:`\geq`\ 350 hPa, and
  <350 hPa, respectively.
- For PBL layer averages, the levels correspond to 6 layers with a thickness of 30 hPa each.
- For flight level, the levels are 30 m, 50 m, 80 m, 100 m, 305 m, 457 m, 610 m, 914 m, 1524 m, 1829 m,
  2134 m, 2743 m, 3658 m, 4572 m, 6000 m, 7010 m.
- For AGL radar reflectivity, the levels are 4000 and 1000 m (see Appendix A for details).
- For surface or shelter-level output, the <level> is not necessary.

Creating the Flat Text File
---------------------------

If the control file requires any modifications, a preprocessing step will be required by the user to
convert the modified xml file :bolditalic:`parm/postcntrl.xml` to a flat text file
:bolditalic:`parm/postxconfig-NT.txt`. The user will first need to edit the :bolditalic:`postcntrl.xml`
file to declare which fields are to be output from the UPP.

In order to ensure that the user-edited xml files are error free, XML stylesheets
(:bolditalic:`parm/EMC\_POST\_CTRL\_Schema.xsd` and :bolditalic:`EMC\_POST\_Avblflds\_Schema.xsd`) can
be used to validate both the :bolditalic:`postcntrl.xml` and :bolditalic:`post\_avblflds.xml` files,
respectively. Confirmation of validation will be given (e.g. postcntrl.xml validates) or otherwise
return errors if it does not match the schema. This step is optional, but acts as a safe-guard to avoid
run-time failures with UPP. To run the validation:

.. code-block:: console

    xmllint --noout --schema EMC_POST_CTRL_Schema.xsd postcntrl.xml
    xmllint --noout --schema EMC_POST_Avblflds_Schema.xsd post_avblflds.xml

Once the xmls are validated, the user will need to generate the flat file. The makefile will call the
perl program :bolditalic:`parm/POSTXMLPreprocessor.pl` to regenerate any post flat files
:bolditalic:`postxconfig-NT.txt` where modifications were made since it was last run. Generate the flat
file:

.. code-block:: console

    make

============
Output Files
============

Upon a successful run, :bolditalic:`upp.x` will generate GRIB2 output files in the post-processor
working directory. These files will include all fields that were requested in the control file.

When running UPP stand-alone, the following Grib2 output files will be generated:

   | **GFS Model**: GFSPRS.HHH
   | **LAM (Limited Area Model)**: NATLEV.HHH and PRSLEV.HHH

When executed with the provided run script, UPP provides log files in the post-processor working directory named
:bolditalic:`upp.fHHH.out`, where :bolditalic:`HHH` is the forecast hour. These log files may be consulted for further
run-time information in the event of an error.
