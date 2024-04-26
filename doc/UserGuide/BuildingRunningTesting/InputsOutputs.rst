.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

.. _input-output:

***********************
UPP Inputs and Outputs
***********************

This section describes the input files used to run the UPP and the resulting output files.

.. _input-files:

===========
Input Files
===========

The UPP requires the following input files:
 - The model forecast file
 - The ``itag`` namelist file
 - The :term:`GRIB2` control file (e.g., ``postxconfig-NT.txt``)
 - Additional data files (e.g., lookup tables, coefficient files for satellites)

.. _model-forecast:

--------------
Model Forecast
--------------

The UPP ingests FV3 :term:`write component` files in parallel :term:`netCDF` format.

The table below is a list of the unified model variables available from the :term:`FV3` model core. Whether a specific variable is able to be read by the UPP relies on dependencies such as physics options and choice of model. This table does not include variables that are diagnosed when running the UPP.

UFS Unified Model Variables
 - :doc:`../tables/UFS_unified_variables_table`

.. _itag:

----
ITAG
----

The file called ``itag`` is a Fortran namelist file that contains two sections: ``&model_inputs`` describing the as well as the 2D decomposition specification in ``&nampgb``. It is read in by the ``upp.x`` executable from standard input (stdin -- unit 5) and is generated automatically within the UFS application workflow or standalone run script based on user-defined options. It should not be necessary to edit this file. For description purposes, the namelists ``&model_inputs`` and ``&nampgb`` (in the ``itag`` file) contain the following lines for FV3:

:bolditalic:`&model_inputs`

model_inputs
^^^^^^^^^^^^

The ``&model_inputs`` section parameterizes choices about the set of model output files that will be used for the UPP 

.. list-table:: *Description of the &model_inputs namelist section.*
   :widths: 25 40 15 10
   :header-rows: 1

   * - Variable Name
     - Description
     - Data Type
     - Default Value
   * - ``datestring``
     - Time stamp being processed (e.g., 2022-08-02_19:00:00).
     - character(len=19)
     - ---
   * - ``fileName``
     - Name of input dynamics file; name of full 3-D model output file.
     - character(len=256)
     - ---
   * - ``filenameflat``
     - Input configuration text file defining the requested fields.
     - character(len=256)
     - postxconfig-NT.txt
   * - ``filenameflux``
     - Name of input physics file; name of 2-D model output file with physics and surface fields.
     - character(len=256)
     - ---
   * -``grib``
     - Grib type (Note that UPP only supports Grib2 currently)
     - character(5)
     - grib2
   * - ``ioform``
     - Input file format. Choices: binary, binarympiio, binarynemsio, netcdf, netcdfpara
     - character(len=20)
     - ---
   * - ``modelname``
     - Model name used by UPP internally (e.g., FV3R for LAM, GFS for GFS, NCAR for WRF)
     - character(len=4)
     - ---

.. list-table:: *Description of the ``&nampgb`` namelist section.*
   :widths: 25 40 15 10
   :header-rows: 1

   * - ``aqf_on bool``
     - Turn on Air Quality Forecasting (CMAQ-based)
     - logical
     - .false.
   * - ``d2d_chem``
     - Turn on option to process the 2D aerosol/chemical tracers
     - logical
     - .false.
   * - ``d3d_on``
     - Turn on option to use dynamic 3D fields from GFS
     - logical
     - .false.
   * - ``fileNameAERi`` string 
     - aerosols file name
     - character(len=256)
     - ''''
   * - ``gccpp_on``
     - Turn on option to process the aerosol/chemical tracers related output from UFS-Chem (CCPP-Chem) model
     - logical
     - .false.
   * - ``gocart_on``
     - Turn on option to process the aerosol/chemical tracers related output from GEFS-Aerosols model (GOCART)
     - logical
     - .false.
   * - ``gtg_on``
     - Turn on GTG (Graphical Turbulence Guidance)
     - logical
     - .false.
   * - ``hyb_sigp``
     - Not used
     - logical
     - .true.
   * - ``kpo``
     - The number of pressure levels, if different than standard one specified by ``SPLDEF`` described below.
     - integer
     - 0
   * - ``kpv``
     - The number of output potential vorticity levels
     - integer
     - 8
   * - ``kth``
     - The number of output isentropic levels
     - integer
     - 6
   * - ``method_blsn``
     - Turn on blowing snow effect on visibility diagnostic (default=true)
     - logical
     - .true.
   * - ``nasa_on``
     - Turn on option to process the aerosol/chemical tracers related output from UFS-Aerosols model (NASA GOCART)
     - logical
     - .false.
   * - ``numx``
     - The number of i regions in a 2D decomposition; Each i row is distibuted to numx ranks; numx=1 is the special case of a 1D decomposition in Y only.
     - integer
     - 1
   * - ``po``
     - List indicating pressure levels in output
     - real,dimension(70)
     - 0
   * - ``popascal``
     - Switch to indicate if pressure levels are in pascals (multply by 100 if false)
     - logical
     - .false.
   * - ``pv``
     - List indicating the potential vorticity level output
     - real,dimension(70)
     - (/0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0,(0.,k=kpv+1,70)/)
   * - ``rdaod``
     - Turn on the option to process the AOD from the GFS scheme
     - logical
     - .false.
   * - ``slrutah_on``
     - Calculate snow to liquid ratio (SLR) using method from University of Utah.(default=false)
     - logical
     - .false.
   * - ``th``
     - List indicating isentropic level output
     - real,dimension(70)
     - (/310.,320.,350.,450.,550.,650.,(0.,k=kth+1,70)/)
   * - ``vtimeunits string 
     - valid time units, default="", Choices: FMIN
     - character(len=4)
     - ''''
   * - ``write_ifi_debug_files``
     - Write debug files for In-Flight Icing (IFI), a restricted option in UPP
     - logical
     - .false.

.. _control-file:

------------
Control File
------------

The user interacts with the UPP through the control file to define what fields and levels to output. It is composed of a header and a body. The header specifies the output file information. The body includes which fields and levels to process.

A default control file, ``postxconfig-NT.txt``, is provided and read by the UPP. Users who wish to customize the control file to add or remove fields and/or levels may do so by modifying ``postcntrl.xml`` and then remaking the text file as described in the later section: :ref:`create_txt_file`.

.. Note::
   The control file names ``postxconfig-NT.txt`` and ``postcntrl.xml`` are generic names and are different depending on the application used. Control files for various operational models are located in the ``UPP/parm`` directory.

.. _control-output:

Selecting Which Variables the UPP Outputs
-------------------------------------------

To output a field, the body of the control file needs to contain an entry for the appropriate variable. If an entry for a particular field is not yet available in the control file, it may be added to the control file with the appropriate entries for that field. For variables found on vertical levels (e.g., isobaric or height levels), the desired levels to be output must be listed (see next section: :ref:`control_levels`). A list of available GRIB2 fields that can be output by UPP can be found in the table :doc:`../tables/UPP_GRIB2_Table_byID`. Please note that some fields are dependent on model, physics, and other fields.

.. _control_levels:

Controlling which levels the UPP outputs
----------------------------------------

The ``<level>`` tag in the ``postcntrl.xml`` file is used to list the desired levels for output. The following levels are currently available for output:

- For isobaric output, 46 levels are possible, from 2 to 1000 hPa (*2, 5, 7, 10, 20, 30, 50, 70 mb and then every 25 mb from 75 to 1000 mb*). The complete list of levels is specified in ``sorc/ncep_post.fd/CTLBLK.f``.

   - Modify specification of variable ``LSMDEF`` to change the number of pressure levels: LSMDEF=47
   - Modify specification of ``SPLDEF`` array to change the values of pressure levels:
     (/200.,500.,700.,1000.,2000.,3000.,5000.,7000.,7500.,10000.,12500.,15000.,17500.,20000., …/)

- For model-level output, all model levels are possible, from the highest to the lowest.
- When using the Noah LSM, the soil layers are 0-10 cm, 10-40 cm, 40-100 cm, and 100-200 cm.
- When using the RUC LSM, the soil levels are 0 cm, 1 cm, 4 cm, 10 cm, 30 cm, 60 cm, 100 cm, 160 cm, and 300 cm. (For the old RUC LSM, there are only 6 layers, and if using this, you will need to change ``NSOIL`` for “RUC LSM” from 9 to 6 in the ``sorc/ncep_post.fd/WRFPOST.f`` routine.)
- When using Pliem-Xiu LSM, there are two layers: 0-1 cm, 1-100 cm
- For low, mid, and high cloud layers, the layers are :math:`\geq`\ 642 hPa, :math:`\geq`\ 350 hPa, and <350 hPa, respectively.
- For PBL layer averages, the levels correspond to 6 layers with a thickness of 30 hPa each.
- For flight level, the levels are 30 m, 50 m, 80 m, 100 m, 305 m, 457 m, 610 m, 914 m, 1524 m, 1829 m, 2134 m, 2743 m, 3658 m, 4572 m, 6000 m, 7010 m.
- For AGL radar reflectivity, the levels are 4000 and 1000 m.
- For surface or shelter-level output, the ``<level>`` is not necessary.

.. _create_txt_file:

Creating the Flat Text File
---------------------------

If the control file requires any modifications, a preprocessing step will be required by the user to convert the modified XML file ``parm/postcntrl.xml`` to a flat text file ``parm/postxconfig-NT.txt``. The user will first need to edit the ``postcntrl.xml`` file to declare which fields are to be output from the UPP.

In order to ensure that the user-edited XML files are error free, XML stylesheets (``parm/EMC_POST_CTRL_Schema.xsd`` and ``EMC_POST_Avblflds_Schema.xsd``) can be used to validate both the ``postcntrl.xml`` and ``post_avblflds.xml`` files respectively. Confirmation of validation will be given (e.g., ``postcntrl.xml`` validates) or otherwise return errors if it does not match the schema. This step is optional, but acts as a safeguard to avoid run-time failures with the UPP. To run the validation:

.. code-block:: console

    xmllint --noout --schema EMC_POST_CTRL_Schema.xsd postcntrl.xml
    xmllint --noout --schema EMC_POST_Avblflds_Schema.xsd post_avblflds.xml

Once the XMLs are validated, the user will need to generate the flat file. The command below will run the Perl program ``parm/PostXMLPreprocessor.pl`` to generate the post flat file:

.. code-block:: console

    /usr/bin/perl PostXMLPreprocessor.pl your_user_defined_xml post_avblflds.xml your_user_defined_flat

where ``your_user_defined_xml`` is your modified XML and ``your_user_defined_flat`` is the output text file.

.. _output-files:

============
Output Files
============

Upon a successful run, ``upp.x`` will generate GRIB2 output files in the post processor working directory. These files will include all fields that were requested in the control file.

When running UPP standalone, the following GRIB2 output files will be generated:

   | **GFS Model**: ``GFSPRS.HHH``
   | **LAM (Limited Area Model)**: ``NATLEV.HHH`` and ``PRSLEV.HHH``

When executed with the provided run script, UPP provides log files in the post-processor working directory named ``upp.fHHH.out``, where ``HHH`` is the forecast hour. These log files may be consulted for further runtime information in the event of an error.
