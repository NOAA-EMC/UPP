.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

.. _running-upp:

***********************
Running UPP Stand-Alone
***********************

This section describes how to prepare model output, clone and build UPP, generate control files, and run UPP in a standalone environment using the provided ``run_upp`` script.

Create a Working directory
===========================

Create and enter a top-level working directory. This directory will be referred to as ``TOP_DIR``.

   .. code-block:: console

      mkdir wrk_dir
      cd wrk_dir

Prepare Forecast Output
========================

1. The UPP needs to use forecast output as input to the UPP. In this example we will use output from a UFS WM forecast run:

   .. code-block:: console

      git clone --recursive https://github.com/ufs-community/ufs-weather-model.git
      cd ufs-weather-model/tests

2. Create a custom configuration file with the test or tests you want to run. 

   .. code-block:: console

      vi mytests.configuration
      i

.. note::
   Users can check `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`_ for a full lists of tests and select the ones they want to run.

3. Paste the following lines, which compile and run two sample cases (``control_c48`` and ``control_p8``)

   .. code-block:: console

      COMPILE | atm_dyn32 | intel | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16,FV3_GFS_v16_flake,FV3_GFS_v17_p8,FV3_GFS_v17_p8_rrtmgp,FV3_GFS_v15_thompson_mynn_lam3km,FV3_WoFS_v0,FV3_GFS_v17_p8_mynn,FV3_GFS_v17_p8_ugwpv1,FV3_GFS_v16_gfdlmpv3,FV3_GFS_v17_p8_ugwpv1_tempo -D32BIT=ON | | fv3 |
       RUN     | control_c48               |                                      | baseline |
       RUN     | control_p8                | - noaacloud                          | baseline |

4. Save and exit:

   .. code-block:: console

      <ESC>
      :wq

5. Edit each test's output frequency by modifying the test files under ``ufs-weather-model/tests/tests``. For example, to output every three hours for 24 hours:

   .. code-block:: console

      export OUTPUT_FH= `3 -1`
 
6. If you do not have access to the ``stmp`` disk space, you made need to alter the ``dprefix`` path for the machine youre running on in ``rt.sh``. For example:

   .. code-block:: console

      dprefix=${dprefix:-"/work2/noaa/epic/gpetro/orion/test-sigma/stmp"}.


7. Run the forecasts tests:

   .. code-block:: console

      nohup ./rt.sh -e -k -a epic -l mytests.conf &

UPP Procedures
================

1. Clone the ``UPP`` into the ``wrk_dir``:

   .. code-block:: console

      cd wrk_dir
      git clone git clone https://github.com/NOAA-EMC/UPP.git
      cd UPP/parm

2. Modifying the postcntrl*.xml file.

Select a ``postcntrl*.xml`` file that is most relevant to your experiment. In this example, for a GFS experiment, navigate to to GFS and choose an XML from that directory to modify. 

If the user wishes to generate model output at user-defined sigma levels for temperature, U, and V values on sigma surfaces, copy the entries from ``post_avblflds.xml``
(206, 208, and 209) and add them to, e.g., ``postcntrl_gfs_f00_two.xml``. Then, add the desired levels to the entries using a <level></level> tag. Users may choose to remove extraneous information, including ``post_avblfldidx`` and ``pname``. For example:

   .. code-block:: console

      <param>
      <post_avblfldidx>206</post_avblfldidx>
      <shortname>TMP_ON_SIGMA_LVLS</shortname>
      <pname>TMP</pname>
      <fixed_sfc1_type>sigma_lvl</fixed_sfc1_type>
      <scale>4.0</scale>
      <level>0.0000 0.0500 0.1000 0.1500 0.2000 0.2500 0.3000 0.3500 0.4000 0.4500 0.5000 0.5500 0.6000 0.6500 0.7000 0.7500 0.8000 0.8500 0.9000 0.9500 1.0000</level>
   </param>

   <param>
      <post_avblfldidx>208</post_avblfldidx>
      <shortname>UGRD_ON_SIGMA_LVLS</shortname>
      <pname>UGRD</pname>
      <fixed_sfc1_type>sigma_lvl</fixed_sfc1_type>
      <scale>4.0</scale>
      <level>0.0000 0.0500 0.1000 0.1500 0.2000 0.2500 0.3000 0.3500 0.4000 0.4500 0.5000 0.5500 0.6000 0.6500 0.7000 0.7500 0.8000 0.8500 0.9000 0.9500 1.0000</level>
   </param>

   <param>
      <post_avblfldidx>209</post_avblfldidx>
      <shortname>VGRD_ON_SIGMA_LVLS</shortname>
      <pname>VGRD</pname>
      <fixed_sfc1_type>sigma_lvl</fixed_sfc1_type>
      <scale>4.0</scale>
      <level>0.0000 0.0500 0.1000 0.1500 0.2000 0.2500 0.3000 0.3500 0.4000 0.4500 0.5000 0.5500 0.6000 0.6500 0.7000 0.7500 0.8000 0.8500 0.9000 0.9500 1.0000</level>
   </param>

3. Generate the flat text file:

   .. code-block:: console

      cd UPP/parm
      /usr/bin/perl PostXMLPreprocessor.pl gfs/postcntrl_gfs_f00_two.xml
      post_avblflds.xml gfs/postxconfig-NT-gfs-f00-two.txt

.. note::
   ``PostXMLPreprocessor.pl`` must be run from the ``parm`` directory or it will produce an error. 




A script (``run_upp``) for running the UPP package is now fetched via ``wget``

:underline:`Before running the script, perform the following instructions:`

1. Create a working directory. This directory will be reffered to as ``TOP_DIR``.



2. Make a directory to put the UPP results in.

   .. code-block:: console

       mkdir postprd

3. Make a directory for staging a copy of the desired control file.

   .. code-block:: console

       mkdir parm

4. Optional: If desired, edit the control XML file(s) in ``/UPP/parm`` to reflect the fields
   and levels you want UPP to output. It is recommended that you make copies of the original
   beforehand.

   | **GFS XMLs**: ``postcntrl_gfs_f00.xml`` (0-hour lead time) and
     ``postcntrl_gfs.xml`` (all other lead times)
   | **LAM (Limited Area Model) XML**: ``fv3lam.xml``

   Remake the flat text file(s) following the steps in the "Control File: Creating the Flat Text File"
   section.

5. Copy the flat text file(s) to the ``/parm`` directory in your ``DOMAINPATH``. These are the files
   that UPP reads directly.

   | **GFS text files**: ``postxconfig-NT-GFS-F00.txt`` (0-hour lead time) and
     ``postxconfig-NT-GFS.txt`` (all other lead times).
   | **LAM text file**: ``postxconfig-NT-fv3lam.txt``

6. Navigate to the ``/postprd`` directory and retrieve the ``./run_upp`` script via ``wget``:

   .. code-block:: console

      cd /postprd
      wget https://raw.githubusercontent.com/wiki/NOAA-EMC/UPP/run_upp
      chmod 755 run_upp


7. Edit the run script as outlined in the :ref:`"Run Script Overview" <run-script-overview>` section below. Once these directories are set
   up and the edits outlined below are complete, the script can be run interactively from the
   ``/postprd`` directory by simply typing the script name on the command line.

.. _run-script-overview:

===================
Run Script Overview
===================

.. note::
   It is recommended that the user refer to the ``run_upp`` script while reading this overview. All user-modified variables are contained at the top of the ``run_upp`` script in the user-edit section, along with a brief description. Descriptions below follow the ``run_upp`` script.

1. Set up basic path variables:

   * ``TOP_DIR``: Top level directory for building and running UPP
   * ``DOMAINPATH``: Working directory for this run
   * ``UPP_HOME``: Location of the **UPP** directory
   * ``POSTEXEC``: Location of the **UPP** executable
   * ``modelDataPath``: Location of the model output data files to be processed by the UPP
   * ``txtCntrlFile``: Name and location of the flat text file that lists desired fields for output.

   .. note::
      For FV3, the scripts are configured such that UPP expects the flat text file to be in ``/parm``,
      and the postprocessor working directory to be called ``/postprd``, all under ``DOMAINPATH``.
      This setup is for user convenience to have a script ready to run; paths may be modified, but be
      sure to check the run script to make sure settings are correct.

2. Specify dynamical core being run:

   * ``model``: Which model is used? ("GFS" or "LAM" - Limited Area Model)

3. Specify the format for the input model files and output UPP files:

   * ``inFormat``: Format of the model data ("netcdfpara")
   * ``outFormat``: Format of output from UPP ("grib2")

4. Specify the forecast cycles to be post-processed:

   * ``startdate``: Forecast start date (YYYYMMDDHH)
   * ``fhr``: First forecast hour to be post-processed
   * ``lastfhr``: Last forecast hour to be post-processed
   * ``incrementhr``: Increment (in hours) between forecast files
       
   .. attention::
         
      Do not set ``incrementhr`` to 0 or the script will loop continuously! 

5. Set/uncomment the run command for your system (e.g., ``mpirun``).

   * ``RUN_COMMAND``: System run commands

       |     - The default execution command in the distributed scripts is for a single processor:
       |       ``./upp.x > upp.${fhr}.out 2>&1``

       |     - To run UPP using :term:`MPI` (dmpar compilation), the command line should be:
       |       >> LINUX-MPI systems: ``mpirun -np N upp.x > outpost 2>&1``
       |          (Note: On some systems a host file also needs to be specified:
                  ``-machinefile "host"``)
       |       >> IBM: ``mpirun.lsf upp.x < itag > outpost``
       |       >> SGI MPT: ``mpiexec_mpt upp.x < itag > outpost``

6. Set the value for ``numx``.

   * ``numx``: The number of subdomains in the x-direction used for decomposition.

       |     - For 1D decomposition, set numx=1 (default)
       |     - For 2D decomposition, set numx>1

7. Set naming convention for prefix and extension of output file name.
   
   * ``comsp`` is the initial string of the output file name. By default, it is not set, and the prefix of the output file will be the string set in the ``postcntrl.xml`` file ``DATSET`` parameter. If set, it will concatenate the setting to the front of the string specified in the XML file ``DATSET`` parameter.
   * ``tmmark`` is used for the file extension (in ``run_upp``, ``tmmark=tm00``; if not set, it is set to ``.GrbF``)

Upon a successful run, UPP will generate output files for each forecast hour in the ``/postprd`` directory.

When executed with the provided run script, UPP provides log files in the post-processor working directory named
``upp.fHHH.out``, where ``HHH`` is the forecast hour. These log files may be consulted for further runtime information in the event of an error.
