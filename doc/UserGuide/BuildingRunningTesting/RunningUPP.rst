.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

.. _running-upp:

***********************
Running UPP Stand-Alone
***********************

A script (``run_upp``) for running the UPP package is included in the ``/scripts`` directory. 

:underline:`Before running the script, perform the following instructions:`

1. ``cd`` to your ``DOMAINPATH`` directory. This is the top working directory for the forecast run.

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

6. Copy the ``/scripts/run_upp`` script to the ``/postprd`` directory.

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
