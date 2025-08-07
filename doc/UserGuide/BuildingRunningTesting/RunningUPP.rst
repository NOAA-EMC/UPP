.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

.. _running-upp:

***********************
Running UPP Standalone
***********************

This section describes how to run UPP in a standalone environment using the ``run_upp`` script. It assumes that users already have forecast model output available and that they have already :ref:`cloned and built the UPP <building-upp>`. For an example of how to generate forecast model output from the UFS WM, see :numref:`Section %s <prepare-forecast>` of the UPP tutorial, or refer to the :doc:`authoritative UFS WM documentation <ufs-wm:index>` and `UFS WM GitHub Discussions Q & A <https://github.com/ufs-community/ufs-weather-model/discussions/categories/q-a>`_. 

#. Optional: If desired, edit the control XML file(s) in ``/UPP/parm`` to reflect the fields
   and levels you want UPP to output. For an example of how to do this, see :numref:`Section %s <modify-xml>`. It is recommended that you make copies of the original beforehand.

   * **GFS** text files are located in ``UPP/parm/gfs`` and include ``postxconfig-NT-gfs-f00-two.txt`` (0-hour lead time) and ``postxconfig-NT-gfs-two.txt`` (all other lead times).
   * **LAM** text files are located in ``UPP/parm/rrfs`` and include ``postxconfig-NT-rrfs.txt``. 

   Remake the flat text file(s) following the steps in :numref:`Section %s: Control File: Creating the Flat Text File <create_txt_file>` or :numref:`Section %s <modify-xml>`. 

#. Navigate to your top-level working directory. This directory will be referred to as ``${TOP_DIR}`` throughout the documentation. 

#. Create a directory to place the UPP results in: 

   .. code-block:: console

      cd ${TOP_DIR}
      mkdir postprd

   The UPP typically assumes that this directory will be named ``postprd``. 

#. Download the UPP utility script for running standalone UPP, and change the permissions:

   .. code-block:: console

      cd postprd
      wget https://raw.githubusercontent.com/wiki/NOAA-EMC/UPP/run_upp
      chmod 755 run_upp

#. Edit the run script as outlined below in :numref:`Section %s: Run Script Overview <run-script-overview>`. Once these directories are set
   up and the edits outlined below are complete, the script can be run from the ``postprd`` directory: 

   .. code-block:: console

      ./run_upp -c <compiler> -m <machine>

   .. note::
      The UPP is supported on Ursa, Orion, and Hercules NOAA RDHPCS. It will likely run on other machines, but users may have to create a modulefile for it and/or modify the ``run_upp`` script if loading the ``upp_common`` module doesn't work.

.. _run-script-overview:

===================
Run Script Overview
===================

It is recommended that the user refer to the ``run_upp`` script while reading this overview. All user-modified variables are contained at the top of the ``run_upp`` script in the user-edit section, along with a brief description. Descriptions below follow the ``run_upp`` script.

#. Set up basic path variables:

   * ``TOP_DIR``: Top level directory for building and running UPP
   * ``POSTPRD_DIR``: Working directory for this run
   * ``UPP_HOME``: Location of the **UPP** directory
   * ``POSTEXEC``: Location of the **UPP** executable
   * ``modelDataPath``: Location of the model output files to be processed by the UPP
   * ``txtCntrlFile``: Full path to the flat text file that lists desired fields for output. 
      
      * **GFS** text files are located in ``UPP/parm/gfs`` and include ``postxconfig-NT-gfs-f00-two.txt`` (0-hour lead time) and ``postxconfig-NT-gfs-two.txt`` (all other lead times).
      * **LAM** text files are located in ``UPP/parm/rrfs`` and include ``postxconfig-NT-rrfs.txt``. 

   * ``CRTMDIR``: Path to simulated synthetic satellite files (required for the LAM implementation)

   .. note::
      For FV3, the scripts are configured such that UPP expects the flat text file to be in ``/parm``,
      and the postprocessor working directory to be called ``/postprd``.
      This setup is for user convenience to have a script ready to run; paths may be modified, but be
      sure to check the run script to make sure settings are correct.

#. Specify dynamical core being run:

   * ``model``: Which model is used? Valid options: Global Forecast System (``GFS``) or Limited Area Model (``LAM``). 

   .. note::
      Forecast output should include ``atmf*``, ``sfcf*``, and ``GFSPRS*`` files when running with model=GFS, and it should include ``phyf*``, ``dynf*``, ``NATLEV*``, and ``PRSLEV*`` files when model=LAM. 

#. Specify the format for the input model files and output UPP files:

   * ``inFormat``: Format of the model data ("netcdfpara")
   * ``outFormat``: Format of output from UPP ("grib2")

   .. COMMENT: Maybe we shouldn't list this, since it's set, and the user can't/shouldn't change it? ^

#. Specify the forecast cycles to be post-processed:

   * ``startdate``: Forecast start date (YYYYMMDDHH)
   * ``fhr``: First forecast hour to be post-processed
   * ``lastfhr``: Last forecast hour to be post-processed
   * ``incrementhr``: Increment (in hours) between forecast files (cannot be set to zero)

#. Set/uncomment the run command for your system:

   * ``RUN_COMMAND``: System run commands

       |     - The default execution command in the distributed scripts is for a single processor:
       |       ``./upp.x > upp.${fhr}.out 2>&1``

       |     - To run UPP using :term:`MPI` (dmpar compilation), the command line should be:
       |       >> NOAA RDHPCS with Slurm-based job scheduler: ``srun -A <account> -n 4 ${POSTEXEC}/upp.x``
       |          (Note: ``<account>`` should be replaced with the actual name of an account where the user can charge 
       |          computational resources.)
       |       >> LINUX-MPI systems: ``mpirun -np N upp.x > outpost 2>&1``
       |          (Note: On some systems a host file also needs to be specified:
                  ``-machinefile "host"``)
       |       >> IBM: ``mpirun.lsf upp.x < itag > outpost``
       |       >> SGI MPT: ``mpiexec_mpt upp.x < itag > outpost``

#. Set the value for ``numx``.

   * ``numx``: The number of subdomains in the x-direction used for decomposition.

       |     - For 1D decomposition, set numx=1 (default)
       |     - For 2D decomposition, set numx>1

Upon a successful run, UPP will generate output files for each forecast hour in ``${POSTPRD_DIR}``.

When executed with the provided run script, UPP provides log files in the post-processor working directory named
``upp.fHHH.out``, where ``HHH`` is the forecast hour. These log files may be consulted for further runtime information in the event of an error.
