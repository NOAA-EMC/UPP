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

Create and navigate into a top-level working directory. This directory will be referred to as ``TOP_DIR`` throughout the documentation.

   .. code-block:: console

      mkdir wrk_dir
      cd wrk_dir

Prepare Forecast Output
========================

#. The UPP needs to use forecast output as input to the UPP. In this example we will use output created from a UFS WM forecast run:

   .. code-block:: console

      git clone --recursive https://github.com/ufs-community/ufs-weather-model.git
      cd ufs-weather-model/tests

#. Create a custom configuration file with the test or tests you want to run. 

   .. code-block:: console

      vi mytests.conf
      i

.. note::
   Users can check `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`_ for a full lists of tests and select the ones they want to run.

.. COMMENT: Add info about files UPP expects? 

#. Paste the following lines, which compile and run two sample cases (``control_p8`` and ``regional_control``), into ``mytests.conf``:

   .. code-block:: console

      COMPILE | atm_dyn32 | intel | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16,FV3_GFS_v16_flake,FV3_GFS_v17_p8,FV3_GFS_v17_p8_rrtmgp,FV3_GFS_v15_thompson_mynn_lam3km,FV3_WoFS_v0,FV3_GFS_v17_p8_mynn,FV3_GFS_v17_p8_ugwpv1,FV3_GFS_v16_gfdlmpv3,FV3_GFS_v17_p8_ugwpv1_tempo -D32BIT=ON | | fv3 |
      RUN | control_p8                                        | - noaacloud                          | baseline |
      RUN | regional_control                                  |                                      | baseline |


#. Save and exit:

   .. code-block:: console

      <ESC>
      :wq

#. Edit each test's output frequency by modifying the test files under ``ufs-weather-model/tests/tests``. For example, to output every three hours for 24 hours:

   .. code-block:: console

      export OUTPUT_FH= `3 -1`
 
#. If you do not have access to the ``stmp`` disk space, you made need to alter the ``dprefix`` path for the machine youre running on in ``rt.sh``. For example:

   .. code-block:: console

      dprefix=${dprefix:-"/work2/noaa/epic/gpetro/orion/test-sigma/stmp"}.


#. Run the forecasts tests:

   .. code-block:: console

      nohup ./rt.sh -e -k -a epic -l mytests.conf &

UPP Procedures
================

1. Clone the ``UPP`` into the ``wrk_dir``:

   .. code-block:: console

      cd wrk_dir
      git clone git clone https://github.com/NOAA-EMC/UPP.git
      cd UPP/parm

2. Modifying the postcntrl*.xml file (optional):

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

4. Modify source code

By default, only certain sigma levels are outputted. These levels are defined in `SET_LVLSXML.f <https://github.com/NOAA-EMC/UPP/blob/develop/sorc/ncep_post.fd/SET_LVLSXML.f>`_ using the ASIGO1 array. Users must review these levels to confirm comatibility with their requirements. If deafult sigma levels are insufficient users must modify ``SET_LVLSXML.f`` to include your desired sigma levels.

For example in ``UPP/sorc/ncep_post.fd/SET_LVLSXML.f`` set:

   .. code-block:: console

      ELSE  ! SPECIFY SIGO
         ASIGO1( 1)=   0.0000
         ASIGO1( 2)=   0.0500
         ASIGO1( 3)=   0.1000
         ASIGO1( 4)=   0.1500
         ASIGO1( 5)=   0.2000
         ASIGO1( 6)=   0.2500
         ASIGO1( 7)=   0.3000
         ASIGO1( 8)=   0.3500
         ASIGO1( 9)=   0.4000
         ASIGO1(10)=   0.4500
         ASIGO1(11)=   0.5000
         ASIGO1(12)=   0.5500
         ASIGO1(13)=   0.6000
         ASIGO1(14)=   0.6500
         ASIGO1(15)=   0.7000
         ASIGO1(16)=   0.7500
         ASIGO1(17)=   0.8000
         ASIGO1(18)=   0.8500
         ASIGO1(19)=   0.9000
         ASIGO1(20)=   0.9500
         ASIGO1(21)=   1.0000

5. Build/Compile the UPP:

   .. code-block:: console

      cd UPP/tests
      ./compile_upp.sh

This will generate the UPP executable in the ``UPP/exec`` directory

#. Create a post-processing output directory:

   .. code-block:: console

      cd $TOP_DIR
      mkdir postprd

.. note::
   This directory can be created anywhere but default settings assume that is named postprd and created inside $TOP_DIR

#. Download the UPP utility script for running standalone UPP and change the permissions:

   .. code-block:: console

      cd postprd
      wget https://raw.githubusercontent.com/wiki/NOAA-EMC/UPP/run_upp
      chmod 755 run_upp

#. Modifying the script

Users will need to edit directory paths and start date for the experiment. It may also be necessary to modify the run command, model type, and I/O file formats. For example:

   .. code-block:: console

      ...
      export TOP_DIR=/work2/noaa/epic/jsmith/hercules/test-sigma/
      export POSTPRD_DIR=${TOP_DIR}/postprd
      export UPP_HOME=${TOP_DIR}/UPP
      export POSTEXEC=${UPP_HOME}/exec
      export modelDataPath=${TOP_DIR}/control_p8_intel
      export txtCntrlFile=${UPP_HOME}/parm/gfs/postxconfig-NT-gfs-f00-two.txt
      export CRTMDIR=${UPP_HOME}/crtm/fix
      # Set date/time information
      export startdate=2021032206
      export fhr=00
      export lastfhr=06
      export incrementhr=03
      # Specify model ("GFS" or "LAM" in upper case)
      export model="GFS"
      # Set input format from model and ouput format from UPP
      export inFormat="netcdfpara"
      export outFormat="grib2"
      # Set run command: 
      # Single processor command example
      export RUN_COMMAND="${POSTEXEC}/upp.x "

      #MPI sample command
      # "-n 4" can be changed to a different number of tasks. 
      export RUN_COMMAND="srun -A epic -n 4 ${POSTEXEC}/upp.x "

      # The number of subdomains in the x-direction (set to >=2 for 2d decomposition)
      export numx=1
      ...

.. _run-script-overview:

===================
Run Script Overview
===================

.. note::
   It is recommended that the user refer to the ``run_upp`` script while reading this overview. All user-modified variables are contained at the top of the ``run_upp`` script in the user-edit section, along with a brief description. Descriptions below follow the ``run_upp`` script.

#. Set up basic path variables:

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

#. Specify dynamical core being run:

   * ``model``: Which model is used? ("GFS" or "LAM" - Limited Area Model)

#. Specify the format for the input model files and output UPP files:

   * ``inFormat``: Format of the model data ("netcdfpara")
   * ``outFormat``: Format of output from UPP ("grib2")

#. Specify the forecast cycles to be post-processed:

   * ``startdate``: Forecast start date (YYYYMMDDHH)
   * ``fhr``: First forecast hour to be post-processed
   * ``lastfhr``: Last forecast hour to be post-processed
   * ``incrementhr``: Increment (in hours) between forecast files
       
   .. attention::
         
      Do not set ``incrementhr`` to 0 or the script will loop continuously! 

#. Set/uncomment the run command for your system (e.g., ``mpirun``).

   * ``RUN_COMMAND``: System run commands

       |     - The default execution command in the distributed scripts is for a single processor:
       |       ``./upp.x > upp.${fhr}.out 2>&1``

       |     - To run UPP using :term:`MPI` (dmpar compilation), the command line should be:
       |       >> LINUX-MPI systems: ``mpirun -np N upp.x > outpost 2>&1``
       |          (Note: On some systems a host file also needs to be specified:
                  ``-machinefile "host"``)
       |       >> IBM: ``mpirun.lsf upp.x < itag > outpost``
       |       >> SGI MPT: ``mpiexec_mpt upp.x < itag > outpost``

#. Set the value for ``numx``.

   * ``numx``: The number of subdomains in the x-direction used for decomposition.

       |     - For 1D decomposition, set numx=1 (default)
       |     - For 2D decomposition, set numx>1

#. Set naming convention for prefix and extension of output file name.
   
   * ``comsp`` is the initial string of the output file name. By default, it is not set, and the prefix of the output file will be the string set in the ``postcntrl.xml`` file ``DATSET`` parameter. If set, it will concatenate the setting to the front of the string specified in the XML file ``DATSET`` parameter.
   * ``tmmark`` is used for the file extension (in ``run_upp``, ``tmmark=tm00``; if not set, it is set to ``.GrbF``)

Upon a successful run, UPP will generate output files for each forecast hour in the ``/postprd`` directory.

When executed with the provided run script, UPP provides log files in the post-processor working directory named
``upp.fHHH.out``, where ``HHH`` is the forecast hour. These log files may be consulted for further runtime information in the event of an error.
