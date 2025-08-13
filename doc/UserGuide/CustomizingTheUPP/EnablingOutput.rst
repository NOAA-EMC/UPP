.. role:: bolditalic
    :class: bolditalic

.. _enabling-output:

***************************************************
Tutorial: Enabling User-Defined Sigma Level Output
***************************************************

This documentation describes the steps required to generate UPP output for temperature, U, and V components at user-defined sigma levels. These correspond to indices 206, 208, and 209 in the :doc:`GRIB2 <../tables/UPP_GRIB2_Table_byID>` table. The instructions address the configuration of XML files and the resolution of potential issues encountered during this process.

Create a Working Directory
===========================

Create and navigate into a top-level working directory. This directory will be referred to as ``${TOP_DIR}`` throughout the tutorial.

.. code-block:: console

   mkdir wrk_dir
   cd wrk_dir

.. note:: 

   Users can optionally save their working directory in an environment variable. For example: 

   .. code-block:: console

      export TOP_DIR=/path/to/wrk_dir

.. _prepare-forecast:

Prepare Forecast Output
========================

#. The UPP requires forecast output as its input. In this example, we will use output created from a UFS Weather Model (WM) forecast run. This subsection is provided for ease of use, but users should refer to the :doc:`authoritative UFS WM documentation <ufs-wm:index>` and `UFS WM GitHub Discussions Q & A <https://github.com/ufs-community/ufs-weather-model/discussions/categories/q-a>`_ for questions:

   .. code-block:: console

      git clone --recursive https://github.com/ufs-community/ufs-weather-model.git
      cd ufs-weather-model/tests

#. Create a custom configuration file with the test or tests you want to run. 

   .. code-block:: console

      vi mytests.conf
      i

   .. note::
      Users can check `rt.conf <https://github.com/ufs-community/ufs-weather-model/blob/develop/tests/rt.conf>`_ for a full list of tests and select the ones they want to run. The UPP expects ``atmf*``, ``sfcf*``, and ``GFSPRS*`` files for the GFS model, and it expects ``phyf*``, ``dynf*``, ``NATLEV*``, and ``PRSLEV*`` files for the LAM model, so users should run and/or adapt tests with these expected outputs.

#. Paste the following lines, which compile and run two sample cases (``control_p8`` and ``regional_control``), into ``mytests.conf``:

   .. code-block:: console

      COMPILE | atm_dyn32 | intel | -DAPP=ATM -DCCPP_SUITES=FV3_GFS_v16,FV3_GFS_v16_flake,FV3_GFS_v17_p8,FV3_GFS_v17_p8_rrtmgp,FV3_GFS_v15_thompson_mynn_lam3km,FV3_WoFS_v0,FV3_GFS_v17_p8_mynn,FV3_GFS_v17_p8_ugwpv1,FV3_GFS_v16_gfdlmpv3,FV3_GFS_v17_p8_ugwpv1_tempo -D32BIT=ON | | fv3 |
      RUN | control_p8                                        | - noaacloud                          | baseline |
      RUN | regional_control                                  |                                      | baseline |

#. Save and exit: Type ``esc`` and ``:wq``. Then hit ``Enter``/``Return``.

#. Edit each test's output frequency by modifying the test files under ``ufs-weather-model/tests/tests``. For example, to output every three hours for 24 hours:

   .. code-block:: console

      export OUTPUT_FH='3 -1'

   .. note:: 

      ``'3 -1'`` means every 3 hours. The ``-1`` indicates that it is not a list of hours. For example, someone could choose to output an arbitrary list of hours like ``export OUTPUT_FH='0 4 5 7 12 16 19'``. 

6. If you are not running on NOAA :term:`RDHPCS` or do not have access to the ``stmp*`` disk space on those systems, alter the ``dprefix`` path for the machine you are running on in ``rt.sh`` to point to a directory (such as ``$TOP_DIR``) where you have write permissions. For example:

   .. code-block:: console

      dprefix=${dprefix:-"path/to/${TOP_DIR}"}

   .. note::

      For users not running on NOAA :term:`RDHPCS`, see the UFS WM documentation for :ref:`building the Weather Model on other systems <ufs-wm:build-wm>`.

#. Run the forecast tests:

   .. code-block:: console

      nohup ./rt.sh -e -k -a epic -l mytests.conf &
   
   See the UFS WM documentation for ``rt.sh`` :ref:`optional arguments <ufs-wm:cmd-line-opts>`. 


UPP Procedures
================

.. _clone:

Clone the UPP
--------------

Clone the ``UPP`` into ``${TOP_DIR}``:

.. code-block:: console

   cd ${TOP_DIR}
   git clone https://github.com/NOAA-EMC/UPP.git
   cd UPP/parm

.. _modify-xml:

Modify ``postcntrl*.xml``
--------------------------

An XML :ref:`control file <control-file>` determines what fields and levels the UPP will output. Control files for various operational models are located in the ``UPP/parm`` directory. The ``post_avblflds.xml`` file contains all fields that the UPP can currently output.
   
Select a ``postcntrl*.xml`` file that is most relevant to your experiment. In this example, for a GFS experiment, navigate to ``UPP/parm/gfs`` and choose an XML from that directory to modify. 

If the user wishes to generate model output at user-defined sigma levels for temperature, U, and V values on sigma surfaces, copy the entries from ``post_avblflds.xml``
(206, 208, and 209), and add them to, e.g., ``postcntrl_gfs_f00_two.xml``. Then, add the desired levels to the entries using a ``<level></level>`` tag. Users may choose to remove extraneous information, including ``post_avblfldidx`` and ``pname``. For example:

.. code-block:: console

   <param>
      <post_avblfldidx>206</post_avblfldidx>
      <shortname>TMP_ON_SIGMA_LVLS</shortname>
      <pname>TMP</pname>
      <fixed_sfc1_type>sigma_lvl</fixed_sfc1_type>
      <level>0. 500. 1000. 1500. 2000. 2500. 3000. 3500. 4000. 4500. 5000. 5500. 6000. 6500. 7000. 7500. 8000. 8500. 9000. 9500. 10000.</level>
      <scale>4.0</scale>
   </param>

   <param>
      <post_avblfldidx>208</post_avblfldidx>
      <shortname>UGRD_ON_SIGMA_LVLS</shortname>
      <pname>UGRD</pname>
      <fixed_sfc1_type>sigma_lvl</fixed_sfc1_type>
      <level>0. 500. 1000. 1500. 2000. 2500. 3000. 3500. 4000. 4500. 5000. 5500. 6000. 6500. 7000. 7500. 8000. 8500. 9000. 9500. 10000.</level>
      <scale>4.0</scale>
   </param>

   <param>
      <post_avblfldidx>209</post_avblfldidx>
      <shortname>VGRD_ON_SIGMA_LVLS</shortname>
      <pname>VGRD</pname>
      <fixed_sfc1_type>sigma_lvl</fixed_sfc1_type>
      <level>0. 500. 1000. 1500. 2000. 2500. 3000. 3500. 4000. 4500. 5000. 5500. 6000. 6500. 7000. 7500. 8000. 8500. 9000. 9500. 10000.</level>
      <scale>4.0</scale>
   </param>

.. _validate-xml:

Validate the XML postcontrol file
-----------------------------------

After modifying the XML postcontrol file, users must validate the file to ensure it was formatted successfully:

.. code-block:: console

   cd UPP/parm
   # For GFS/Global
   xmllint --noout --schema EMC_POST_CTRL_Schema.xsd gfs/postcntrl_gfs_f00_two.xml
   # For LAM/regional
   xmllint --noout --schema EMC_POST_CTRL_Schema.xsd rrfs/rrfs_postcntrl.xml

If the XML file is formatted correctly, it will print a message like:

.. code-block:: console

   gfs/postcntrl_gfs_f00_two.xml validates

.. _gen-text:

Generate the flat text file
----------------------------

If the control file requires any modifications, the user must convert the modified XML file to a flat text file. The command below will run the Perl program ``parm/PostXMLPreprocessor.pl`` to generate the ``postxconfig*`` flat file:

   .. code-block:: console

      cd UPP/parm
      # Global/GFS
      /usr/bin/perl PostXMLPreprocessor.pl gfs/postcntrl_gfs_f00_two.xml post_avblflds.xml gfs/postxconfig-NT-gfs-f00-two.txt
      # OR Regional/RRFS
      /usr/bin/perl PostXMLPreprocessor.pl rrfs/rrfs_postcntrl.xml post_avblflds.xml rrfs/postxconfig-NT-rrfs.txt

.. include:: ../doc-snippets/run-parm.rst  

.. include:: ../doc-snippets/ursa-workaround.rst

.. _modify-sorc:

Modify UPP source code
------------------------

By default, only certain sigma levels are outputted. These levels are defined in `SET_LVLSXML.f <https://github.com/NOAA-EMC/UPP/blob/develop/sorc/ncep_post.fd/SET_LVLSXML.f>`_ using the ``ASIGO1`` array. Users must review these levels to confirm compatibility with their requirements. If default sigma levels are insufficient, users must modify ``SET_LVLSXML.f`` to include their desired sigma levels.

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

.. _build-upp:



Build/Compile the UPP
----------------------

After modifying the source code, users must (re)build/(re)compile the UPP to ensure that changes take effect: 

.. code-block:: console

   cd UPP/tests
   ./compile_upp.sh

This will generate the UPP executable in the ``UPP/exec`` directory. 

.. _postprd:

Prepare the post-processing working directory
----------------------------------------------

#. Create a post-processing output directory:

   .. code-block:: console

      cd ${TOP_DIR}
      mkdir postprd

   .. note::
      This directory can be created anywhere, but default settings assume that it is named ``postprd`` and created inside ``${TOP_DIR}``.

#. Download the UPP utility script for running standalone UPP and change the permissions:

   .. code-block:: console

      cd postprd
      wget https://raw.githubusercontent.com/wiki/NOAA-EMC/UPP/run_upp
      chmod 755 run_upp

#. Modify the script:

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
      # export RUN_COMMAND="${POSTEXEC}/upp.x "

      #MPI sample command
      # "-n 4" can be changed to a different number of tasks. 
      export RUN_COMMAND="srun -A epic -n 4 ${POSTEXEC}/upp.x "

      # The number of subdomains in the x-direction (set to >=2 for 2d decomposition)
      export numx=1
      ...

   .. include:: ../doc-snippets/expected-files.rst

#. For regional (LAM) post-processing only:
   
   Regional/LAM post-processing requires satellite files, which can be downloaded from the UPP repository:

   .. code-block:: console

      mkdir crtm && cd crtm
      wget https://github.com/NOAA-EMC/UPP/releases/download/upp_v11.0.0/fix.tar.gz
      tar -xzf fix.tar.gz

   Make sure to adjust the ``CRTMDIR`` path in ``run_upp`` to point to the location of these files.

.. _run-runupp:

Run UPP
--------

#. Run the script:

   .. include:: ../doc-snippets/run-upp.rst

   Users can add the ``-v`` option to see more output.
   Check the ``upp.f{fhr}.out`` files to see if there were any errors when the script ran. 

#. Checking output:
   
   To view information about a file, use the *wgrib2* utility. Users may need to load the *wgrib2* module first (e.g., via ``module load wgrib2``). For example, to see summary information from ``GFSPRS.000``, run: 

   .. code-block:: console

      wgrib2 -s GFSPRS.000

   For more detailed information, use the ``-V`` flag instead, but note that this is computationally intensive and may require allocating and using a compute node on some systems.

   To see information about a specific variable, run:

   .. code-block:: console

      wgrib2 GFSPRS.000 -match ":0.995 sigma level:"
   
   To obtain statistics for a specific variable, use the ``-stats`` flag. To print the values of the variables, use the ``-text`` flag along with the name of a file where the values can be printed:

   .. code-block:: console

      wgrib2 GFSPRS.000 -match ":UGRD:0.995 sigma level:" -stats
      wgrib2 GFSPRS.000 -match ":UGRD:0.995 sigma level:" -text gfsprs.000.txt
