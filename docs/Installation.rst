.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

********************
Building Stand-Alone
********************

The UPP uses a CMake-based build system to integrate all the required components for building the UPP.
Once built, the UPP can be run stand-alone (outside the UFS Applications) to post-process model output.

=====================
Software Requirements
=====================

The UPP is tested on a variety of research platforms, including NOAA HPC systems (e.g. Hera, Orion) and
the NCAR HPC Cheyenne. These supported platforms are pre-configured for building and running the UPP and already
have the required libraries available via `HPC-Stack<https://github.com/NOAA-EMC/hpc-stack>`_ in a centralized
location. The HPC-Stack is a script-based build system that builds the software stack required by UFS components.

Users working on unsupported platforms will need to install the HPC-Stack on their system and can do so following
the instructions in the `HPC-Stack User's Guide<https://hpc-stack.readthedocs.io/en/latest/>`_.

============================
Obtaining and Installing UPP
============================

Building and running UPP V10.0.12 has been tested and is supported on the following pre-configured platforms.

+---------------+----------------------+
| System        | Compiler and Version |
+===============+======================+
| NCAR Cheyenne | Intel 2021.2         |
|               +----------------------+
|               | GNU 10.1.0           |
+---------------+----------------------+
| NOAA Hera     | Intel 18.0.5.274     |
+---------------+----------------------+
| NOAA Orion    | Intel 2018.4         |
+---------------+----------------------+

Move to the directory where you want to install UPP and clone the repository.

.. code-block:: console

    git clone -b branch-or-tag-name https://github.com/NOAA-EMC/UPP

where, ``branch-or-tag-name`` is the release branch or tag you wish to clone.


.. code-block:: console

    cd UPP/tests

    ./compile_upp.sh

.. note::
   To build in debug mode, you can add :bolditalic:`-DCMAKE_BUILD_TYPE=Debug` to the *cmake_opts*
   parameter in the :bolditalic:`compile_upp.sh` script.
   This removes compiler optimization flags and adds -g to the fortran compilation. You can also use
   :bolditalic:`-DCMAKE_BUILD_TYPE=RELWITHDEBINFO`, which gives the -g, but keeps the -O2 optimization
   for the fortran compilation.

Move back to the top level UPP directory and create a directory for the CRTM fix files to be unpacked
in. Download the fix files from the Github `release page
<https://github.com/NOAA-EMC/UPP/releases/tag/upp_v10.0.12>`_ or use the wget command. Unpack the
tar file.

.. code-block:: console

    cd ../
    mkdir crtm && cd crtm
    wget https://github.com/NOAA-EMC/UPP/releases/download/upp_v10.0.12/fix.tar.gz
    tar -xzf fix.tar.gz

.. note::
   To make a clean build, simply remove both the **tests/build** and **tests/install** directories and the
   :bolditalic:`exec/upp.x` executable and then rerun the :bolditalic:`compile_upp.sh script. This is
   recommended if a mistake is made during the installation process.
   
=======================
UPP Directory Structure
=======================

Under the main directory **UPP** reside the following relevant subdirectories (The * indicates a
directory that exists only after the build is complete):

     | **exec***: Contains the :bolditalic:`upp.x` executable after successful compilation

     | **modulefiles**: Contains modulefiles for specific platforms and compilers for building on pre-
       configured machines.

     | **parm**: Contains parameter files, which can be modified by the user to control how the post
       processing is performed.

     | **scripts**: Contains a sample run script to process fv3 history files.
     |   - **run_upp**: runs :bolditalic:`upp.x`.

     | **sorc**: Contains source codes for:
     |   - **ncep_post.fd**: Source code for the UPP

     | **tests**: Contains the scripts used to install UPP
     |   - **build***: Contains the UPP build
     |   - **install***: Contains the installed executable, modules, and libraries
