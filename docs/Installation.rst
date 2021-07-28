.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

********************
Building Stand-Alone
********************

=====================
Software Requirements
=====================

Before installing the UPP code, it is necessary to ensure that you have the required libraries
available on your system. These libraries include:

  - The external NCEP libraries
    https://github.com/NOAA-EMC/NCEPLIBS-external

  - The NCEP libraries
    https://github.com/NOAA-EMC/NCEPLIBS

An introduction of each can be found in their respective top level :bolditalic:`README.md` files.
Detailed instructions for building the libraries on various platforms can be found in the
**NCEPLIBS-external/doc** directory.

Certain machines do have the NCEP libraries in a pre-installed location for use to build UPP. Paths to
these pre-installed libraries are available on the
`UFS-SRW wiki <https://github.com/ufs-community/ufs-srweather-app/wiki/Supported-Platforms-and-Compilers>`_
and include platform name and compiler version.

============================
Obtaining and Installing UPP
============================

Building and running UPP V9.0.0 has been tested on the following platforms using pre-configured libraries.

+---------------+----------------------+
| System        | Compiler and Version |
+===============+======================+
| NCAR Cheyenne | Intel 19.1.1         |
|               +----------------------+
|               | GNU 9.1.0            |
|               | GNU 10.1.0           |
+---------------+----------------------+
| NOAA Hera     | Intel 18.0.5.274     |
+---------------+----------------------+

Move to the directory where you want to clone and build UPP and clone the repository into the directory
EMC_post.

.. code-block:: console

    git clone -b release-tag-name --recurse-submodules https://github.com/NOAA-EMC/EMC_post

where, ``release-tag-name`` is the release tag you wish to clone (e.g. for stand-alone UPP version 9, use
the release tag :bolditalic:`upp-v9.0.0`).

Move into the top level UPP directory and create and move into the build directory. Then build the UPP code
using the cmake utility.
The path ``INSTALL_PREFIX`` should point to the location of the pre-installed NCEP libraries.

.. code-block:: console

    cd EMC_post
    mkdir build && cd build

    cmake .. -DCMAKE_INSTALL_PREFIX=.. -DCMAKE_PREFIX_PATH=${INSTALL_PREFIX}
    make install

.. note::
   To build in debug mode, you can add :bolditalic:`-DCMAKE_BUILD_TYPE=Debug` to the cmake command.
   This removes compiler optimization flags and adds -g to the fortran compilation. You can also use
   :bolditalic:`-DCMAKE_BUILD_TYPE=RELWITHDEBINFO`, which gives the -g, but keeps the -O2 optimization
   for the fortran compilation.

Move back to the top level UPP directory and create a directory for the CRTM fix files to be unpacked
in. Download the fix files from the Github `release page
<https://github.com/NOAA-EMC/EMC_post/releases/tag/upp-v9.0.0>`_ or use the wget command. Unpack the
tar file.

.. code-block:: console

    cd ../
    mkdir crtm && cd crtm
    wget https://github.com/NOAA-EMC/EMC_post/releases/download/upp-v9.0.0/fix.tar.gz
    tar -xzf fix.tar.gz

.. note::
   To make a clean build, simply remove both the **/build** directory and the
   :bolditalic:`bin/upp.x` executable and then re-create the build from step #2. This is recommended if a
   mistake is made during the installation process. If a simple change is made to the code, you can simply
   type :bolditalic:`make install` again in the pre-existing build directory.
   
=======================
UPP Directory Structure
=======================

Under the main directory **EMC_post** reside the following relevant subdirectories (The * indicates a
directory that exists only after the build is complete):

     | **bin***: Contains the :bolditalic:`upp.x` executable after successful compilation

     | **build**: Contains the UPP build

     | **include***: Contains include modules built/used during compilation of UPP

     | **lib***: Libraries built/used by UPP that are separate from NCEPlibs

     | **parm**: Contains parameter files, which can be modified by the user to control how the post
       processing is performed.

     | **scripts**: Contains a sample run script to process fv3 history files.
     |   - **run_upp**: runs :bolditalic:`upp.x`.

     | **sorc**: Contains source codes for:
     |   - **ncep_post.fd**: Source code for the UPP
