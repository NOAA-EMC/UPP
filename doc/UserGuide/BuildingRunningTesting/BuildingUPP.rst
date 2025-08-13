.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

.. _building-upp:

*************************
Building UPP Standalone
*************************

The UPP uses a CMake-based build system to integrate all the required components for building the UPP.
Once built, the UPP can be run standalone (outside the :term:`UFS` Applications) to post-process model output.

=====================
Software Requirements
=====================

The UPP is tested on a variety of research platforms, including NOAA HPC systems (e.g., Ursa, Orion). These supported platforms are preconfigured for building and running the UPP and already
have the required libraries available via `spack-stack <https://github.com/JCSDA/spack-stack>`__ in a centralized
location. The :term:`spack-stack` is a :term:`Spack`-based method for installing UFS prerequisite software libraries.

Users working on unsupported platforms will need to install spack-stack on their system and can do so following
the instructions in the :doc:`spack-stack User's Guide <spack-stack:index>`.

----------------
Common Modules
----------------

As of August 8, 2025, the UPP uses the following `common modules <https://github.com/NOAA-EMC/UPP/blob/develop/modulefiles/upp_common.lua>`_ from spack-stack: 

.. code-block:: console

   hdf5 1.14.3
   netcdf-c 4.9.2
   netcdf-fortran 4.6.1
   jasper 2.0.32
   libpng 1.6.37
   zlib 1.2.13
   g2 3.5.1
   g2tmpl 1.13.0
   bacio 2.4.1
   ip 5.1.0
   crtm 2.4.0.1
   w3emc 2.10.0
   nemsio 2.5.4
   sigio 2.3.3
   wrf-io 1.2.0

Individual machines may subsequently load slightly different versions. The most updated list of modules for a given machine can be viewed `in each machine's modulefile <https://github.com/NOAA-EMC/UPP/tree/develop/modulefiles>`__. 
Users on non-Tier-1 systems should look at the modulefile for the system 
whose architecture most closely resembles their own system's architecture to determine which modules they may need.

============================
Obtaining and Installing UPP
============================

Building and running UPP v11.0.0 has been tested and is supported on the following pre-configured platforms.

+---------------+----------------------+
| System        | Compiler and Version |
+===============+======================+
| NOAA Hera     | Intel 18.0.5.274     |
+---------------+----------------------+
| NOAA Orion    | Intel 2018.4         |
+---------------+----------------------+

To install the UPP, navigate to the directory where you want to install UPP and clone the repository. This directory will be referred to as ``${TOP_DIR}`` throughout the documentation.

.. include:: ../doc-snippets/clone.rst

where, ``branch-or-tag-name`` is the release branch or tag you wish to clone (e.g., ``upp_v11.0.0``). (Leaving off the ``-b`` argument will clone the default ``develop`` branch of the repository.)

Move to the directory with the build script and build the UPP.

.. code-block:: console

   cd UPP/tests
   ./compile_upp.sh

.. note::
   To build in debug mode, you can activate ``-DCMAKE_BUILD_TYPE=Debug`` by adding the ``-d`` flag when running the ``compile_upp.sh`` script.
   This removes compiler optimization flags and adds the ``-g`` flag to the Fortran compilation. You can also use
   ``-DCMAKE_BUILD_TYPE=RELWITHDEBINFO``, which adds the ``-g`` flag and uses the ``-O2`` optimization
   for the Fortran compilation.

---------------------------------------
Post-processing CRTM Files (optional)
---------------------------------------

Move back to the top-level UPP directory and create a directory where the CRTM fix files will be unpacked. Download the fix files from the GitHub `release page
<https://github.com/NOAA-EMC/UPP/releases/tag/upp_v11.0.0>`_ or use the ``wget`` command. Unpack the tar file.

.. code-block:: console

    cd ../
    mkdir crtm && cd crtm
    wget https://github.com/NOAA-EMC/UPP/releases/download/upp_v11.0.0/fix.tar.gz
    tar -xzf fix.tar.gz

.. hint::
   To make a clean build, simply remove both the ``tests/build`` and ``tests/install`` directories and the
   ``exec/upp.x`` executable and then rerun the ``compile_upp.sh`` script. This is recommended if a mistake is made during the installation process.
   
