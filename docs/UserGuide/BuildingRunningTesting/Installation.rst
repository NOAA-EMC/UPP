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

The UPP is tested on a variety of research platforms, including NOAA HPC systems (e.g., Hera, Orion). These supported platforms are pre-configured for building and running the UPP and already
have the required libraries available via `spack-stack <https://github.com/JCSDA/spack-stack>`__ in a centralized
location. The :term:`spack-stack` is a :term:`Spack`-based method for installing UFS prerequisite software libraries.

Users working on unsupported platforms will need to install spack-stack on their system and can do so following
the instructions in the :ref:`spack-stack User's Guide <spack-stack:index>`.

.. note::

   Users can install :term:`HPC-Stack` instead of spack-spack by following the instructions in the :ref:`HPC-Stack User's Guide <hpc-stack:index>`. However, support for HPC-Stack is being deprecated, and limited assistance is available for use of HPC-Stack with the UPP. 

----------------
Common Modules
----------------

As of February 1, 2024, the UPP uses the following `common modules <https://github.com/NOAA-EMC/UPP/blob/develop/modulefiles/upp_common.lua>`__ from spack-stack: 

.. code-block:: console

   cmake 3.16.1+
   hdf5/1.14.0
   netcdf-c 4.9.2
   netcdf-fortran 4.6.1
   jasper 2.0.32
   libpng 1.6.37 / png 1.6.35
   zlib 1.2.13
   g2 3.4.5
   g2tmpl 1.10.2
   parallelio 2.5.10
   bacio 2.4.1
   ip 4.3.0
   sp 2.5.0
   crtm 2.4.0.1
   w3emc 2.10.0
   nemsio 2.5.4
   sigio 2.3.2
   sfcio 1.4.1
   wrf_io 1.2.0

Individual machines may subsequently load slightly different versions. The most updated list of modules for a given machine can be viewed `in each machine's modulefile <https://github.com/NOAA-EMC/UPP/tree/develop/modulefiles>`__. 
Users on a non-Tier-1 system should look at the modulefile for the system 
whose architecture most closely resembles their own system's architecture.

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

To install the UPP, navigate to the directory where you want to install UPP and clone the repository.

.. code-block:: console

    git clone -b branch-or-tag-name https://github.com/NOAA-EMC/UPP

where, ``branch-or-tag-name`` is the release branch or tag you wish to clone (e.g., ``upp_v11.0.0``). (Leaving off the ``-b`` argument will clone all branches of the repository.)

Move to the directory with the build script and build the UPP.

.. code-block:: console

    cd UPP/tests

    ./compile_upp.sh

.. note::
   To build in debug mode, you can add ``-DCMAKE_BUILD_TYPE=Debug`` to the *cmake_opts* parameter in the :bolditalic:`compile_upp.sh` script.
   This removes compiler optimization flags and adds ``-g`` to the fortran compilation. You can also use
   ``-DCMAKE_BUILD_TYPE=RELWITHDEBINFO``, which gives the ``-g``, but keeps the ``-O2`` optimization
   for the fortran compilation.

Move back to the top level UPP directory and create a directory where the CRTM fix files will be unpacked. Download the fix files from the GitHub `release page
<https://github.com/NOAA-EMC/UPP/releases/tag/upp_v11.0.0>`__ or use the ``wget`` command. Unpack the tar file.

.. code-block:: console

    cd ../
    mkdir crtm && cd crtm
    wget https://github.com/NOAA-EMC/UPP/releases/download/upp_v11.0.0/fix.tar.gz
    tar -xzf fix.tar.gz

.. note::
   To make a clean build, simply remove both the **tests/build** and **tests/install** directories and the
   :bolditalic:`exec/upp.x` executable and then rerun the :bolditalic:`compile_upp.sh` script. This is
   recommended if a mistake is made during the installation process.
   
