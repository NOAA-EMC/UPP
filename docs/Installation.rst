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

The UPP is tested on a variety of research platforms, including NOAA HPC systems (e.g., Hera, Orion) and
the NCAR HPC Cheyenne. These supported platforms are pre-configured for building and running the UPP and already
have the required libraries available via `HPC-Stack <https://github.com/NOAA-EMC/hpc-stack>`__ in a centralized
location. The :term:`HPC-Stack` is a script-based build system that builds the software stack required by UFS components.

Users working on unsupported platforms will need to install the HPC-Stack on their system and can do so following
the instructions in the `HPC-Stack User's Guide <https://hpc-stack.readthedocs.io/en/latest/>`__.

.. note::

   UFS applications are gradually shifting to :term:`spack-stack`, which is a :term:`Spack`-based method for installing the same UFS prerequisite software libraries installed by HPC-Stack. The spack-stack is currently used on NOAA Cloud platforms and in containerized UFS applications, while HPC-Stack is still used on other Level 1 systems and is the software stack validated by the UFS Weather Model. Users are encouraged to check out `spack-stack <https://github.com/NOAA-EMC/spack-stack>`__ to prepare for the upcoming shift in support from HPC-Stack to spack-stack. Users can install spack-spack instead of HPC-Stack by following the instructions in the `spack-stack User's Guide <https://spack-stack.readthedocs.io/en/latest/>`.

----------------
Common Modules
----------------

As of June 14, 2023, the UPP uses the following common modules from HPC-Stack: 

.. code-block:: console

   cmake 3.16.1+
   hdf5/1.10.6
   netcdf 4.7.4
   jasper 2.0.22+
   libpng 1.6.37 / png 1.6.35
   zlib 1.2.11
   g2 3.4.1+
   g2tmpl 1.10.0+
   bacio 2.4.1
   ip 3.3.3
   sp 2.3.3
   crtm 2.3.0
   w3emc 2.9.2
   nemsio 2.5.2+
   sigio 2.3.2
   sfcio 1.4.1
   wrf_io 1.1.1+


The most updated list of modules can be viewed in each machine's modulefile 
`here <https://github.com/NOAA-EMC/UPP/tree/develop/modulefiles>`__. 
Users on a non-Tier-1 system should look at the modulefile for the system 
whose architecture most closely resembles their own system's architecture.

============================
Obtaining and Installing UPP
============================

Building and running UPP v11.0.0 has been tested and is supported on the following pre-configured platforms.

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

To install the UPP, navigate to the directory where you want to install UPP and clone the repository.

.. code-block:: console

    git clone -b branch-or-tag-name https://github.com/NOAA-EMC/UPP

where, ``branch-or-tag-name`` is the release branch or tag you wish to clone. (Leaving this off the ``-b`` argument will clone all branches of the repository.)

Move to the directory with the build script and build the UPP.

.. code-block:: console

    cd UPP/tests

    ./compile_upp.sh

.. note::
   To build in debug mode, you can add ``-DCMAKE_BUILD_TYPE=Debug`` to the *cmake_opts*
   parameter in the :bolditalic:`compile_upp.sh` script.
   This removes compiler optimization flags and adds ``-g`` to the fortran compilation. You can also use
   ``-DCMAKE_BUILD_TYPE=RELWITHDEBINFO``, which gives the ``-g``, but keeps the ``-O2`` optimization
   for the fortran compilation.

Move back to the top level UPP directory and create a directory where the CRTM fix files will be unpacked. Download the fix files from the GitHub `release page
<https://github.com/NOAA-EMC/UPP/releases/tag/upp_v11.0.0>`__ or use the ``wget`` command. Unpack the
tar file.

.. code-block:: console

    cd ../
    mkdir crtm && cd crtm
    wget https://github.com/NOAA-EMC/UPP/releases/download/upp_v11.0.0/fix.tar.gz
    tar -xzf fix.tar.gz

.. note::
   To make a clean build, simply remove both the **tests/build** and **tests/install** directories and the
   :bolditalic:`exec/upp.x` executable and then rerun the :bolditalic:`compile_upp.sh` script. This is
   recommended if a mistake is made during the installation process.
   
