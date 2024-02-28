.. _FAQ:

******************************
Frequently Asked Questions
******************************

* :ref:`Is UPP compatible with NetCDF4? <netcdf4>`
* :ref:`How do I compile on another platform/compiler? <compile-other-platforms>`
* :ref:`How can I output satellite fields with the Unified Post Processor (UPP)? <satellite-fields>`
* :ref:`How do I add a new variable to UPP output? <new-variable>`
* :ref:`Why is the variable I requested not present in the UPP output? <var-not-present>`
* :ref:`If the UPP fails, how do I troubleshoot the problem? <troubleshooting>`
* :ref:`How do I regrid UPP output to another domain or projection? <regridding-faq>`
* :ref:`I am running UPP in parallel, but it fails. <upp-parallel>`
* :ref:`My FV3GFS unipost output is on a Gaussian grid. How can I process it to another grid such as a lat-lon grid or other user-defined grid? <process-grid>`
* :ref:`What does this warning mean in my compile.log? <nemsio-error>` ``libnemsio.a(nemsio_module_mpi.o): In function '__nemsio_module_mpi_MOD_readmpi4': nemsio_module_mpi.f90:(.text+0x1088): undefined reference to 'mpi_type_create_indexed_block_'``
* :ref:`Why do I see ** FATAL ERROR: Statistical processing bad n=0 ** when using the wgrib2 utility on my UPP output? <wgrib2-error>`


.. _netcdf4:

Is UPP compatible with NetCDF4?
=================================

The UPP is compatible with NetCDF4 when used on UFS model output.

.. _compile-other-platforms:

How do I compile on another platform/compiler?
================================================

We are not able to support all platform and compiler combinations out there but will try to help with specific issues when able. Users may request support on the UPP `GitHub Discussions <https://github.com/NOAA-EMC/UPP/discussions/categories/q-a>`__ page. We always welcome and are grateful for user-contributed configurations.

.. _satellite-fields:

How can I output satellite fields with the Unified Post Processor (UPP)?
==========================================================================

Currently, the standalone release of the UPP can be utilized to output satellite fields if desired. The UPP documentation :ref:`lists the grib2 fields <grib2-fields-by-id>`, including satellite fields, produced by the UPP. After selecting which fields to output, the user must :ref:`adjust the control file <control-file>` according to the instructions in the UPP documentation to output the desired fields. When outputting satellite products, users should note that not all physics options are supported for outputting satellite products. Additionally, for regional runs, users must ensure that the satellite field of view overlaps some part of their domain. 

Most UFS application releases do not currently support this capability, although it is available in the Short-Range Weather (SRW) Application. `This SRW App pull request (PR) <https://github.com/ufs-community/regional_workflow/pull/682>`__ added the option for users to output satellite fields using the SRW App. The capability is documented in the :ref:`SRW App User’s Guide <srw:SatelliteProducts>`.

.. _new-variable:

How do I add a new variable to UPP output?
============================================

If the desired variable is already available in the UPP code, then the user can simply add that variable to the ``postcntrl.xml`` file and :ref:`remake the postxconfig-NT.txt file <create_txt_file>` that the UPP reads. Please note that some variables may be dependent on the model and/or physics used.

If the desired variable is not already available in the UPP code, it can be added following the instructions for :ref:`adding a new variable <add-new-var>` in the UPP User’s Guide.

.. _var-not-present:

Why is the variable I requested not present in the UPP output?
================================================================

There are a few possible reasons why a requested variable might not appear in the UPP output:

#. The variable may be dependent on the model. 
#. Certain variables are dependent on the model configuration. For example, if a variable depends on a particular physics suite, it may not appear in the output when a different physics suite is used. 
#. The requested variable may depend on output from a different field that was not included in the model output.

.. _troubleshooting:

If the UPP fails, how do I troubleshoot the problem?
======================================================

If the user suspects that the UPP failed (e.g., no UPP output was produced or console output includes an error message like ``mv: cannot stat `GFSPRS.GrbF00`: No such file or directory``), the best way to diagnose the issue is to consult the UPP runtime log file for errors. When using the standalone UPP with the ``run_upp`` script, this log file will be located in the ``postprd`` directory under the name ``upp.fHHH.out``, where ``HHH`` refers to the 3-digit forecast hour being processed. When the UPP is used with the SRW App, the UPP log files can be found in the experiment directory under ``log/run_post_fHHH.log``.

.. _regridding-faq:

How do I regrid UPP output to another domain or projection?
=============================================================

UPP output is in standard grib2 format and can be interpolated to another grid using the third-party utility `wgrib2 <https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/new_grid.html>`__. Some basic examples can also be found in :numref:`Section %s <regridding>`.

.. _upp-parallel:

I am running UPP in parallel, but it fails.
==================================================================

This may be a memory issue; try increasing the number of CPUs or spreading them out across nodes (e.g., increase ``ptiles``). We also know of one version of MPI (mpich v3.0.4) that does not work with UPP. A work-around was found by modifying the ``UPP/sorc/ncep_post.fd/WRFPOST.f`` routine to change all ``unit 5`` references (which is standard I/O) to ``unit 4`` instead.

.. _process-grid:

My FV3GFS unipost output is on a Gaussian grid. How can I process it to another grid such as a lat-lon grid or other user-defined grid?
=============================================================================================================================================

For regridding grib2 unipost output, the wgrib2 utility can be used. See `complete documentation on grid specification with examples of regridding for all available grid definitions <https://www.cpc.ncep.noaa.gov/products/wesley/wgrib2/new_grid.html>`__. The :ref:`Regridding section <regridding>` of this UPP User’s Guide also gives examples (including an example from operations) of using wgrib2 to interpolate to various common grids.

.. _nemsio-error:

What does this warning mean in my compile.log? ``libnemsio.a(nemsio_module_mpi.o): In function '__nemsio_module_mpi_MOD_readmpi4': nemsio_module_mpi.f90:(.text+0x1088): undefined reference to 'mpi_type_create_indexed_block_'``
====================================================================================================================================================================================================================================

This warning appears for some platforms/compilers because a call in the *nemsio* library is never used or referenced for a serial build. This is just a warning and should not hinder a successful build of UPP or negatively impact your UPP run.

.. _wgrib2-error:

Why do I see ``** FATAL ERROR: Statistical processing bad n=0 **`` when using the wgrib2 utility on my UPP output?
=====================================================================================================================

This error message is displayed when using more recent versions of the wgrib2 utility on files for forecast hour zero that contain accumulated or time-averaged fields. This is due to the newer versions of wgrib2 no longer allowing the ``n`` parameter to be zero or empty. 

Users should consider using a separate control file (e.g., ``postcntrl_gfs_f00.xml``) for forecast hour zero that does not include accumulated or time-averaged fields, since they are zero anyway. Users can also continue to use an older version of *wgrib2*; v2.0.4 is the latest known version that does not result in this error.