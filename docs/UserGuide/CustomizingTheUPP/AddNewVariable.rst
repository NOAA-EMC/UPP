*********************
Adding a New Variable
*********************

This chapter provides general procedures and an example of how to add a new variable to the UPP code.
Please keep in mind that it may not be an exhaustive step-by-step process depending on your particular situation.
While we can provide general assistance for adding a new variable, users should be aware that this
requires good knowledge of Fortran and a thorough understanding of the code.

NOAA UPP developers who wish to add new variables to the UPP will need to:

1.  Read and follow procedures on the `UPP wiki page <https://github.com/NOAA-EMC/UPP/wiki/UPP-Code-Development>`__ on how to contribute your code changes to the UPP main development branch. Doing so will ensure your changes are merged
    to the UPP development branch quickly.

2.  Submit your pull request with small incremental changes. Advantages of doing this include avoiding conflicts with other UPP developers in terms of using the UPP internal index and variables.

3.  Please do not modify existing algorithms without coordinating with UPP code managers (Wen Meng and Hui-Ya Chuang). UPP supports many NOAA operational models, and we cannot change operational products without coordination and advanced notice.

We encourage non-NOAA UPP developers to contact EPIC via
`GitHub Discussions <https://github.com/NOAA-EMC/UPP/discussions>`_ to make them aware of modifications you
are making. In some cases, if they determine the changes you are making may be relevant for operational
and/or community purposes, they will be interested in incorporating your changes into the code base for
support and future release. We would then work with you to make this possible.

.. _add-var-process:

=========================================
Process Overview: Adding a New Variable
=========================================

The following steps outline the process for adding a new variable. This description is followed by a detailed
example in :numref:`Section %s <add-var-example>` below.

#. Check whether your new variable has been defined in the file ``parm/post_avblflds.xml`` in your UPP working
   directory. This file defines all available :term:`GRIB2` fields in the UPP. Users may also check the table showing 
   :doc:`../tables/UPP_GRIB2_Table_byID`.

   A. If NO (not available in ``post_avblflds.xml``), check whether your new variable has been defined in the
      `NCEP Grib2 Table <https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2.shtml>`__
      (Product Discipline and Parameter Category).

      \i. If NO (not available in the NCEP Grib2 Table):

         a. NOAA users can email Andrew.Benjamin@noaa.gov with the following information for your new
            variable: variable definition, unit, and what Grib2 discipline and category you think this
            variable should belong to. Andrew will define your new variable in the `NCEP Grib2 Table
            <https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2.shtml>`_ and
            inform you of the Grib2 discipline and category numbers you should use.

         b. Contact Andrew to update ``parm/params_grib2_tbl_new.text`` with your new variable and
            generate a ``params_grib2_tbl_new`` that lists variables in alphabetical order to improve post-processing
            efficiency.

         c. Save new ``params_grib2_tbl_new.text`` and ``params_grib2_tbl_new`` under ``parm/`` of your UPP
            working directory.

         d. Non-NOAA users should coordinate through EPIC for the above three steps. Users may post a
            `GitHub Discussions <https://github.com/NOAA-EMC/UPP/discussions/categories/enhancements>`__ 
            topic and tag @FernandoAndrade-NOAA and @gspetro-NOAA for directions in steps a-c. 

         e. Add a new entry in ``post_avblflds.xml`` with your new variable; then follow step B below, then step 2 and beyond. You should assign a new UPP ID for your new variable.

      \ii. If YES (variable is available in the NCEP Grib2 Table):

          a. Add a new entry in ``post_avblflds.xml`` with your new variable, then follow step B below, then step 2 and beyond. You should assign a new UPP ID for your new variable.

   B. If YES (variable is in ``post_avblflds.xml``), then your new variable is already available in the UPP. 
      Follow steps i) and ii), make a test UPP run, and then look for your new variable in your output.
      You can skip the remaining steps about modifying the source code.

      \i. Add a new entry in your application’s control xml file (e.g., ``fv3lam.xml`` for the FV3LAM application, ``postcntrl_gfs.xml`` for the FV3GFS application). This file lets users control which variables to output from UPP for Grib2.

      \ii. Generate ``your_user_defined_flat`` file (e.g., ``postxconfig-NT-fv3lam.txt`` for FV3LAM application) by executing:

         .. code-block:: console

            /usr/bin/perl PostXMLPreprocessor.pl your_user_defined_xml post_avblflds.xml your_user_defined_flat

         This flat file (instead of the xml file) is read in by the UPP because it is much faster to read a text file than an XML file.

#. Allocate and initialize the field in ``sorc/ncep_post.fd/ALLOCATE_ALL.f``.

   This file contains the instantiation or allocation of each variable. Note that the variables are defined
   based on the parallel processing capability of UPP. Use an example from the file.

#. Deallocate the field in ``sorc/ncep_post.fd/DEALLOCATE.f``.

   All good programmers give back their resources when they are done. Please update this routine to
   return your resource to the system.

#. Declare the new variable in ``VRBLS2D_mod.f``, ``VRBLS3D_mod.f``, or ``VRBLS4D_mod.f``.
    
   The variable is declared in one of these module-defining files depending on its dimension.

#. Read model output if necessary using ``INITPOST_NETCDF.f``.

   Check first to see if all variables needed to derive your new variable are already available in the UPP. If not,
   you will need to use this file (or another appropriate ``INITPOST_*.f`` file) for reading the model output files. 
   The appropriate one should be chosen based on the model and the model output format.

#. Add to appropriate routine(s) for filling the new variable (e.g., ``SURFCE.f``, ``MDLFLD.f``, ``MDL2P.f``).

   This is the place where you will derive your new variable and then fill the Grib2 array with the data to be
   written out later on.

#. Build or rebuild the code for changes to take effect before running your UPP run script.

.. _add-var-example:

===========================================================
Example Procedure: Steps for adding a new variable ‘TG3’
===========================================================

This example adds TG3 to the UPP. TG3 is the averaged climatology of surface temperature, which the land surface models (LSMs) use to specify bottom soil temperature, where the depth of the bottom is LSM-dependent. For this example, a depth of 500cm is used.

- This example illustrates adding a new variable from GFS output that will be read into UPP
  and directly output into the Grib2 output files (i.e., no additional computations/calculations
  are needed for the field).
- Additions to each of the routines are highlighted. 
- Locations of routines are in ``UPP/sorc/ncep_post.fd`` unless specified otherwise.
- The new variable, TG3, added in this example is found in the ``gfs.t00z.sfcf006.nc`` file; however, both the
  ``gfs.t00z.sfcf006.nc`` and ``gfs.t00z.atmf006.nc`` output files are required to run UPP for GFS.

  New variable to add::

   float tg3(time, grid_yt, grid_xt) ;
         tg3:long_name = "deep soil temperature" ;
         tg3:units = "K" ;
         tg3:missing_value = 9.99e+20 ;
         tg3:cell_methods = "time: point" ;
         tg3:output_file = "sfc" ;

1. Check whether your new variable has been defined in the file ``parm/post_avblflds.xml`` in your UPP working
   version.

   A. This variable is not available in ``parm/post_avblflds.xml``.

      \i. Check whether your new variable has been defined in the NCEP Grib2 Table.

         1) This variable is not defined in the NCEP Grib2 Table.

            a)-d) For the purpose of this example alone, steps a) - d) are not executed as instructed.
               Instead, manual instructions are provided here for adding to the ``params_grib2_table_new`` in order
               to create a working example. 

               For this example, the variable will be added to ``parm/params_grib2_tbl_new`` manually. You would only
               do this if you had no plans to contribute your addition to the UPP ``develop`` branch; otherwise, follow the
               instructions as a NOAA or Other user in steps a) - d). 
 
               For all current UPP output fields, the ``params_grib2_table_new`` lists, in order, the following attributes:
                - Discipline (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table0-0.shtml)
                - Category (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-1.shtml)
                - Parameter Number (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2.shtml)
                - Table information (0 for parameters from the WMO table; 1 for parameters from the local NCEP table)
                - Abbreviated Variable Name (from the parameters table)

               User Procedure
                - Add this variable as TG3.
                - TG3 is a land surface product (discipline=2)
                - TG3 is a vegetation/biomass product (category=0)
                - Pick an unused parameter number from the table defined by discipline=2 and category=0
                  (Table 4.2-0-0: https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2-2-0.shtml). 
                  The parameter number should not be in use in Table 4.2 or the current ``params_grib2_tbl_new``.
                  In this case, the unused parameter number 251 was chosen.
                - Add using the NCEP local table (table=1)
                - Choose an abbreviated parameter name to describe your field (e.g., TG3)
                - Add alphabetically (by variable name) to the table as:
      
                ::

                 2 0 251 1 TG3

            e) **Add the new variable to** ``UPP/parm/post_avblflds.xml``, **which lists all fields available
               for output in GRIB2 format.** This file is generally not modified unless adding a new field or
               modifying an existing one. Users should indicate the following variable attributes in the XML file:

                - ``post_avblfldidx``: the unique array index number used to store this variable. The number chosen here
                  is just an example, and it is important to pick one that is not yet in use.
                - ``shortname``: name describing the variable and level type
                - ``pname``: the abbreviation for your variable (should match what is used in ``params_grib2_tbl_new``)
                - ``table info``: table used if not standard WMO
                - ``fixed_sfc1_type``: level type
                - ``level``: generally only used here if it is a fixed level specific to the variable (e.g., T2m, TSOIL5m)
                - ``scale``: precision of data written out to Grib2 file

               User procedure
                - Add as:
      
                ::

                 <param>
                   <post_avblfldidx>1063</post_avblfldidx>
                   <shortname>DEEP_TSOIL_ON_DEPTH_BEL_LAND_SFC</shortname>
                   <pname>TG3</pname>
                   <fixed_sfc1_type>depth_bel_land_sfc</fixed_sfc1_type>
                   <table_info>NCEP</table_info>
                   <level>500.</level>
                   <scale>3.0</scale>
                 </param>

   B. Add the variable to the user-defined control file.

      i. Add a new entry in your application's control XML file (e.g., ``fv3lam.xml`` for the FV3LAM application,
         ``postcntrl_gfs.xml`` for the ``FV3GFS`` application). This file lets users control which variables to output
         from the UPP for Grib2.

         User procedure
          - Add as:

          ::

           <param>
             <shortname>DEEP_TSOIL_ON_DEPTH_BEL_LAND_SFC</shortname>
             <scale>4.0</scale>
           </param>

      ii. Generate ``your_user_defined_flat`` file (e.g., ``postxconfig-NT-fv3lam.txt`` for the FV3LAM application) by
          executing:

          ::

           >> /usr/bin/perl PostXMLPreprocessor.pl your_user_defined_xml post_avblflds.xml your_user_defined_flat

          This flat file (instead of the XML file) is read in by the UPP.

2. Allocate and initialize the new variable in ``ALLOCATE_ALL.f`` using an example from the file.
   Note that the variables are defined based on the parallel processing capability of the UPP. 

   User Procedure
    - Allocate in the *VRBLS2D* GFS section of ``ALLOCATE_ALL.f`` as:

    ::

      allocate(tg3(ista_2l:iend_2u,jsta_2l:jend_2u))
      
    - Initialize TG3 in the initialization section that comes after the allocation section you added to.

    ::

      tg3(i,j)=spval

3. De-allocate the variable to give the resources back in ``DEALLOCATE.f``.
   Updating this routine returns your resources to the system.

   User procedure
    - Add in *VRBLS2D* GFS section of ``DEALLOCATE.f`` as:
      
    ::

     deallocate(tg3)

4. Declare the new variable in the appropriate file (e.g., ``VRBLS2D_mod.f``, 
   ``VRBLS3D_mod.f``, or ``VRBLS4D_mod.f``) depending on its dimensions.

   User procedure
    - TG3 is a 2-dimensional field, so declare it in ``VRBLS2D_mod.f``.
    - Add to the GFS section as:
      
    ::

     tg3(:,:)

5. Read the field from the GFS model output file by adding the new variable into ``INITPOST_NETCDF.f``.
   This file is used for reading the GFS model FV3 output files in parallel netCDF format.

   User procedure
    - Add to top section of the routine in the ‘use vrbls2d’ section to initiate the new variable as:
      
    ::

     tg3

    - Read in the new variable in the section for reading the 2D netCDF file using another 2D variable
      as an example, such as ``hpbl``. Add as:
      
    ::

     ! deep soil temperature
           VarName='tg3'
           call read_netcdf_2d_para(ncid2d,ista,ista_2l,iend,iend_2u,jsta,jsta_2l,jend,jend_2u, &
           spval,VarName,tg3)

6. Determine the appropriate routine to add the new variable to (e.g., ``SURFCE.f``, ``MDLFLD.f``,
   ``MDL2P.f``). The appropriate routine will depend on what your field is. 
   For example, if you have a new diagnostic called *foo*,
   and you want it interpolated to pressure levels, you would need to add it to ``MDL2P.f``. If *foo* were only a
   surface variable, you would add it to ``SURFCE.f``. If you wanted *foo* on native model levels, you
   would add it to ``MDLFLD.f``. If you are not sure which routine to add the new variable to, choose a
   similar variable as a template, and add it in the same places.

   .. note:: 
      
      This is also where you would add any calculations needed for your new variable, should they
      be required.

   User procedure
    - Treat TG3 like a surface field, similar to the other soil fields, and add it to ``SURFCE.f``.
    - Use another 2D variable, such as 'SNOW WATER EQUIVALENT' as a template. This variable is also
      being read through and output, similar to what we want.
    - Add to top section in ‘use vrbls2d, only’ to initiate the new variable as:
      
    ::

     tg3

    - Add in main section using a template variable as a guide.

    ::

     ! DEEP SOIL TEMPERATURE
     IF ( IGET(1063).GT.0 ) THEN
       ID(1:25) = 0
       If(grib=='grib2') then
         cfld=cfld+1
         fld_info(cfld)%ifld=IAVBLFLD(IGET(1063))
     !$omp parallel do private(i,j,jj)
         do j=1,jend-jsta+1
           jj = jsta+j-1
           do i=1,iend-ista+1
           ii = ista+i-1
             datapd(i,j,cfld) = TG3(ii,jj)
           enddo
         enddo
       endif
     ENDIF

7. Build or rebuild the code for changes to take effect before running your UPP run script.
   
   User procedure for building on pre-configured machines: 

    ::

    >> cd UPP/tests
    >> ./compile_upp.sh

   Assuming the modified code built successfully, and you were able to produce Grib2 output, you can check the Grib2
   file for your new variable.

   **GRIB2 output of the new variable from this example procedure (using the wgrib2 utility if available on your system):**

    ::

     wgrib2 -V GFSPRS.006

     716:37731711:vt=2019061506:500 m underground:6 hour fcst:var discipline=2 center=7 local_table=1 parmcat=0 parm=251:
         ndata=73728:undef=0:mean=278.383:min=215.47:max=302.4
         grid_template=40:winds(N/S):
         Gaussian grid: (384 x 192) units 1e-06 input WE:NS output WE:SN
         number of latitudes between pole-equator=96 #points=73728
         lat 89.284225 to -89.284225
         lon 0.000000 to 359.062500 by 0.937500

   - For this example, since the new variable was not added to the NCEP Grib2 table, it will not be defined by the
     variable name. Instead it will be defined using the Grib2 parameter information entered into ``params_grib2_tbl_new``
     from step 1 of this procedure.