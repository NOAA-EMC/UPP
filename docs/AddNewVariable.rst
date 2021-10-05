*********************
Adding a New Variable
*********************

This document provides general procedures and an example of how to add a new variable to the UPP code.
Please keep in mind it may not be an exhaustive step-by-step depending on your particular situation.
While we can provide general assistance for adding a new variable, users should be aware that this
requires good knowledge of Fortran and thorough understanding of the code.

NOAA UPP developers who wish to add new variables to the UPP, please follow the following procedures:

1.  Read and follow procedures on the `UPP wiki page <https://github.com/NOAA-EMC/UPP/wiki/UPP-Code-Development>`_
    on how to contribute your code changes to the UPP main development. Doing so will ensure your changes are merged
    to the UPP development quickly.

2.  Submit your pull request with small incremental changes. Advantages of doing this include avoiding
    conflicts with other UPP developers in terms of using the UPP internal Index and variables.

3.  Please do not modify existing algorithms without coordinating with UPP code managers (Wen Meng and
    Hui-Ya Chuang). UPP supports many NOAA operational models and we can not change operational products
    without coordination and advanced notices.

We encourage non NOAA UPP developers to contact the Developmental Testbed Center (DTC) via the UPP
`forum <https://forums.ufscommunity.org/forum/post-processing>`_ to make them aware of modifications you
are making. In some cases, if they determine the changes you are making may be relevant for operational
and/or community purposes, they will be interested in incorporating your changes into the code base for
support and future release. We would then work with you to make this possible.

The following outlines a brief description of the steps to be taken and are described in more detail
with examples in the sections below.

1.  Check whether your new variable has been defined in the file parm/post_avblflds.xml in your UPP working
    version. This file defines all available GRIB2 fields in UPP.

    A.  If no (not available in post_avblflds.xml)

        i.  Check whether your new variable has been defined in the
            `NCEP Grib2 Table <https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2.shtml>`_
            (Product Discipline and Parameter Category).

            1.  If no (not available in the NCEP Grib2 Table)

                a.  NOAA users can email Boi.Vuong@noaa.gov with the following information for your new
                    variable: variable definition, unit, and what Grib2 discipline and category you think this
                    variable should belong to. Boi will define your new variable in the `NCEP Grib2 Table
                    <https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2.shtml>`_ and
                    inform you of the Grib2 discipline and category numbers you should use.

                b.  Contact with Boi to update parm/params_grib2_tbl_new.text with your new variable and
                    generate a params_grib2_tbl_new which lists in alphabetical order to improve post
                    processing efficiency.

                c.  Save new params_grib2_tbl_new.text and params_grib2_tbl_new under parm/ of your UPP
                    working version.

                d.  Other users please coordinate through the DTC for the above three steps.

                e.  Add a new entry in post_avblflds.xml with your new variable, then follow step B), then step 2)
                    and beyond. You should assign a new UPP ID for your new variable.

            2.  If yes (available in the NCEP Grib2 Table)

                a.  Add a new entry in post_avblflds.xml with your new variable, then follow step B), then step 2)
                    and beyond. You should assign a new UPP ID for your new variable.

    B.  If yes (available in post_avblflds.xml), then your new variable is already available in UPP. Follow the
        following two steps i) and ii), make a test UPP run, and then look for your new variable in your output.
        You can skip the remaining steps about modifying the source codes.

        i.  Add a new entry in your application’s control xml file (e.g. fv3lam.xml for FV3LAM application,
            postcntrl_gfs.xml for FV3GFS application). This file lets users control which variables to output
            from UPP for Grib2.

        ii. Generate your_user_defined_flat file (e.g. postxconfig-NT-fv3lam.txt for FV3LAM application) by
            executing:

            /usr/bin/perl PostXMLPreprocessor.pl your_user_defined_xml post_avblflds.xml your_user_defined_flat

            This flat file (instead of the xml file) is read in by UPP as it was much faster to read a text file
            than an xml file.

2.  Define the new field: RQSTFLD.F

    This file contains a list of all possible fields to be output by UPP, corresponding Grib1 key-word character
    strings, UPP ID (Index) for internal code, and Grib1 IDs. Note that as of December 2020, EMC removed the Grib1
    option from its repository as part of its re-reginnering effort. Users will continue to use UPP ID (Index) for
    defining the new variable to be added. However, the character string and Grib1 ID in RQSTFLD.F will no longer
    be used by UPP.

3.  Allocate the field: ALLOCATE.f

    This file is the instantiation or allocation of the variable. Note that the variables are defined
    based on the parallel processing capability of UPP - use an example from the file.

4.  Deallocate the field: DEALLOCATE.f

    All good programmers give back their resources when they are done. Please update this routine to
    return your resource to the system.

5.  Declare the new variable: VRBLS2D_mod.f, VRBLS3D_mod.f, or VRBLS4D_mod.f
    
    The variable is declared in one of these modules defining files depending on its dimension.

6.  Read model output if necessary: INITPOST_GFS_NETCDF_PARA.f (current operational netcdf output with GFS V16),
    INITPOST_NETCDF.f (LAM FV3 netcdf)

    Check first to see if all variables needed to derive your new variable are already available in UPP. If not,
    you will need to use one of these files for reading the model output files. The appropriate one will need to
    be chosen based on the model and model output format.

7.  Add to appropriate routine for filling the new variable: e.g. SURFCE.f, MDLFLD.f, MDL2P.f, etc

    This is the place where you will derive your new variable and then fill the Grib2 array with the data to be
    written out later on.

8. Build or rebuild the code for changes to take effect before running your UPP run script.


**Example Procedure: Steps for adding a new variable ‘TG3’**

- This example illustrates adding a new variable from GFS output that will be read into UPP
  and directly output into the Grib2 output files (i.e. in this case no additional computations/calculations
  are needed for the field).
- Additions to each of the routines are highlighted. 
- Locations of routines are in UPP/sorc/ncep_post.fd unless specified otherwise.
- The new variable, TG3, added in this example is found in the gfs.t00z.sfcf006.nc; however, both the
  gfs.t00z.sfcf006.nc and gfs.t00z.atmf006.nc output files are required to run UPP for GFS.
- TG3 is the averaged climatology of surface temperature, which the LSMs use to specify bottom soil T,
  where the depth of the bottom is LSM dependent. For this example, a depth of 500cm is used.

  New variable to add::

   float tg3(time, grid_yt, grid_xt) ;
         tg3:long_name = "deep soil temperature" ;
         tg3:units = "K" ;
         tg3:missing_value = 9.99e+20 ;
         tg3:cell_methods = "time: point" ;
         tg3:output_file = "sfc" ;

1. Check whether your new variable has been defined in the file parm/post_avblflds.xml in your UPP working
   version.

   A. This variable is not available in parm/post_avblflds.xml.

      i. Check whether your new variable has been defined in the NCEP Grib2 Table.

         1) This variable is not defined in the NCEP Grib2 Table.

            a)-d) For the purpose of this example alone, steps a) - d) are not executed as instructed.
               Instead, manual instructions are provided here for adding to the params_grib2_table_new in order
               to create a working example. 

               For this example, the variable will be added to parm/params_grib2_tbl_new manually. You would only
               do this if you had no plans to contribute your addition to UPP develop, otherwise, follow the
               instructions as a NOAA or Other user in steps a) - d). 
 
               For all current UPP output fields, the params_grib2_table_new lists, in order, the:
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
                  The parameter number should not be in use in table 4.2 or the current params_grib2_tbl_new.
                  In this case, the unused parameter number 231 was chosen.
                - Add using the NCEP local table (table=1)
                - Choose an abbreviated parameter name to describe your field (e.g. TG3)
                - Add alphabetically to the table as:
      
                ::

                 2 0 231 1 TG3

            e) Add the new variable to the UPP/parm/post_avblflds.xml, which lists all fields available
               for output in GRIB2 format. This file is generally not modified unless adding a new field or
               modifying an existing one.
                - Post_avblfldidx: the unique array number given in the RQSTFLD.f routine.
                - Shortname: name describing the variable and level type
                - Pname: the abbreviation for your variable (should match what is used in params_grib2_tbl_new)
                - Table info: table used if not standard WMO
                - Fixed_sfc1_type: level type
                - Level: Generally only used here if it's a fixed level specific to the variable (e.g. T2m, TSOIL5m)
                - Scale: precision of data written out to Grib2 file

               User procedure
                - Add as:
      
                ::

                 <param>
                   <post_avblfldidx>999</post_avblfldidx>
                   <shortname>DEEP_TSOIL_ON_DEPTH_BEL_LAND_SFC</shortname>
                   <pname>TG3</pname>
                   <fixed_sfc1_type>depth_bel_land_sfc</fixed_sfc1_type>
                   <table_info>NCEP</table_info>
                   <level>500.</level>
                   <scale>3.0</scale>
                 </param>

   B. Add the variable to the user defined control file.

      i. Add a new entry in your application’s control xml file (e.g. fv3lam.xml for FV3LAM application,
         postcntrl_gfs.xml for FV3GFS application). This file lets users control which variables to output
         from UPP for Grib2.

         User procedure
          - Add as:

          ::

           <param>
             <shortname>DEEP_TSOIL_ON_DEPTH_BEL_LAND_SFC</shortname>
             <scale>4.0</scale>
           </param>

      ii. Generate your_user_defined_flat file (e.g. postxconfig-NT-fv3lam.txt for FV3LAM application) by
          executing:

          ::

           >> /usr/bin/perl PostXMLPreprocessor.pl your_user_defined_xml post_avblflds.xml your_user_defined_flat

          This flat file (instead of the xml file) is read in by UPP as it was much faster to read a text file
          than an xml file.

2. Define the new variable in RQSTFLD.F which includes a list of all possible fields to be output by
   UPP, corresponding Grib1 key-word character strings, UPP ID (Index) for internal code, and Grib1 IDs.
   Ensure your code is up-to-date and pick a unique identifier that is not already used for the new variable.
   Currently, the 900's are being used for new contributions.

   Example Entry

       | ! HWRF addition for v_flux as pass through variable:

       |   DATA IFILV(901),AVBL(901),IQ(901),IS(901),AVBLGRB2(901) &
       |   &            /1,'MODEL SFC V WIND STR’,125,001,         &
       |   &            'V_FLX ON surface’/

   Where:
     - **IFILV** Identifies field as MASS/VELOCITY point (e.g. 1)
     - **AVBL** is the model output character string variable name for Grib1 (e.g. MODEL SFC V WIND STR)
     - **IQ** is the GRIB PDS OCTET 9 (table 2) - Indicator of parameter and units (e.g. 125)
     - **IS** is the GRIB PDS OCTET 10 (table 3&3a) - Indicator of type of level or layer (e.g. 001)
     - **AVBLGRB2** is the model output character string variable name for Grib2 (e.g. V_FLX ON surface)
     - A UNIQUE array Index UPP uses to store this variable in parallel arrays (e.g. **901**)

   User procedure
    - Soil temperature (TSOIL) is found in the Grib1 parameter tables as parameter number 085, so this
      can be used for the Grib1 ID.
      http://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html
    - Use level type 'depth below land surface', which is 111.
      http://www.nco.ncep.noaa.gov/pmb/docs/on388/table3.html
    - New variables are continuously being added to UPP, so be sure to check that the UPP Index 999 is
      still available before using it to add your new variable. If it is already in use, pick the next
      available Index.
    - Add as:

    ::

     DATA IFILV(999),AVBL(999),IQ(999),IS(999),AVBLGRB2(999) &
     &          /1,'DEEP SOIL TMP',085,111,                  &
     &          'DEEP TSOIL ON depth_bel_land_sfc'/

   .. note::
      Since Grib1 is no longer supported, the variable character strings and Grib IDs for Grib1 are not
      important, but still need to be included here for correct formatting.

3. Allocate the new variable in ALLOCATE_ALL.f
   This file is the instantiation or allocation of the variable. Note that the variables are defined
   based on the parallel processing capability of UPP - use an example from the file.

   User Procedure
    - Add in VRBLS2D GFS section as:

    ::

      allocate(tg3(im,jsta_2l:jend_2u))

4. De-allocate the variable to give the resources back in DEALLOCATE.f
   All good programmers give back their resources when they are done. Please update this
   routine to return your resources to the system.

   User procedure
    - Add in VRBLS2D GFS section as:
      
    ::

     deallocate(tg3)

5. Declare the new variable in the appropriate file depending on its dimensions;
   VRBLS2D_mod.f, VRBLS3D_mod.f or VRBLS4D_mod.f

   User procedure
    - tg3 is a 2-dimensional field, so declare it in VRBLS2D_mod.f
    - Add to the GFS section for adding new fields as:
      
    ::

     tg3(:,:)

6. Read the field from the GFS model output file by adding the new variable into INITPOST_GFS_NETCDF_PARA.f.
   This file is used for reading the GFS model FV3 output files in netcdf format.

   User procedure
    - Add to top section of the routine in ‘use vrbls2d’ to initiate the new variable as:
      
    ::

     tg3

    - Read in the new variable in the section for reading the 2D netcdf file using another 2D variable
      as an example, such as 'hpbl'. Add as:
      
    ::

     ! deep soil temperature
           VarName='tg3'
           call read_netcdf_2d_para(ncid2d,im,jsta,jsta_2l,jend,jend_2u, &
           spval,VarName,tg3)

7. Determine the appropriate routine to add the new variable to (e.g. SURFCE.f, MDLFLD.f,
   MDL2P.f, etc). This is the place that you will fill the Grib2 array with the data to be written out later on.
   The appropriate routine will depend on what your field is. For example, if you have a new diagnostic called foo,
   and you want it interpolated to pressure levels, you would need to add it to MDL2P.f. If foo was only a
   surface variable, you would add it to SURFCE.f. If you wanted foo on native model levels, you
   would add it to MDLFLD.f. If you’re not sure which routine to add the new variable to, choose a
   similar variable as a template.

   Note: This is also where you would add any calculations needed for your new variable, should it
   be required.

   User procedure
    - Treat tg3 like a surface field (SURFCE.f), similar to the other soil fields.
    - Use another 2D variable, such as 'SNOW WATER EQUIVALENT' as a template. This variable is also
      being read through and output, similar to what we want.
    - Add to top section in ‘use vrbls2d, only’ to initiate the new variable as:
      
    ::

     tg3

    - Add in main section using a template variable as a guide.

    ::

     ! DEEP SOIL TEMPERATURE
     IF ( IGET(999).GT.0 ) THEN
       ID(1:25) = 0
       If(grib=='grib2') then
         cfld=cfld+1
         fld_info(cfld)%ifld=IAVBLFLD(IGET(999))
     !$omp parallel do private(i,j,jj)
         do j=1,jend-jsta+1
           jj = jsta+j-1
           do i=1,im
             datapd(i,j,cfld) = TG3(i,jj)
           enddo
         enddo
       endiF
     ENDIF

8. Build or rebuild the code for changes to take effect before running your UPP run script.
   
   User procedure IF you already have the code built. Otherwise, see the User's Guide for instructions on building.

    ::

    >> cd UPP/build
    >> make install

   Assuming the modified code built successfully and you were able to produce Grib2 output, you can check the Grib2
   file for your new variable.

   GRIB2 output of the new variable from this example procedure (using the wgrib2 utility if available on your system).
    - For this example, since the new variable was not added to the NCEP Grib2 table, it will not be defined by the
      variable name. Instead it will be defined using the Grib2 parameter information entered into params_grib2_tbl_new
      from step 1 of this procedure.

    ::

     wgrib2 -V GFSPRS.006

     716:37731711:vt=2019061506:500 m underground:6 hour fcst:var discipline=2 center=7 local_table=1 parmcat=0 parm=231:
         ndata=73728:undef=0:mean=278.383:min=215.47:max=302.4
         grid_template=40:winds(N/S):
         Gaussian grid: (384 x 192) units 1e-06 input WE:NS output WE:SN
         number of latitudes between pole-equator=96 #points=73728
         lat 89.284225 to -89.284225
         lon 0.000000 to 359.062500 by 0.937500

