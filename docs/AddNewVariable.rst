*********************
Adding a New Variable
*********************

This document provides general procedures and an example of how to add a new variable to the UPP code.
Please keep in mind it may not be an exhaustive step-by-step depending on your particular situation.
While we can provide general assistance for adding a new variable, users should be aware that this
requires good knowledge of Fortran and thorough understanding of the code.

We encourage users to contact us via the UPP `forum <https://forums.ufscommunity.org/forum/post-processing>`_
to make us aware of modifications you are making. In some cases, if we determine the changes you are
making may be relevant for operational and/or community purposes, we will be interested in incorporating
your changes into the code base for support and future release. We would then work with you to make this
possible.

The following outlines a brief description of the steps to be taken and are described in more detail
with examples in the sections below.

1.  Allocate the field: ALLOCATE.f

    *This file is the instantiation or allocation of the variable. Note that the variables are defined
    based on the parallel processing capability of UPP - use an example from the file.*

2.  Deallocate the field: DEALLOCATE.f

    *All good programmers give back their resources when they are done. Please update this routine to
    return your resource to the system.*

3.  Declare the new variable: VRBLS2D_mod.f, VRBLS3D_mod.f, or VRBLS4D_mod.f
    
    *The variable is declared in one of these modules defining files depending on its dimension.*

4. Define field for grib1: RQSTFLD.f

   *This file contains a list of all possible fields to be output by UPP, corresponding
   key-word character string user places in wrf_cntrl.parm file, UPP ID for internal
   code, and grib IDs.*

5.  Read model output: INITPOST_GFS_NEMS_MPIIO.f (GFS nemsio), INITPOST_GFS_NETCDF (GFS netcdf),
    INITPOST_NETCDF (LAM netcdf)

    *These files are used for reading the model output files. The appropriate one will need to be
    chosen based off the model and model output format.*

6.  Add to appropriate routine for filling the new variable: e.g. SURFCE.f, MDLFLD.f, MDL2P.f, etc

    *This is the place that you will fill the array with the data and output the field.*

7.  Define table/grib2 parameters for grib2 output: params_grib2_tbl_new

    *This table contains the necessary parameter information for grib2 fields.*

8.  Define the field for grib2 output: post_avlbflds.xml

    *This file is used for defining all available grib2 fields.*

9.  Define control file entry for output: postcntrl.xml & postxconfig-NT.txt

    *These files are used for controlling which fields are output by UPP for grib2.*


**Example Procedure: Steps for adding a new variable ‘TG3’**

- This example illustrates adding a new variable from GFS output that will be read into UPP
  and directly output into the Grib2 output files (i.e. no additional computations/calculations
  are needed for the field).
- Additions to each of the routines are highlighted. 
- Locations of routines are in /EMC_post/sorc/ncep_post.fd unless specified otherwise.
- Sample GFS files for the following procedures are available for download
  `here <https://dtcenter.org/sites/default/files/community-code/upp/AddNewVar_GFSdata.tar.gz>`_.
 - This data is the 6-hr forecast of a GFS initialization of 2019-06-15_00:00:00
 - The new variable, TG3, added in this example is found in the sfcf006.nc; however, both the sfcf006.nc
   and atmf006.nc output files are required to run UPP for GFS.
 - TG3 is the averaged climatology of surface temperature, which the LSMs use to specify bottom soil T,
   where the depth of the bottom is LSM dependent. For this example, a depth of 500cm is used.

   New variable to add::

    float tg3(time, grid_yt, grid_xt) ;
          tg3:long_name = "deep soil temperature" ;
          tg3:units = "K" ;
          tg3:missing_value = 9.99e+20 ;
          tg3:cell_methods = "time: point" ;
          tg3:output_file = "sfc" ;

1. Allocate the new variable in ALLOCATE_ALL.f
   This file is the instantiation or allocation of the variable. Note that the variables are defined
   based on the parallel processing capability of UPP - use an example from the file.

   User Procedure
    - Add in VRBLS2D GFS section as:

    ::

      allocate(tg3(im,jsta_2l:jend_2u))

2. De-allocate the variable to give the resources back in DEALLOCATE.f
   All good programmers give back their resources when they are done. Please update this
   routine to return your resources to the system.

   User procedure
    - Add in VRBLS2D GFS section as:
      
    ::

     deallocate(tg3)

3. Declare the new variable in the appropriate file depending on its dimensions;
   VRBLS2D_mod.f, VRBLS3D_mod.f or VRBLS4D_mod.f

   User procedure
    - tg3 is a 2-dimensional field, so declare it in VRBLS2D_mod.f
    - Add to the GFS section for adding new fields as:
      
    ::

     tg3(:,:)

4. List the new variable in RQSTFLD.F which includes a list of all possible fields to be output by
   UPP, as well as the corresponding UPP ID for internal code, variable character stings, and grib IDs.
   Be sure to pick a unique identifier that is not already used for the new variable. Right now, the
   900's are being used for community contributions.

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
     - A UNIQUE array location UPP uses to store this variable in parallel arrays (e.g. **901**)

   User procedure
    - Soil temperature (TSOIL) is found in the Grib1 parameter tables as parameter number 085, so this
      can be used for the Grib1 ID.
      http://www.nco.ncep.noaa.gov/pmb/docs/on388/table2.html
    - Use level type 'depth below land surface', which is 111.
      http://www.nco.ncep.noaa.gov/pmb/docs/on388/table3.html
    - Add as:

    ::

     DATA IFILV(979),AVBL(979),IQ(979),IS(979),AVBLGRB2(979) &
     &          /1,'DEEP SOIL TMP',085,111,                  &
     &          'DEEP TSOIL ON depth_bel_land_sfc'/

   .. note::
      Since Grib1 is no longer supported, the variable character strings and Grib IDs for Grib1 are not
      important, but still need to be included here for correct formatting.

5. Read the field from the GFS model output file by adding the new variable into INITPOST_GFS_NETCDF.f.
   This file is used for reading the GFS model output files in netcdf format.

   User procedure
    - Add to top section of the routine in ‘use vrbls2d’ to initiate the new variable as:
      
    ::

     tg3

    - Read in the new variable in the section 'start reading 2D netcdf file' using another 2D variable
      as an example, such as 'hpbl'. Add as:
      
    ::

     ! deep soil temperature
           VarName='tg3'
           call read_netcdf_2d_scatter(me,ncid2d,1,im,jm,jsta,jsta_2l &
            ,jend_2u,MPI_COMM_COMP,icnt,idsp,spval,VarName,tg3)

6. Determine the correct routine to add the new variable to (e.g. SURFCE.f, MDLFLD.f,
   MDL2P.f, etc). You will need to determine the correct routine to add your field into; this is the
   place that you will fill the array with the data and output the field. The correct
   routine will depend on what your field is. For example, if you have a new diagnostic called foo, and
   you want it interpolated to pressure levels, you would need to add it to MDL2P.f. If foo was only a
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
     IF ( IGET(979).GT.0 ) THEN
       ID(1:25) = 0
       If(grib=='grib2') then
         cfld=cfld+1
         fld_info(cfld)%ifld=IAVBLFLD(IGET(979))
     !$omp parallel do private(i,j,jj)
         do j=1,jend-jsta+1
           jj = jsta+j-1
           do i=1,im
             datapd(i,j,cfld) = TG3(i,jj)
           enddo
         enddo
       endiF
     ENDIF

   .. note::
      Since Grib1 is no longer supported, the if-statement for filling the grid for this output type is
      removed here and is only filled for Grib2 output.

7. Add the new variable to /EMC_post/parm/params_grib2_tbl_new.
   For all current UPP output fields, this table lists, in order, the:
    - Discipline (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table0-0.shtml)
    - Category (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-1.shtml)
    - Parameter Number (https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2.shtml)
    - Table information (0 for parameters from the WMO table; 1 for parameters from the local
      NCEP table)
    - Abbreviated Variable Name (from the parameters table)

   User Procedure
    - Here we could just use TSOIL, which is already in the table; howerver, instead we will add this
      using a new name, TG3, to demonstrate this step.
    - TG3 is a land surface product (discipline=2)
    - TG3 is a vegetation/biomass product (category=0)
    - Pick an unused parameter number from the table defined by discipline=2 and
      category=0 (Table 4.2-0-0: https://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_doc/grib2_table4-2-2-0.shtml). 
      The parameter number should not be in use in table 4.2 or the current
      params_grib2_tbl_new. In this case, the unused parameter number 231 was chosen.
    - Add using the NCEP local table (table=1)
    - Choose an abbreviated parameter name to describe your field (e.g. TG3)
    - Add as:
      
    ::

     2 0 231 1 TG3

8. Add the new variable to the /EMC_post/parm/post_avblflds.xml, which lists all fields available
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
       <post_avblfldidx>979</post_avblfldidx>
       <shortname>DEEP_TSOIL_ON_DEPTH_BEL_LAND_SFC</shortname>
       <pname>TG3</pname>
       <fixed_sfc1_type>depth_bel_land_sfc</fixed_sfc1_type>
       <table_info>NCEP</table_info>
       <level>500.</level>
       <scale>3.0</scale>
     </param>

9. Add the new variable to the /EMC_post/parm/postcntrl_gfs.xml file, which lists all fields and levels
   you wish to output for GRIB2. Remake the /EMC_post/parm/postxconfig-NT-GFS.txt file, which is read by
   UPP and contains the information from the xml.
    - See the User’s guide on steps for creating the text control file
   
   User procedure
    - Add as:
      
    ::

     <param>
       <shortname>DEEP_TSOIL_ON_DEPTH_BEL_LAND_SFC</shortname>
       <scale>4.0</scale>
     </param>

10. Build or rebuild the code to include the changes before running your UPP run script.
   
    User procedure IF you already have the code built. Otherwise, see the User's Guide for instructions
    on building.

    ::

      >> cd EMC_post/build
      >> make install

11. Assuming the modified code built successfully and you were able to produce Grib2
    output, you can check the Grib2 file for your new variable.

    GRIB2 output of the new variable from this example procedure (using the wgrib2 utility if
    available on your system).
     - The new variable will not be defined by the variable name. Instead it will be defined
       using the Grib2 parameter information you entered into params_grib2_tbl_new from
       step 7 of this procedure.

  ::

    wgrib2 -V GFSPRS.006

    716:37731711:vt=2019061506:500 m underground:6 hour fcst:var discipline=2 center=7 local_table=1 parmcat=0 parm=231:
        ndata=73728:undef=0:mean=278.383:min=215.47:max=302.4
        grid_template=40:winds(N/S):
        Gaussian grid: (384 x 192) units 1e-06 input WE:NS output WE:SN
        number of latitudes between pole-equator=96 #points=73728
        lat 89.284225 to -89.284225
        lon 0.000000 to 359.062500 by 0.937500

