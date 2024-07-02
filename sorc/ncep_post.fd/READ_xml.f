!> @file
!> @brief read_postxconfig() reads the post available field XML file and post control XML file. 
!> Each set of output fields going to one output file will be saved and processed later. 
!> In other words, post control file will be read in whole once. 
!> 
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!>   2012-01-27 | Jun Wang | INITIAL CODE
!>   2015-03-10 | Lin Gan  | Replace XML file with flat file implementation with parameter marshalling
!>   2016-07-08 | J. Carley | Clean up prints 

      SUBROUTINE READ_xml()

       use xml_perl_data,only: post_avblflds,paramset,read_postxconfig
       use grib2_module, only: num_pset
       use rqstfld_mod,only: num_post_afld,MXLVL,lvlsxml
       use CTLBLK_mod, only: me
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     DECLARE VARIABLES.
!     
!******************************************************************************
!     START READCNTRL_XML HERE.
!     
!> @brief Read post available field table

      call read_postxconfig()
      num_post_afld=size(paramset(1)%param)
      num_pset=size(paramset)

! LinGan below line removed because now we only read one flat file
!
!      if(size(post_avblflds%param)==0) then
!        call read_xml_file_post_t( 'post_avblflds.xml')
!        num_post_afld=size(post_avblflds%param)
!        allocate(lvlsxml(MXLVL,num_post_afld))
!      write(*,*)'in readxml, aft read post_avblflds.xml,num_post_afld=',num_post_afld
!      endif
!
!     READ post cntrl file
!      write(*,*)'in readxml,bf readxml,size(paramset)=',size(paramset)
!      if(size(paramset)==0) then
!        call read_xml_file_post_t( 'postcntrl.xml')
!        num_pset=size(paramset)
!        write(*,*)'in readxml, aft read postcntrl.xml,num_pset=',num_pset
!      endif
!

      RETURN
      end subroutine read_xml

