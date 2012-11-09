      SUBROUTINE READ_xml()
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    READCNTRLgrb2_xml  READS POST xml CONTROL FILE
!   PRGRMMR: J. WANG         ORG: NCEP/EMC   DATE: 12-01-27       
!     
! ABSTRACT:
!     THIS ROUTINE READS THE POST AVAILABLE FIELD XML FILE and 
!       POST CONTROL XML FILE. EACH SET OF OUTPUT FIELDS GOING TO ONE 
!       OUTPUT FILE WILL WILL BE SAVED AND PROCESSED LATER. IN OTHER
!       WORDS, POST CONTROL FILE WILL BE READ IN WHOLE ONCE.
!     
! PROGRAM HISTORY LOG:
!   01_27_2012  Jun Wang - INITIAL CODE
!     
! USAGE:    CALL READ_XML()
!   INPUT ARGUMENT LIST:
!     NONE
!
!   OUTPUT ARGUMENT LIST: 
!     NONE     - 
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!
!     LIBRARY:
!       COMMON   - RQSTFLDGRB2
!                  CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : IBM      
!$$$  
!
!     
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
       use xml_data_post_t,only: post_avblflds,paramset,read_xml_file_post_t
       use grib2_module, only: num_pset
       use rqstfld_mod,only: num_post_afld,MXLVL,lvlsxml
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     DECLARE VARIABLES.
!     
!******************************************************************************
!     START READCNTRL_XML HERE.
!     
!     READ post available field table
      print *,'in readxml,bf readxml,size(post_avblflds%param)=',size(post_avblflds%param)
!
      if(size(post_avblflds%param)==0) then
        call read_xml_file_post_t( 'post_avblflds.xml')
        num_post_afld=size(post_avblflds%param)
        allocate(lvlsxml(MXLVL,num_post_afld))
      print *,'in readxml, aft read post_avblflds.xml,num_post_afld=',num_post_afld
      endif
!
!     READ post cntrl file
      print *,'in readxml,bf readxml,size(paramset)=',size(paramset)
      if(size(paramset)==0) then
        call read_xml_file_post_t( 'postcntrl.xml')
        num_pset=size(paramset)
        print *,'in readxml, aft read postcntrl.xml,num_pset=',num_pset
      endif

      RETURN
      end subroutine read_xml

