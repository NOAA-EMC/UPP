!> @file
!> @brief module: intio_tags_mod defines variables related to integer
!input/output
module intio_tags_mod
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
  INTEGER, PARAMETER :: 
  INTEGER, PARAMETER ::  int_ioexit			  	               =10   !< integer code for an I/O operation exit
  INTEGER, PARAMETER ::  int_open_for_write_begin	  	     =20   !< integer code for beginning an open operation for writing
  INTEGER, PARAMETER ::  int_open_for_write_commit	  	   =30   !< integer code for committing an open operation for writing
  INTEGER, PARAMETER ::  int_open_for_read 		  	         =40   !< integer code for an open operation for reading
  INTEGER, PARAMETER ::  int_intio_nextrec 		  	         =50   !< integer code for the next record in an I/O operation 
  INTEGER, PARAMETER ::  int_inquire_opened 		  	       =60   !< integer code for inquiring whether a file is opened 
  INTEGER, PARAMETER ::  int_inquire_filename 		  	     =70   !< integer code for inquiring about a filename 
  INTEGER, PARAMETER ::  int_iosync 			  	             =80   !< integer code for synchronizing I/O
  INTEGER, PARAMETER ::  int_ioclose 			 	               =90   !< integer code for closing an I/O operation
  INTEGER, PARAMETER ::  int_next_time 			  	           =100  !< integer code for getting the next time
  INTEGER, PARAMETER ::  int_set_time 			  	           =110  !< integer code for setting time 
  INTEGER, PARAMETER ::  int_next_var 			  	           =120  !< integer code for getting the next variable
  INTEGER, PARAMETER ::  int_dom_ti_real 		  	           =140  !< integer code for domain type information (real)
  INTEGER, PARAMETER ::  int_dom_ti_double 		  	         =160  !< integer code for domain type information (double)
  INTEGER, PARAMETER ::  int_dom_ti_integer 		  	       =180  !< integer code for domain type information (integer) 
  INTEGER, PARAMETER ::  int_dom_ti_logical 		 	         =200  !< integer code for domain type information (logical)
  INTEGER, PARAMETER ::  int_dom_ti_char 		 	             =220  !< integer code for domain type information (character)
  INTEGER, PARAMETER ::  int_dom_td_real 		  	           =240  !< integer code for domain type descriptors (real)
  INTEGER, PARAMETER ::  int_dom_td_double 		  	         =260  !< integer code for domain type descriptors (double) 
  INTEGER, PARAMETER ::  int_dom_td_integer 		  	       =280  !< integer code for domain type descriptors (integer) 
  INTEGER, PARAMETER ::  int_dom_td_logical 		  	       =300  !< integer code for domain type descriptors (logical)
  INTEGER, PARAMETER ::  int_dom_td_char 		 	             =320  !< integer code for domain type descriptors (character) 
  INTEGER, PARAMETER ::  int_var_ti_real 		  	           =340  !< integer code for variable information (real) 
  INTEGER, PARAMETER ::  int_var_ti_double 		  	         =360  !< integer code for variable information (double)
  INTEGER, PARAMETER ::  int_var_ti_integer 		 	         =380  !< integer code for variable information (integer)
  INTEGER, PARAMETER ::  int_var_ti_logical 		           =400  !< integer code for variable information (logical) 
  INTEGER, PARAMETER ::  int_var_ti_char 			             =420  !< integer code for variable information (character) 
  INTEGER, PARAMETER ::  int_var_td_real 			             =440  !< integer code for variable type descriptors (real) 
  INTEGER, PARAMETER ::  int_var_td_double 		 	           =460  !< integer code for variable type descriptors (double) 
  INTEGER, PARAMETER ::  int_var_td_integer 		 	         =480  !< integer code for variable type descriptors (integer) 
  INTEGER, PARAMETER ::  int_var_td_logical 		 	         =500  !< integer code for variable type descriptors (logical)
  INTEGER, PARAMETER ::  int_var_td_char 		 	             =520  !< integer code for variable type descriptors (character) 
  INTEGER, PARAMETER ::  int_field 			  	               =530  !< integer code for field-related operations
  INTEGER, PARAMETER ::  int_var_info 			  	           =540  !< integer code for variable information operations 
  INTEGER, PARAMETER ::  int_noop 			  	               =550  !< integer code for a no-operation or placeholder operation

end module intio_tags_mod
