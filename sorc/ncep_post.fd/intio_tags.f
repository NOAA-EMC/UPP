!> @file
!> @brief module: intio_tags_mod defines variables related to integer
!input/output
module intio_tags_mod
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
  INTEGER, PARAMETER :: &
  int_ioexit  (:,:) &                   !< 10
  ,int_open_for_write_begin  (:,:) &	 	!< 20
  ,int_open_for_write_commit  (:,:) &	  !< 30
  ,int_open_for_read  (:,:) & 		 	    !< 40
  ,int_intio_nextrec  (:,:) & 		      !< 50
  ,int_inquire_opened  (:,:) & 		  	  !< 60
  ,int_inquire_filename  (:,:) & 		  	!< 70
  ,int_iosync  (:,:) & 			  	        !< 80
  ,int_ioclose  (:,:) & 			  	      !< 90
  ,int_next_time  (:,:) & 			  	    !< 100
  ,int_set_time  (:,:) & 			  	      !< 110
  ,int_next_var  (:,:) & 			  	      !< 120
  ,int_dom_ti_real  (:,:) & 		  	    !< 140
  ,int_dom_ti_double  (:,:) & 		  	  !< 160
  ,int_dom_ti_integer  (:,:) &		  	  !< 180
  ,int_dom_ti_logical  (:,:) & 		  	  !< 200
  ,int_dom_ti_char  (:,:) & 		  	    !< 220
  ,int_dom_td_real  (:,:) & 		  	    !< 240
  ,int_dom_td_double  (:,:) & 		  	  !< 260
  ,int_dom_td_integer  (:,:) & 		  	  !< 280
  ,int_dom_td_logical  (:,:) & 		  	  !< 300
  ,int_dom_td_char  (:,:) & 		 	      !< 320
  ,int_var_ti_real  (:,:) & 		  	    !< 340
  ,int_var_ti_double  (:,:) & 		  	  !< 360
  ,int_var_ti_integer  (:,:) &		  	  !< 380
  ,int_var_ti_logical  (:,:) & 		  	  !< 400
  ,int_var_ti_char  (:,:) & 		 	      !< 420
  ,int_var_td_real  (:,:) & 		  	    !< 440
  ,int_var_td_double  (:,:) & 		  	  !< 460
  ,int_var_td_integer  (:,:) & 		  	  !< 480
  ,int_var_td_logical  (:,:) & 		  	  !< 500
  ,int_var_td_char  (:,:) & 		  	    !< 520
  ,int_field  (:,:) & 			  	        !< 530
  ,int_var_info  (:,:) & 			  	      !< 540
  ,int_noop  (:,:) &			  	          !< 550

end module intio_tags_mod
