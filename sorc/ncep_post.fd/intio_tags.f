!> @file
!> @brief intio_tags_mod defines variables related to integer input/output

module intio_tags_mod

       implicit none
!
  INTEGER, PARAMETER ::  int_ioexit                        =10   !< Assigns ID 10 to int_ioexit - exit
  INTEGER, PARAMETER ::  int_open_for_write_begin          =20   !< Assigns ID 20 to int_open_for_write_begin - open for write operation
  INTEGER, PARAMETER ::  int_open_for_write_commit         =30   !< Assigns ID 30 to int_open_for_write_commit - open for write operation & commit
  INTEGER, PARAMETER ::  int_open_for_read                 =40   !< Assigns ID 40 to int_open_for_read - open for read
  INTEGER, PARAMETER ::  int_intio_nextrec                 =50   !< Assigns ID 50 to int_intio_nextrec - next record request
  INTEGER, PARAMETER ::  int_inquire_opened                =60   !< Assigns ID 60 to int_inquire_opened - check if opened
  INTEGER, PARAMETER ::  int_inquire_filename              =70   !< Assigns ID 70 to int_inquire_filename - filename request
  INTEGER, PARAMETER ::  int_iosync                        =80   !< Assigns ID 80 to int_iosync - sync data in buffer ?
  INTEGER, PARAMETER ::  int_ioclose                       =90   !< Assigns ID 90 to int_ioclose - close
  INTEGER, PARAMETER ::  int_next_time                     =100  !< Assigns ID 100 to int_next_time - get next timestamp
  INTEGER, PARAMETER ::  int_set_time                      =110  !< Assigns ID 110 to int_set_time - set time request
  INTEGER, PARAMETER ::  int_next_var                      =120  !< Assigns ID 120 to int_next_var - next variable request
  INTEGER, PARAMETER ::  int_dom_ti_real                   =140  !< Assigns ID 140 to int_dom_ti_real - time-independent domain metadata - real type
  INTEGER, PARAMETER ::  int_dom_ti_double                 =160  !< Assigns ID 160 to int_dom_ti_double - time-independent domain metadata - double type
  INTEGER, PARAMETER ::  int_dom_ti_integer                =180  !< Assigns ID 180 to int_dom_ti_integer - time-independent domain metadata - integer type
  INTEGER, PARAMETER ::  int_dom_ti_logical                =200  !< Assigns ID 200 to int_dom_ti_logical - time-independent domain metadata - logical/boolean type
  INTEGER, PARAMETER ::  int_dom_ti_char                   =220  !< Assigns ID 220 to int_dom_ti_char - time-independent domain metadata - character type
  INTEGER, PARAMETER ::  int_dom_td_real                   =240  !< Assigns ID 240 to int_dom_td_real - time-dependent domain metadata - real type
  INTEGER, PARAMETER ::  int_dom_td_double                 =260  !< Assigns ID 260 to int_dom_td_double - time-dependent domain metadata - double type
  INTEGER, PARAMETER ::  int_dom_td_integer                =280  !< Assigns ID 280 to int_dom_td_integer - time-dependent domain metadata - integer type
  INTEGER, PARAMETER ::  int_dom_td_logical                =300  !< Assigns ID 300 to int_dom_td_logical - time-dependent domain metadata - logical/boolean type
  INTEGER, PARAMETER ::  int_dom_td_char                   =320  !< Assigns ID 320 to int_dom_td_char - time-dependent domain metadata - character type
  INTEGER, PARAMETER ::  int_var_ti_real                   =340  !< Assigns ID 340 to int_var_ti_real - time-independent variable metadata - real type
  INTEGER, PARAMETER ::  int_var_ti_double                 =360  !< Assigns ID 360 to int_var_ti_double - time-independent variable metadata - double type
  INTEGER, PARAMETER ::  int_var_ti_integer                =380  !< Assigns ID 380 to int_var_ti_integer - time-independent variable metadata - integer type
  INTEGER, PARAMETER ::  int_var_ti_logical                =400  !< Assigns ID 400 to int_var_ti_logical - time-independent variable metadata - logical/boolean type
  INTEGER, PARAMETER ::  int_var_ti_char                   =420  !< Assigns ID 420 to int_var_ti_char - time-independent variable metadata - character type
  INTEGER, PARAMETER ::  int_var_td_real                   =440  !< Assigns ID 440 to int_var_td_real - time-dependent variable metadata - real type
  INTEGER, PARAMETER ::  int_var_td_double                 =460  !< Assigns ID 460 to int_var_td_double - time-dependent variable metadata - double type
  INTEGER, PARAMETER ::  int_var_td_integer                =480  !< Assigns ID 480 to int_var_td_integer - time-dependent variable metadata - integer type
  INTEGER, PARAMETER ::  int_var_td_logical                =500  !< Assigns ID 500 to int_var_td_logical - time-dependent variable metadata - logical/boolean type
  INTEGER, PARAMETER ::  int_var_td_char                   =520  !< Assigns ID 520 to int_var_td_char - time-dependent variable metadata - character type
  INTEGER, PARAMETER ::  int_field                         =530  !< Assigns ID 530 to int_field - write field request
  INTEGER, PARAMETER ::  int_var_info                      =540  !< Assigns ID 540 to int_var_info - variable info request ?
  INTEGER, PARAMETER ::  int_noop                          =550  !< Assigns ID 550 to int_noop - do nothing/no operation request

end module intio_tags_mod
