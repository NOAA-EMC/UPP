!> @file
!> @brief wrf_io_flags declares variables related to WRF input/output.
    module wrf_io_flags_mod
      implicit none
      integer, parameter  :: WRF_FILE_NOT_OPENED                  = 100     !< Assigns ID 100 to WRF_FILE_NOT_OPENED
      integer, parameter  :: WRF_FILE_OPENED_NOT_COMMITTED        = 101     !< Assigns ID 101 to WRF_FILE_OPENED_NOT_COMMITTED
      integer, parameter  :: WRF_FILE_OPENED_AND_COMMITTED        = 102     !< Assigns ID 102 to WRF_FILE_OPENED_AND_COMMITTED
      integer, parameter  :: WRF_FILE_OPENED_FOR_READ             = 103     !< Assigns ID 103 to WRF_FILE_OPENED_FOR_READ
      integer, parameter  :: WRF_REAL                             = 104     !< Assigns ID 104 to WRF_REAL
      integer, parameter  :: WRF_REAL8                            = 105     !< Assigns ID 105 to WRF_REAL8
      integer, parameter  :: WRF_INTEGER                          = 106     !< Assigns ID 106 to WRF_INTEGER
      integer, parameter  :: WRF_LOGICAL                          = 107     !< Assigns ID 107 to WRF_LOGICAL
   end module wrf_io_flags_mod

