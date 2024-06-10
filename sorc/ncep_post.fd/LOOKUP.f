!> @file
!> @brief module: lookup_mod defines variables used to create lookup tables for pressure, temperature, and specific humidity.
  module  lookup_mod
!
  implicit none
!
  integer,parameter :: ITB=076,   &  !< Table horizontal size (i index)
  JTB=134,   &    !< Table vertical size (j index)
  ITBQ=152,  &    !< _____
  JTBQ=440       !< _____
  
  
  real :: PL   &  !< Lower bound of pressure range
  ,THL  &      !< Lower bound of potential temperature range
  ,RDQ &       !< Scaling factor for specific humidity
  ,RDTH  &     !< Scaling factor for potential temperature
  ,RDP &       !< Scaling factor for pressure
  ,RDTHE &     !< Scaling factor for equivalent potential temperature
  ,PLQ   &     !< _____ Lower bound of pressure range for specific humidity ?
  ,RDPQ &      !< _____ Scaling factor for pressure and specific humidity ?
  ,RDTHEQ      !< _____ Scaling factor for equivalent potential temperature and specific humidity ?

  real,dimension(JTB)  :: QS0, &    !< Base for specific humidity
  SQS                               !< Scaling factor for specific humidity
  real,dimension(ITB)  :: THE0, &   !< Base for equivalent potential temperature
  STHE                              !< Range for equivalent potential temperature
  real,dimension(ITBQ) :: THE0Q, &  !< _____
  STHEQ                             !< _____
  real,dimension(ITB,JTB) :: PTBL          !< Saturation pressure table
  real,dimension(JTB,ITB) :: TTBL          !< Temperature table
  real,dimension(JTBQ,ITBQ) :: TTBLQ       !< _____
!
  end module lookup_mod
