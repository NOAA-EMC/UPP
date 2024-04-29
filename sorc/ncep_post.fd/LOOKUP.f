!> @file
!> @brief module: lookup_mod defines _____???
  module  lookup_mod
!
  implicit none
!
  integer,parameter :: ITB=076   &  <! _____
  ,JTB=134   &    <! _____
  ,ITBQ=152  &    <! _____
  ,JTBQ=440       <! _____
  
  
  real :: PL   &  !< Pressure Level
  ,THL  &      !< Temperature
  ,RDQ &       !< Relative Humidity
  ,RDTH  &     !< Mixing Ratio
  ,RDP &       !< Density
  ,RDTHE &     !< Potential Temperature___?
  ,PLQ   &     !< Pressure Level for Q
  ,RDPQ &      !< Density for Q
  ,RDTHEQ      !< Potential Temperature for Q

  real,dimension(JTB)  :: QS0,SQS &         !< Specific humidity and its storage array
  real,dimension(ITB)  :: THE0,STHE &       !< Potential temperature and its storage array
  real,dimension(ITBQ) :: THE0Q,STHEQ &     !< Pressure-temperature table___?
  real,dimension(ITB,JTB) :: PTBL &         !< Pressure-temperature table___?
  real,dimension(JTB,ITB) :: TTBL &         !< Temperature-pressure table___?
  real,dimension(JTBQ,ITBQ) :: TTBLQ &      !< Temperature-pressure table for Q___?
!
  end module lookup_mod
