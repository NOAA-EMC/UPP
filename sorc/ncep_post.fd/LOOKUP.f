!> @file
!> @brief module: lookup_mod defines _____???
  module  lookup_mod
!
  implicit none
!
  integer,parameter :: &
  ITB (:,:) &      !< 076
  ,JTB (:,:) &     !< 134
  ,ITBQ (:,:) &    !< 152
  ,JTBQ (:,:) &    !< 440
  
  real :: &        
  PL (:,:) &       !< Pressure Level___?
  ,THL (:,:) &     !< Temperature___?
  ,RDQ (:,:) &     !< Relative Humidity___?
  ,RDTH (:,:) &    !< Mixing Ratio___?
  ,RDP (:,:) &     !< Density___?
  ,RDTHE (:,:) &   !< Potential Temperature___?
  ,PLQ (:,:) &     !< Pressure Level for Q___?
  ,RDPQ (:,:) &    !< Density for Q___?
  ,RDTHEQ (:,:) &  !< Potential Temperature for Q___?

  real,dimension(JTB)  :: QS0,SQS (:,:) &         !< Specific humidity and its storage array___?
  real,dimension(ITB)  :: THE0,STHE (:,:) &       !< Potential temperature and its storage array___?
  real,dimension(ITBQ) :: THE0Q,STHEQ (:,:) &     !< Pressure-temperature table___?
  real,dimension(ITB,JTB) :: PTBL (:,:) &         !< Pressure-temperature table___?
  real,dimension(JTB,ITB) :: TTBL (:,:) &         !< Temperature-pressure table___?
  real,dimension(JTBQ,ITBQ) :: TTBLQ (:,:) &      !< Temperature-pressure table for Q___?
!
  end module lookup_mod
