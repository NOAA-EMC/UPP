!> @file
!> @brief msfps() computes the map scale factor for a polar stereographic grid at a give latitude.
!>
!> This subroutine computes the map scale factor for a polar stereographic grid at a give latitude.
!>
!> @param[in] LAT Latitude at which map factor is valid.
!> @param[in] TRUELAT1 TRUELAT 1.
!> @param[out] MSF Map scale factor.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2006-11-01 | Rozumalski | Swiped from WRF si package
!>
!> @author Rozumalski @date 2006-11-01

! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!> msfps() computes the map scale factor for a polar stereographic grid at a give latitude.
!>
!> This subroutine computes the map scale factor for a polar stereographic grid at a give latitude.
!>
!> @param[in] LAT Latitude at which map factor is valid.
!> @param[in] TRUELAT1 TRUELAT 1.
!> @param[out] MSF Map scale factor.

      SUBROUTINE MSFPS(LAT,TRUELAT1,MSF)


! Computes the map scale factor for a Polar Stereographic grid at a given
! latitude.

      IMPLICIT NONE

! Define some private constants
!
      REAL, PARAMETER   :: pi = 3.1415927
      REAL, PARAMETER   :: rad_per_deg = pi / 180.

      REAL, INTENT(IN)           :: lat  ! latitude where msf is requested
      REAL, INTENT(IN)           :: truelat1
      REAL, INTENT(OUT)          :: msf

      REAL                       :: psi1, psix, pole

      IF (truelat1 >= 0.) THEN
        psi1 = (90. - truelat1) * rad_per_deg
        pole =90.
      ELSE
        psi1 = (90. + truelat1) * rad_per_deg
        pole = -90.
      ENDIF

      psix = (pole - lat)*rad_per_deg
      msf = ((1.+COS(psi1))/(1.0 + COS(psix)))
      RETURN

      END SUBROUTINE MSFPS

