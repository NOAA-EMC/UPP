SUBROUTINE wrf_debug( level , str ) 
  USE module_wrf_error 
  IMPLICIT NONE 
  CHARACTER*(*) str 
  INTEGER , INTENT (IN) :: level 
  CALL wrf_message( str ) 
  RETURN
END SUBROUTINE wrf_debug 
