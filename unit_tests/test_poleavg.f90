! This is a test program for UPP.
!
! This program tests the POLEAVG() subroutine.
!
! Alyson Stahl, 2/2026
program test_poleavg
    implicit none

    real, parameter :: tol = 1.0e-8

    interface 
        subroutine POLEAVG(IM, JM, JSTA, JEND, SMALL, COSL, SPVAL, VAR)
            integer, intent(in) :: IM, JM, JSTA, JEND
            real, intent(in) :: SMALL, SPVAL
            real, dimension(IM,JSTA:JEND), intent(in) :: COSL
            real, dimension(IM,JSTA:JEND), intent(inout) :: VAR
        end subroutine POLEAVG
    end interface

    print *, "SUCCESS!"
end program test_poleavg