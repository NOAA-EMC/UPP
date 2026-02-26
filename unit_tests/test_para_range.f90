! This is a test program for UPP.
!
! This program tests subroutines PARA_RANGE() and PARA_RANGE2() in PARA_RANGE.f.
!
! Alyson Stahl, 2/2026
program test_para_range
    implicit none

    real, parameter :: tol = 1.0e-8

    interface 
        subroutine PARA_RANGE(N1, N2, NPROCS, IRANK, ISTA, IEND)
            integer, intent(in) :: N1, N2, NPROCS, IRANK
            integer, intent(out) :: ISTA, IEND
        end subroutine PARA_RANGE
        subroutine PARA_RANGE2(IM, JM, NX, NY, NRANK, ISTA, IEND, JSTA, JEND)
            integer, intent(in) :: IM, JM, NX, NY, NRANK
            integer, intent(out) :: ISTA, IEND, JSTA, JEND
        end subroutine PARA_RANGE2
    end interface

    print *, "SUCCESS!"
end program test_para_range