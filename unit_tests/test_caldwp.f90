! This is a test program for UPP.
!
! This program tests the CALDWP() subroutine.
!
! Alyson Stahl, 4/2026
program test_caldwp
    use params_mod, only: eps, oneps, d001, h1m12
    use ctlblk_mod, only: jsta, jend, spval, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 6
    integer :: i, j, res
    real :: P1D(1, npts), Q1D(1, npts), T1D(1, npts)
    real :: TDWP(1, npts), EXP_TDWP(1, npts)

    interface
        subroutine CALDWP(P1D,Q1D,TDWP,T1D)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in)    ::  P1D,Q1D,T1D
            real, dimension(ista:iend,jsta:jend), intent(inout) ::  TDWP
        end subroutine CALDWP
    end interface

    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    spval = 9.9e10

    ! Test Case 1: Standard case without any spvals or clipped values. Initialize input arrays with these values.
    P1D = 100000.0
    Q1D = 0.003
    T1D = 280.0
    EXP_TDWP(1,1) = 269.9227295

    ! Test Case 2: TDWP clipped to ambient temperature.
    Q1D(1, 2) = 0.99
    T1D(1, 2) = 250.0
    EXP_TDWP(1,2) = 250.0

    ! Test Case 3: MAX(H1M12,EVP(I,J)*D001) == H1M12
    P1D(1,3) = 0.0
    EXP_TDWP(1,3) = 192.1247864

    ! Test Case 4: P1D and Q1D have spvals.
    P1D(1, 4) = spval
    Q1D(1, 4) = spval
    EXP_TDWP(1,4) = 280.0

    ! Test Case 5: P1D has a spval.
    P1D(1, 5) = spval
    EXP_TDWP(1,5) = 280.0

    ! Test Case 6: Q1D has a spval.
    Q1D(1, 6) = spval
    EXP_TDWP(1,6) = 280.0
    
    call CALDWP(P1D, Q1D, TDWP, T1D)

    res = 0
    do i = 1, npts
        if (abs(TDWP(1,i) - EXP_TDWP(1,i)) > tol) then
            print *, "ERROR: TDWP(1,", i, ") = ", TDWP(1,i), &
                     " does not match expected value ", EXP_TDWP(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10
    
    print *, "SUCCESS!"

end program test_caldwp
