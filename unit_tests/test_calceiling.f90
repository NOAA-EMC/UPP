! This is a test program for UPP.
!
! This program tests the CALCEILING() subroutine.
!
! Alyson Stahl, 4/2026
program test_calceiling
    use vrbls2d, only: fis
    use params_mod, only: small, gi
    use ctlblk_mod, only: jsta, jend, spval, im, modelname, ista, iend
    implicit none

    real, parameter :: tol = 1.0e-8
    integer, parameter :: npts = 10
    integer :: i, res
    real :: CLDZ(1, npts), TCLD(1, npts)
    real :: CEILING(1, npts), EXP_CEILING_RAPR(1, npts), EXP_CEILING(1, npts)
    
    interface
        subroutine CALCEILING(CLDZ,TCLD,CEILING)
            use ctlblk_mod, only: jsta, jend, ista, iend
            real, dimension(ista:iend,jsta:jend), intent(in) :: CLDZ, TCLD
            real, dimension(ista:iend,jsta:jend), intent(inout) :: CEILING 
        end subroutine CALCEILING
    end interface

    ! Grid parameters
    ista = 1
    iend = 1
    jsta = 1
    jend = npts
    spval = 9.9e10
    modelname = ""

    allocate(fis(1,1:npts))
    fis = 1500.0 / GI

    ! Test Case 1: TCLD > 50 and CLDZ(I,J) - FIS(I,J)*GI > 0 (general case)
    CLDZ = 2000.0
    TCLD = 60.0
    EXP_CEILING_RAPR = 500.0
    EXP_CEILING = 2000.0

    ! Test Case 2: TCLD == 50 and CLDZ(I,J) - FIS(I,J)*GI > 0 (boundary case)
    TCLD(1,2) = 50.0

    ! Test Case 3: TCLD > 50 and CLDZ(I,J) - FIS(I,J)*GI < 0
    CLDZ(1,3) = 1000.0
    EXP_CEILING_RAPR(1,3) = 20000.0
    EXP_CEILING(1,3) = 1000.0

    ! Test Case 4: TCLD < 50
    TCLD(1,4) = 40.0
    EXP_CEILING_RAPR(1,4) = 20000.0
    EXP_CEILING(1,4) = 20000.0

    ! Test Case 5: TCLD == spval
    TCLD(1,5) = spval
    EXP_CEILING_RAPR(1,5) = spval
    EXP_CEILING(1,5) = spval

    ! First testing with modelname != "RAPR"
    CEILING = -999.0
    call CALCEILING(CLDZ, TCLD, CEILING)

    res = 0
    do i = 1, npts
        if (abs(CEILING(1,i) - EXP_CEILING(1,i)) > tol) then
            print *, "Test Case ", i, " failed: Expected CEILING = ", EXP_CEILING(1,i), &
                     " but got ", CEILING(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) then 
        deallocate(fis)
        stop 10
    end if

    ! Now testing with modelname = "RAPR"
    modelname = "RAPR"
    CEILING = -999.0
    call CALCEILING(CLDZ, TCLD, CEILING)

    res = 0
    do i = 1, npts
        if (abs(CEILING(1,i) - EXP_CEILING_RAPR(1,i)) > tol) then
            print *, "Test Case ", i, " failed for modelname == RAPR: Expected CEILING = ", &
                EXP_CEILING_RAPR(1,i), " but got ", CEILING(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) then 
        deallocate(fis)
        stop 20
    end if

    deallocate(fis)
    print *, "SUCCESS!"
end program test_calceiling
