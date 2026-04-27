! This is a test program for UPP.
!
! This program tests the GETRECN() and ASSIGNNEMSIOVAR() subroutines in ASSIGNNEMSIOVAR.f.
!
! Alyson Stahl, 4/2026
program test_assignnemsiovar
    use ctlblk_mod, only: me
    implicit none

    real, parameter :: tol = 1e-8
    integer, parameter :: MAX_LEN = 50, NREC = 15, FLDSIZE = 16
    integer, parameter :: nx = 4, ny = 4
    integer :: i, j, k, res
    integer :: IM, JSTA, JEND, JSTA_2L, JEND_2U
    integer :: RECLEV(NREC), FLDLEV, OTHER_FLDLEV
    real :: SPVAL, TMP(FLDSIZE * NREC)
    character(MAX_LEN) :: RECNAME(NREC), RECLEVTYP(NREC)
    character(MAX_LEN) :: FLDNAME, FLDLEVTYP
    character(MAX_LEN) :: OTHER_FLDNAME, OTHER_FLDLEVTYP
    integer :: RECN, EXP_RECN
    real :: BUF(nx, ny), EXP_BUF(nx, ny)

    interface
        subroutine GETRECN(RECNAME, RECLEVTYP, RECLEV, NREC, FLDNAME, &
                            FLDLEVTYP, FLDLEV, RECN)
            integer, intent(in) :: NREC, FLDLEV
            integer, intent(in) :: RECLEV(NREC)
            character(*), intent(in) :: FLDNAME, FLDLEVTYP
            character(*), intent(in) :: RECNAME(NREC), RECLEVTYP(NREC)
            integer, intent(out) :: RECN
        end subroutine GETRECN
        subroutine ASSIGNNEMSIOVAR(IM, JSTA, JEND, JSTA_2L, JEND_2U, &
                                    L, NREC, FLDSIZE, SPVAL, TMP, RECNAME, &
                                    RECLEVTYP, RECLEV, VARNAME, VCOORDNAME, BUF)
            integer, intent(in) :: IM, JSTA, JEND, JSTA_2L, JEND_2U, L, NREC, FLDSIZE
            integer, intent(in) :: RECLEV(NREC)
            real, intent(in) :: SPVAL, TMP(FLDSIZE * NREC)
            character(*), intent(in) :: RECNAME(NREC), RECLEVTYP(NREC)
            character(*), intent(in) :: VARNAME, VCOORDNAME
            real, intent(out) :: BUF(IM, JSTA_2L:JEND_2U)
        end subroutine ASSIGNNEMSIOVAR
    end interface

    ! Running on a single process
    me = 0

    IM = nx
    JSTA = 1
    JEND = ny
    JSTA_2L = 1
    JEND_2U = ny
    SPVAL = -999.0

    ! Test Case 1: Standard case with a full match.

    ! Record to match
    FLDNAME = 'UGRD'
    FLDLEVTYP = 'mid layer'
    FLDLEV = 500

    ! Different record values
    OTHER_FLDNAME = 'VGRD'
    OTHER_FLDLEVTYP = 'bot layer'
    OTHER_FLDLEV = 250

    RECNAME(:) = OTHER_FLDNAME
    RECLEVTYP(:) = OTHER_FLDLEVTYP
    RECLEV(:) = OTHER_FLDLEV

    ! Condition 1: name match, levtyp match, lev no match
    RECNAME(1) = FLDNAME
    RECLEVTYP(1) = FLDLEVTYP
    
    ! Condition 2: name match, levtyp no match, lev match
    RECNAME(2) = FLDNAME
    RECLEV(2) = FLDLEV

    ! Condition 3: name match, levtyp no match, lev no match
    RECNAME(3) = FLDNAME

    ! Condition 4: name no match, levtyp match, lev match
    RECLEVTYP(4) = FLDLEVTYP
    RECLEV(4) = FLDLEV

    ! Condition 5: name no match, levtyp match, lev no match
    RECLEVTYP(5) = FLDLEVTYP

    ! Condition 6: name no match, levtyp no match, lev match
    RECLEV(6) = FLDLEV

    ! Full match at i = 9 (RECN = 9)
    RECNAME(9) = FLDNAME
    RECLEVTYP(9) = FLDLEVTYP
    RECLEV(9) = FLDLEV

    ! Duplicate match after i = 9 (verifies first-match return)
    RECNAME(12) = FLDNAME
    RECLEVTYP(12) = FLDLEVTYP
    RECLEV(12) = FLDLEV

    EXP_RECN = 9

    ! First testing that GETRECN() returns the correct record number for the first full match in the list of records.
    call GETRECN(RECNAME, RECLEVTYP, RECLEV, NREC, FLDNAME, FLDLEVTYP, FLDLEV, RECN)

    if (RECN .ne. EXP_RECN) then
        print *, "ERROR: GETRECN() returned ", RECN, " but expected ", EXP_RECN
        stop 10
    end if

    BUF = 0.0
    TMP = [(real(k), k=1, FLDSIZE*NREC)]

    do j = JSTA, JEND
        do i = 1, IM
            EXP_BUF(i,j) = TMP(i + (j-JSTA)*IM + (EXP_RECN-1)*FLDSIZE)
        end do
    end do

    call ASSIGNNEMSIOVAR(IM, JSTA, JEND, JSTA_2L, JEND_2U, FLDLEV, NREC, FLDSIZE, &
                        SPVAL, TMP, RECNAME, RECLEVTYP, RECLEV, FLDNAME, FLDLEVTYP, BUF)

    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(BUF(i,j) - EXP_BUF(i,j)) > tol) then
                print *, "ERROR: BUF(", i, ",", j, ") = ", BUF(i,j), &
                         " does not match expected value ", EXP_BUF(i,j)
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 20

    ! Test Case 2: No matching record exists.

    ! First checking that GETRECN() returns 0 when there is no match.
    FLDLEV = 999

    call GETRECN(RECNAME, RECLEVTYP, RECLEV, NREC, FLDNAME, FLDLEVTYP, FLDLEV, RECN)

    if (RECN .ne. 0) then
        print *, "ERROR: GETRECN() returned ", RECN, " but expected 0 when there is no match"
        stop 30
    end if

    ! Now checking that ASSIGNNEMSIOVAR() returns the SPVAL for all grid points when there is no match.
    call ASSIGNNEMSIOVAR(IM, JSTA, JEND, JSTA_2L, JEND_2U, FLDLEV, NREC, FLDSIZE, &
                        SPVAL, TMP, RECNAME, RECLEVTYP, RECLEV, FLDNAME, FLDLEVTYP, BUF)

    res = 0
    do i = 1, nx
        do j = 1, ny
            if (abs(BUF(i,j) - SPVAL) > tol) then
                print *, "ERROR: BUF(", i, ",", j, ") = ", BUF(i,j), &
                         " does not match expected value ", SPVAL, " when there is no match"
                res = 1
            end if
        end do
    end do

    if (res .ne. 0) stop 40

    print *, "SUCCESS!"
end program test_assignnemsiovar
