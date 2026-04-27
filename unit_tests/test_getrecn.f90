! This is a test program for UPP.
!
! This program tests the GETRECN() subroutine.
!
! Alyson Stahl, 2/2026
program test_getrecn
    use ctlblk_mod, only: me
    implicit none

    integer, parameter :: MAX_LEN = 50
    integer :: NREC = 15
    integer :: RECLEV(NREC), FLDLEV, OTHER_FLDLEV
    character(MAX_LEN) :: RECNAME(NREC), RECLEVTYP(NREC)
    character(MAX_LEN) :: FLDNAME, FLDLEVTYP
    character(MAX_LEN) :: OTHER_FLDNAME, OTHER_FLDLEVTYP
    integer :: RECN, EXP_RECN

    interface
        subroutine GETRECN(RECNAME, RECLEVTYP, RECLEV, NREC, FLDNAME, &
                            FLDLEVTYP, FLDLEV, RECN)
            integer, intent(in) :: NREC, FLDLEV
            integer, intent(in) :: RECLEV(NREC)
            character(*), intent(in) :: FLDNAME, FLDLEVTYP
            character(*), intent(in) :: RECNAME(NREC), RECLEVTYP(NREC)
            integer, intent(out) :: RECN
        end subroutine GETRECN
    end interface

    ! Running on a single process
    me = 0

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

    call GETRECN(RECNAME, RECLEVTYP, RECLEV, NREC, FLDNAME, FLDLEVTYP, FLDLEV, RECN)

    if (RECN .ne. EXP_RECN) then
        print *, "ERROR: GETRECN() returned ", RECN, " but expected ", EXP_RECN
        call exit(10)
    end if

    ! Test Case 2: No matching record exists.
    FLDLEV = 999

    call GETRECN(RECNAME, RECLEVTYP, RECLEV, NREC, FLDNAME, FLDLEVTYP, FLDLEV, RECN)

    if (RECN .ne. 0) then
        print *, "ERROR: GETRECN() returned ", RECN, " but expected 0 when there is no match"
        call exit(20)
    end if

    print *, "SUCCESS!"
end program test_getrecn
