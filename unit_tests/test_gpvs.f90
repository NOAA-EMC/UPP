! This is a test program for UPP.
!
! This program tests the functions and subroutines in GPVS.f.
!
! The function FPVS0() is intentionally skipped because it
! no longer being used.
!
! Alyson Stahl, 1/2025
program test_gpvs
    use svptbl_mod, only: nx, c1xpvs, c2xpvs, c1xpvs0, c2xpvs0, tbpvs, tbpvs0
    implicit none

    
    print *, 'SUCCESS!'
end program test_gpvs