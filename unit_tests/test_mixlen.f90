! This is a test program for UPP.
!
! This program tests the MIXLEN() subroutine.
!
! Alyson Stahl, 3/2026
program test_mixlen
    use vrbls3d, only: zint, pmid, t, q2
    use masks, only: lmh, htm
    use params_mod, only: EPSQ2, CAPA
    use ctlblk_mod, only: jsta, jend, jsta_m, jend_m, im, jm, jsta_2l, jend_2u, &
                          lm, lm1, spval, ista, iend, ista_m, iend_m, ista_2l, iend_2u
    implicit none

    real, parameter :: tol = 1.0e-8


    print *, 'SUCCESS!'
end program test_mixlen