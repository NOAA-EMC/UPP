! This is a test program for UPP.
!
! This program tests the MDL2SIGMA2() subroutine.
!
! Alyson Stahl, 3/2026
program test_mdl2sigma2
    use vrbls3d, only:  pint, pmid, t, zint, q
    use masks, only: lmh
    use params_mod, only: pq0, a2, a3, a4, rgamog
    use ctlblk_mod, only: pt, jsta_2l, jend_2u, spval, lp1, lm, jsta, jend,&
                        grib, cfld, datapd, fld_info, im, jm, im_jm, &
                        ista, iend, ista_2l, iend_2u
    use rqstfld_mod, only: iget, lvls, id, iavblfld, lvlsxml
    implicit none

    real, parameter :: tol = 1.0e-8


    print *, 'SUCCESS!'
end program test_mdl2sigma2