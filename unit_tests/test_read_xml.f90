! This is a test program for UPP.
!
! This program tests the READ_xml() subroutine.
!
! Alyson Stahl, 7/2026

program test_read_xml
    use READ_XML_UPP_MOD
    use grib2_module, only: num_pset
    use rqstfld_mod, only: num_post_afld
    use CTLBLK_mod, only: filenameflat
    implicit none

    integer :: res
    integer :: EXP_NUM_PSET, EXP_NUM_POST_AFLD

    ! Also used in test_xml_perl_data.f90
    filenameflat = "data/ref_test_xml_perl_data_case2.txt"

    EXP_NUM_PSET = 1
    EXP_NUM_POST_AFLD = 2
    res = 0
    
    call READ_xml()

    if (num_pset .ne. EXP_NUM_PSET) then
        print *, "Test Failed: num_pset = ", num_pset, " Expected = ", EXP_NUM_PSET
        res = 1
    end if

    if (num_post_afld .ne. EXP_NUM_POST_AFLD) then
        print *, "Test Failed: num_post_afld = ", num_post_afld, " Expected = ", EXP_NUM_POST_AFLD
        res = 1
    end if

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_read_xml