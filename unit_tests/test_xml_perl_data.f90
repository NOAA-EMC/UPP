! This is a test program for UPP.
!
! This program tests the subroutines in the xml_perl_data module.
!
! Alyson Stahl, 7/2026
program test_xml_perl_data
    use xml_perl_data
    use rqstfld_mod, only: num_post_afld
    use CTLBLK_mod, only: tprec, tclod, trdlw, trdsw, tsrfc, tmaxmin, td3d, filenameflat
    implicit none

    integer, parameter :: ntests = 3
    real, parameter :: tol = 1.0e-5
    !
    integer :: i, j, res
    character :: inpchar
    !
    integer :: EXP_NUM_POST_AFLD(ntests), EXP_PARAMSET_COUNT(ntests)
    real :: EXP_TPREC, EXP_DEFAULT
    type(paramset_t) :: EXP_PARAMSET(ntests)

    print *, "Testing filter_char_inp() subroutine..."

    ! Test Case 1: Input is "?", should return empty string.
    inpchar = "?"
    call filter_char_inp(inpchar)

    if (inpchar .ne. " ") then
        print *, "Test Case 1 Failed: Input is '?', Result: ", inpchar
        stop 10
    end if

    ! Test Case 2: Input should remain unchanged.
    inpchar = "a"
    call filter_char_inp(inpchar)
    if (inpchar .ne. "a") then
        print *, "Test Case 2 Failed: Input is 'a', Result: ", inpchar
        stop 20
    end if

    print *, "Testing read_postxconfig() subroutine..."
    
    ! These values should only be updated in Test Case 3.
    tprec = 2.
    tclod = 1
    trdlw = 1.
    trdsw = 1.
    tsrfc = 1.
    tmaxmin = 1.
    td3d = 1.

    EXP_TPREC = 2.0
    EXP_DEFAULT = 1.0

    EXP_NUM_POST_AFLD(1) = 0
    EXP_PARAMSET_COUNT(1) = 1

    EXP_NUM_POST_AFLD(2) = 2
    EXP_PARAMSET_COUNT(2) = 1

    EXP_NUM_POST_AFLD(3) = 1
    EXP_PARAMSET_COUNT(3) = 1

    call setup_expected_paramsets()

    ! Test Case 1: param_count = 0
    res = 0
    filenameflat = "data/ref_test_xml_perl_data_case1.txt"
    call read_postxconfig()

    if (num_post_afld .ne. EXP_NUM_POST_AFLD(1)) then
        print *, "Test Case 1 Failed for num_post_afld: Expected ", EXP_NUM_POST_AFLD(1), &
            " but got ", num_post_afld
        res = 1
    end if

    if (tprec .ne. EXP_TPREC) then
        print *, "Test Case 1 Failed for tprec: Expected ", EXP_TPREC, " but got ", tprec
        res = 1
    end if

    if (tclod .ne. EXP_DEFAULT .OR. trdlw .ne. EXP_DEFAULT .OR. trdsw .ne. EXP_DEFAULT .OR. &
        tsrfc .ne. EXP_DEFAULT .OR. tmaxmin .ne. EXP_DEFAULT .OR. td3d .ne. EXP_DEFAULT) then
        print *, "Test Case 1 Failed: tclod, trdlw, trdsw, tsrfc, tmaxmin, or td3d updated unexpectedly."
        res = 1
    end if

    if (size(paramset) .ne. EXP_PARAMSET_COUNT(1)) then
        print *, "Test Case 1 Failed for paramset_count: Expected ", EXP_PARAMSET_COUNT(1), &
            " but got ", size(paramset)
        res = 1
    end if

    ! We expect the paramset to only have one element for all cases. The index of EXP_PARAMSET corresponds 
    ! to the test case number.
    call compare_paramset(paramset(1), EXP_PARAMSET(1), res)

    if (res .ne. 0) stop 30

    ! Test Case 2: param_count > 0 & gen_proc_type != "ens_fcst"
    ! One parameter has cc = cv = level_array_count = level2_array_count = scale_array_count = 0
    filenameflat = "data/ref_test_xml_perl_data_case2.txt"
    call read_postxconfig()

    if (num_post_afld .ne. EXP_NUM_POST_AFLD(2)) then
        print *, "Test Case 2 Failed for num_post_afld: Expected ", EXP_NUM_POST_AFLD(2), &
            " but got ", num_post_afld
        res = 1
    end if

    if (tprec .ne. EXP_TPREC) then
        print *, "Test Case 2 Failed for tprec: Expected ", EXP_TPREC, " but got ", tprec
        res = 1
    end if

    if (tclod .ne. EXP_DEFAULT .OR. trdlw .ne. EXP_DEFAULT .OR. trdsw .ne. EXP_DEFAULT .OR. &
        tsrfc .ne. EXP_DEFAULT .OR. tmaxmin .ne. EXP_DEFAULT .OR. td3d .ne. EXP_DEFAULT) then
        print *, "Test Case 2 Failed: tclod, trdlw, trdsw, tsrfc, tmaxmin, or td3d updated unexpectedly."
        res = 1
    end if

    call compare_paramset(paramset(1), EXP_PARAMSET(2), res)

    if (res .ne. 0) stop 40
    
    ! Test Case 3: gen_proc_type = 'ens_fcst'
    filenameflat = "data/ref_test_xml_perl_data_case3.txt"
    call read_postxconfig()

    if (num_post_afld .ne. EXP_NUM_POST_AFLD(3)) then
        print *, "Test Case 3 Failed for num_post_afld: Expected ", EXP_NUM_POST_AFLD(3), &
            " but got ", num_post_afld
        res = 1
    end if

    if (tprec .ne. EXP_TPREC) then
        print *, "Test Case 3 Failed for tprec: Expected ", EXP_TPREC, " but got ", tprec
        res = 1
    end if

    if (tclod .ne. EXP_TPREC .OR. trdlw .ne. EXP_TPREC .OR. trdsw .ne. EXP_TPREC .OR. &
        tsrfc .ne. EXP_TPREC .OR. tmaxmin .ne. EXP_TPREC .OR. td3d .ne. EXP_TPREC) then
        print *, "Test Case 3 Failed: tclod, trdlw, trdsw, tsrfc, tmaxmin, or td3d were not updated correctly."
        res = 1
    end if

    call compare_paramset(paramset(1), EXP_PARAMSET(3), res)

    if (res .ne. 0) stop 50

    print *, "SUCCESS!"

contains

    subroutine compare_paramset(paramset, expected_paramset, res)
        implicit none
        type(paramset_t), intent(in) :: paramset
        type(paramset_t), intent(in) :: expected_paramset
        integer, intent(inout) :: res
        integer :: i, j

        if (paramset%datset .ne. expected_paramset%datset) then
            print *, "Test Failed for paramset%datset: Expected ", expected_paramset%datset, &
                " but got ", paramset%datset
            res = 1
        end if
        if (paramset%grid_num .ne. expected_paramset%grid_num) then
            print *, "Test Failed for paramset%grid_num: Expected ", expected_paramset%grid_num, &
                " but got ", paramset%grid_num
            res = 1
        end if
        if (paramset%sub_center .ne. expected_paramset%sub_center) then
            print *, "Test Failed for paramset%sub_center: Expected ", expected_paramset%sub_center, &
                " but got ", paramset%sub_center
            res = 1
        end if
        if (paramset%version_no .ne. expected_paramset%version_no) then
            print *, "Test Failed for paramset%version_no: Expected ", expected_paramset%version_no, &
                " but got ", paramset%version_no
            res = 1
        end if
        if (paramset%local_table_vers_no .ne. expected_paramset%local_table_vers_no) then
            print *, "Test Failed for paramset%local_table_vers_no: Expected ", &
                expected_paramset%local_table_vers_no, " but got ", paramset%local_table_vers_no
            res = 1
        end if
        if (paramset%sigreftime .ne. expected_paramset%sigreftime) then
            print *, "Test Failed for paramset%sigreftime: Expected ", expected_paramset%sigreftime, &
                " but got ", paramset%sigreftime
            res = 1
        end if
        if (paramset%prod_status .ne. expected_paramset%prod_status) then
            print *, "Test Failed for paramset%prod_status: Expected ", expected_paramset%prod_status, &
                " but got ", paramset%prod_status
            res = 1
        end if
        if (paramset%data_type .ne. expected_paramset%data_type) then
            print *, "Test Failed for paramset%data_type: Expected ", expected_paramset%data_type, &
                " but got ", paramset%data_type
            res = 1
        end if
        if (paramset%gen_proc_type .ne. expected_paramset%gen_proc_type) then
            print *, "Test Failed for paramset%gen_proc_type: Expected ", expected_paramset%gen_proc_type, &
                " but got ", paramset%gen_proc_type
            res = 1
        end if
        if (paramset%time_range_unit .ne. expected_paramset%time_range_unit) then
            print *, "Test Failed for paramset%time_range_unit: Expected ", &
                expected_paramset%time_range_unit, " but got ", paramset%time_range_unit
            res = 1
        end if
        if (paramset%orig_center .ne. expected_paramset%orig_center) then
            print *, "Test Failed for paramset%orig_center: Expected ", expected_paramset%orig_center, &
                " but got ", paramset%orig_center
            res = 1
        end if
        if (paramset%gen_proc .ne. expected_paramset%gen_proc) then
            print *, "Test Failed for paramset%gen_proc: Expected ", expected_paramset%gen_proc, &
                " but got ", paramset%gen_proc
            res = 1
        end if
        if (paramset%packing_method .ne. expected_paramset%packing_method) then
            print *, "Test Failed for paramset%packing_method: Expected ", &
                expected_paramset%packing_method, " but got ", paramset%packing_method
            res = 1
        end if
        if (paramset%order_of_sptdiff .ne. expected_paramset%order_of_sptdiff) then
            print *, "Test Failed for paramset%order_of_sptdiff: Expected ", &
                expected_paramset%order_of_sptdiff, " but got ", paramset%order_of_sptdiff
            res = 1
        end if
        if (paramset%field_datatype .ne. expected_paramset%field_datatype) then
            print *, "Test Failed for paramset%field_datatype: Expected ", &
                expected_paramset%field_datatype, " but got ", paramset%field_datatype
            res = 1
        end if
        if (paramset%comprs_type .ne. expected_paramset%comprs_type) then
            print *, "Test Failed for paramset%comprs_type: Expected ", expected_paramset%comprs_type, &
                " but got ", paramset%comprs_type
            res = 1
        end if
        if (paramset%type_ens_fcst .ne. expected_paramset%type_ens_fcst) then
            print *, "Test Failed for paramset%type_ens_fcst: Expected ", &
                expected_paramset%type_ens_fcst, " but got ", paramset%type_ens_fcst
            res = 1
        end if
        if (paramset%type_derived_fcst .ne. expected_paramset%type_derived_fcst) then
            print *, "Test Failed for paramset%type_derived_fcst: Expected ", &
                expected_paramset%type_derived_fcst, " but got ", paramset%type_derived_fcst
            res = 1
        end if

        if (.not. associated(paramset%param) .and. .not. associated(expected_paramset%param)) then
            return
        end if

        if (associated(paramset%param) .neqv. associated(expected_paramset%param)) then
            print *, "Test Failed for paramset%param association status."
            res = 1
            return
        end if

        if (size(paramset%param) .ne. size(expected_paramset%param)) then
            print *, "Test Failed for size(paramset%param): Expected ", size(expected_paramset%param), &
                " but got ", size(paramset%param)
            res = 1
            return
        end if

        do i = 1, size(paramset%param)
            if (paramset%param(i)%post_avblfldidx .ne. expected_paramset%param(i)%post_avblfldidx) then
                print *, "Test Failed for paramset%param(i)%post_avblfldidx at i=", i
                res = 1
            end if
            if (paramset%param(i)%shortname .ne. expected_paramset%param(i)%shortname) then
                print *, "Test Failed for paramset%param(i)%shortname at i=", i
                res = 1
            end if
            if (paramset%param(i)%longname .ne. expected_paramset%param(i)%longname) then
                print *, "Test Failed for paramset%param(i)%longname at i=", i
                res = 1
            end if
            if (paramset%param(i)%mass_windpoint .ne. expected_paramset%param(i)%mass_windpoint) then
                print *, "Test Failed for paramset%param(i)%mass_windpoint at i=", i
                res = 1
            end if
            if (paramset%param(i)%pdstmpl .ne. expected_paramset%param(i)%pdstmpl) then
                print *, "Test Failed for paramset%param(i)%pdstmpl at i=", i
                res = 1
            end if
            if (paramset%param(i)%pname .ne. expected_paramset%param(i)%pname) then
                print *, "Test Failed for paramset%param(i)%pname at i=", i
                res = 1
            end if
            if (paramset%param(i)%table_info .ne. expected_paramset%param(i)%table_info) then
                print *, "Test Failed for paramset%param(i)%table_info at i=", i
                res = 1
            end if
            if (paramset%param(i)%stats_proc .ne. expected_paramset%param(i)%stats_proc) then
                print *, "Test Failed for paramset%param(i)%stats_proc at i=", i
                res = 1
            end if
            if (paramset%param(i)%fixed_sfc1_type .ne. expected_paramset%param(i)%fixed_sfc1_type) then
                print *, "Test Failed for paramset%param(i)%fixed_sfc1_type at i=", i
                res = 1
            end if
            if (associated(paramset%param(i)%scale_fact_fixed_sfc1) .neqv. &
                associated(expected_paramset%param(i)%scale_fact_fixed_sfc1)) then
                print *, "Test Failed for scale_fact_fixed_sfc1 association at i=", i
                res = 1
            else if (associated(paramset%param(i)%scale_fact_fixed_sfc1)) then
                if (size(paramset%param(i)%scale_fact_fixed_sfc1) .ne. &
                    size(expected_paramset%param(i)%scale_fact_fixed_sfc1)) then
                    print *, "Test Failed for scale_fact_fixed_sfc1 size at i=", i
                    res = 1
                else
                    do j = 1, size(paramset%param(i)%scale_fact_fixed_sfc1)
                        if (paramset%param(i)%scale_fact_fixed_sfc1(j) .ne. &
                            expected_paramset%param(i)%scale_fact_fixed_sfc1(j)) then
                            print *, "Test Failed for scale_fact_fixed_sfc1 at i=", i, " j=", j
                            res = 1
                        end if
                    end do
                end if
            end if
            if (associated(paramset%param(i)%level) .neqv. associated(expected_paramset%param(i)%level)) then
                print *, "Test Failed for level association at i=", i
                res = 1
            else if (associated(paramset%param(i)%level)) then
                if (size(paramset%param(i)%level) .ne. size(expected_paramset%param(i)%level)) then
                    print *, "Test Failed for level size at i=", i
                    res = 1
                else
                    do j = 1, size(paramset%param(i)%level)
                        if (abs(paramset%param(i)%level(j) - expected_paramset%param(i)%level(j)) .gt. tol) then
                            print *, "Test Failed for level at i=", i, " j=", j
                            res = 1
                        end if
                    end do
                end if
            end if
            if (paramset%param(i)%fixed_sfc2_type .ne. expected_paramset%param(i)%fixed_sfc2_type) then
                print *, "Test Failed for paramset%param(i)%fixed_sfc2_type at i=", i
                res = 1
            end if
            if (associated(paramset%param(i)%scale_fact_fixed_sfc2) .neqv. &
                associated(expected_paramset%param(i)%scale_fact_fixed_sfc2)) then
                print *, "Test Failed for scale_fact_fixed_sfc2 association at i=", i
                res = 1
            else if (associated(paramset%param(i)%scale_fact_fixed_sfc2)) then
                if (size(paramset%param(i)%scale_fact_fixed_sfc2) .ne. &
                    size(expected_paramset%param(i)%scale_fact_fixed_sfc2)) then
                    print *, "Test Failed for scale_fact_fixed_sfc2 size at i=", i
                    res = 1
                else
                    do j = 1, size(paramset%param(i)%scale_fact_fixed_sfc2)
                        if (paramset%param(i)%scale_fact_fixed_sfc2(j) .ne. &
                            expected_paramset%param(i)%scale_fact_fixed_sfc2(j)) then
                            print *, "Test Failed for scale_fact_fixed_sfc2 at i=", i, " j=", j
                            res = 1
                        end if
                    end do
                end if
            end if
            if (associated(paramset%param(i)%level2) .neqv. associated(expected_paramset%param(i)%level2)) then
                print *, "Test Failed for level2 association at i=", i
                res = 1
            else if (associated(paramset%param(i)%level2)) then
                if (size(paramset%param(i)%level2) .ne. size(expected_paramset%param(i)%level2)) then
                    print *, "Test Failed for level2 size at i=", i
                    res = 1
                else
                    do j = 1, size(paramset%param(i)%level2)
                        if (abs(paramset%param(i)%level2(j) - expected_paramset%param(i)%level2(j)) .gt. tol) then
                            print *, "Test Failed for level2 at i=", i, " j=", j
                            res = 1
                        end if
                    end do
                end if
            end if
            if (paramset%param(i)%aerosol_type .ne. expected_paramset%param(i)%aerosol_type) then
                print *, "Test Failed for paramset%param(i)%aerosol_type at i=", i
                res = 1
            end if
            if (paramset%param(i)%prob_type .ne. expected_paramset%param(i)%prob_type) then
                print *, "Test Failed for paramset%param(i)%prob_type at i=", i
                res = 1
            end if
            if (paramset%param(i)%typ_intvl_size .ne. expected_paramset%param(i)%typ_intvl_size) then
                print *, "Test Failed for paramset%param(i)%typ_intvl_size at i=", i
                res = 1
            end if
            if (paramset%param(i)%scale_fact_1st_size .ne. expected_paramset%param(i)%scale_fact_1st_size) then
                print *, "Test Failed for paramset%param(i)%scale_fact_1st_size at i=", i
                res = 1
            end if
            if (abs(paramset%param(i)%scale_val_1st_size - expected_paramset%param(i)%scale_val_1st_size) .gt. tol) then
                print *, "Test Failed for paramset%param(i)%scale_val_1st_size at i=", i
                res = 1
            end if
            if (paramset%param(i)%scale_fact_2nd_size .ne. expected_paramset%param(i)%scale_fact_2nd_size) then
                print *, "Test Failed for paramset%param(i)%scale_fact_2nd_size at i=", i
                res = 1
            end if
            if (abs(paramset%param(i)%scale_val_2nd_size - expected_paramset%param(i)%scale_val_2nd_size) .gt. tol) then
                print *, "Test Failed for paramset%param(i)%scale_val_2nd_size at i=", i
                res = 1
            end if
            if (paramset%param(i)%typ_intvl_wvlen .ne. expected_paramset%param(i)%typ_intvl_wvlen) then
                print *, "Test Failed for paramset%param(i)%typ_intvl_wvlen at i=", i
                res = 1
            end if
            if (paramset%param(i)%scale_fact_1st_wvlen .ne. expected_paramset%param(i)%scale_fact_1st_wvlen) then
                print *, "Test Failed for paramset%param(i)%scale_fact_1st_wvlen at i=", i
                res = 1
            end if
            if (abs(paramset%param(i)%scale_val_1st_wvlen - expected_paramset%param(i)%scale_val_1st_wvlen) .gt. tol) then
                print *, "Test Failed for paramset%param(i)%scale_val_1st_wvlen at i=", i
                res = 1
            end if
            if (paramset%param(i)%scale_fact_2nd_wvlen .ne. expected_paramset%param(i)%scale_fact_2nd_wvlen) then
                print *, "Test Failed for paramset%param(i)%scale_fact_2nd_wvlen at i=", i
                res = 1
            end if
            if (abs(paramset%param(i)%scale_val_2nd_wvlen - expected_paramset%param(i)%scale_val_2nd_wvlen) .gt. tol) then
                print *, "Test Failed for paramset%param(i)%scale_val_2nd_wvlen at i=", i
                res = 1
            end if
            if (paramset%param(i)%scale_fact_lower_limit .ne. expected_paramset%param(i)%scale_fact_lower_limit) then
                print *, "Test Failed for paramset%param(i)%scale_fact_lower_limit at i=", i
                res = 1
            end if
            if (abs(paramset%param(i)%scale_val_lower_limit - expected_paramset%param(i)%scale_val_lower_limit) .gt. tol) then
                print *, "Test Failed for paramset%param(i)%scale_val_lower_limit at i=", i
                res = 1
            end if
            if (paramset%param(i)%scale_fact_upper_limit .ne. expected_paramset%param(i)%scale_fact_upper_limit) then
                print *, "Test Failed for paramset%param(i)%scale_fact_upper_limit at i=", i
                res = 1
            end if
            if (abs(paramset%param(i)%scale_val_upper_limit - expected_paramset%param(i)%scale_val_upper_limit) .gt. tol) then
                print *, "Test Failed for paramset%param(i)%scale_val_upper_limit at i=", i
                res = 1
            end if
            if (associated(paramset%param(i)%scale) .neqv. associated(expected_paramset%param(i)%scale)) then
                print *, "Test Failed for scale association at i=", i
                res = 1
            else if (associated(paramset%param(i)%scale)) then
                if (size(paramset%param(i)%scale) .ne. size(expected_paramset%param(i)%scale)) then
                    print *, "Test Failed for scale size at i=", i
                    res = 1
                else
                    do j = 1, size(paramset%param(i)%scale)
                        if (abs(paramset%param(i)%scale(j) - expected_paramset%param(i)%scale(j)) .gt. tol) then
                            print *, "Test Failed for scale at i=", i, " j=", j
                            res = 1
                        end if
                    end do
                end if
            end if
            if (paramset%param(i)%stat_miss_val .ne. expected_paramset%param(i)%stat_miss_val) then
                print *, "Test Failed for paramset%param(i)%stat_miss_val at i=", i
                res = 1
            end if
            if (paramset%param(i)%leng_time_range_prev .ne. expected_paramset%param(i)%leng_time_range_prev) then
                print *, "Test Failed for paramset%param(i)%leng_time_range_prev at i=", i
                res = 1
            end if
            if (paramset%param(i)%time_inc_betwn_succ_fld .ne. expected_paramset%param(i)%time_inc_betwn_succ_fld) then
                print *, "Test Failed for paramset%param(i)%time_inc_betwn_succ_fld at i=", i
                res = 1
            end if
            if (paramset%param(i)%type_of_time_inc .ne. expected_paramset%param(i)%type_of_time_inc) then
                print *, "Test Failed for paramset%param(i)%type_of_time_inc at i=", i
                res = 1
            end if
            if (paramset%param(i)%stat_unit_time_key_succ .ne. expected_paramset%param(i)%stat_unit_time_key_succ) then
                print *, "Test Failed for paramset%param(i)%stat_unit_time_key_succ at i=", i
                res = 1
            end if
            if (paramset%param(i)%bit_map_flag .ne. expected_paramset%param(i)%bit_map_flag) then
                print *, "Test Failed for paramset%param(i)%bit_map_flag at i=", i
                res = 1
            end if
        end do
    end subroutine compare_paramset

    subroutine setup_expected_paramsets()
        implicit none

        ! Test Case 1: param_count = 0 & gen_proc_type != "ens_fcst"
        EXP_PARAMSET(1)%datset = "ps0"
        EXP_PARAMSET(1)%grid_num = 1
        EXP_PARAMSET(1)%sub_center = "sc0"
        EXP_PARAMSET(1)%version_no = "v0"
        EXP_PARAMSET(1)%local_table_vers_no = "lt0"
        EXP_PARAMSET(1)%sigreftime = "sr0"
        EXP_PARAMSET(1)%prod_status = "ps0"
        EXP_PARAMSET(1)%data_type = "dt0"
        EXP_PARAMSET(1)%gen_proc_type = "fcst"
        EXP_PARAMSET(1)%time_range_unit = "hour"
        EXP_PARAMSET(1)%orig_center = "oc0"
        EXP_PARAMSET(1)%gen_proc = "gp0"
        EXP_PARAMSET(1)%packing_method = "pm0"
        EXP_PARAMSET(1)%field_datatype = "fd0"
        EXP_PARAMSET(1)%comprs_type = "ct0"
        allocate(EXP_PARAMSET(1)%param(0))

        ! Test Case 2: param_count > 0 & gen_proc_type != "ens_fcst"
        ! One parameter has cc = cv = level_array_count = level2_array_count = scale_array_count = 0
        EXP_PARAMSET(2)%datset = "ps1"
        EXP_PARAMSET(2)%grid_num = 2
        EXP_PARAMSET(2)%sub_center = "sc1"
        EXP_PARAMSET(2)%version_no = "v1"
        EXP_PARAMSET(2)%local_table_vers_no = "lt1"
        EXP_PARAMSET(2)%sigreftime = "sr1"
        EXP_PARAMSET(2)%prod_status = "ps1"
        EXP_PARAMSET(2)%data_type = "dt1"
        EXP_PARAMSET(2)%gen_proc_type = "fcst"
        EXP_PARAMSET(2)%time_range_unit = "hour"
        EXP_PARAMSET(2)%orig_center = "oc1"
        EXP_PARAMSET(2)%gen_proc = "gp1"
        EXP_PARAMSET(2)%packing_method = "pm1"
        EXP_PARAMSET(2)%field_datatype = "fd1"
        EXP_PARAMSET(2)%comprs_type = "ct1"
        allocate(EXP_PARAMSET(2)%param(2))
        EXP_PARAMSET(2)%param(1)%post_avblfldidx = 101
        EXP_PARAMSET(2)%param(1)%shortname = "sn1"
        EXP_PARAMSET(2)%param(1)%longname = "longname1"
        EXP_PARAMSET(2)%param(1)%pname = "pn1"
        EXP_PARAMSET(2)%param(1)%table_info = "tb1"
        EXP_PARAMSET(2)%param(1)%stats_proc = "sp1"
        EXP_PARAMSET(2)%param(1)%fixed_sfc1_type = "fs1"
        allocate(EXP_PARAMSET(2)%param(1)%scale_fact_fixed_sfc1(1))
        EXP_PARAMSET(2)%param(1)%scale_fact_fixed_sfc1 = 0
        allocate(EXP_PARAMSET(2)%param(1)%level(1))
        EXP_PARAMSET(2)%param(1)%level = 0.0
        EXP_PARAMSET(2)%param(1)%fixed_sfc2_type = "fs2"
        allocate(EXP_PARAMSET(2)%param(1)%scale_fact_fixed_sfc2(1))
        EXP_PARAMSET(2)%param(1)%scale_fact_fixed_sfc2 = 0
        EXP_PARAMSET(2)%param(1)%aerosol_type = "at1"
        EXP_PARAMSET(2)%param(1)%prob_type = "pt1"
        EXP_PARAMSET(2)%param(1)%typ_intvl_size = "tis1"
        EXP_PARAMSET(2)%param(1)%scale_fact_1st_size = 1
        EXP_PARAMSET(2)%param(1)%scale_val_1st_size = 1.0
        EXP_PARAMSET(2)%param(1)%scale_fact_2nd_size = 2
        EXP_PARAMSET(2)%param(1)%scale_val_2nd_size = 2.0
        EXP_PARAMSET(2)%param(1)%typ_intvl_wvlen = "tw1"
        EXP_PARAMSET(2)%param(1)%scale_fact_1st_wvlen = 3
        EXP_PARAMSET(2)%param(1)%scale_val_1st_wvlen = 3.0
        EXP_PARAMSET(2)%param(1)%scale_fact_2nd_wvlen = 4
        EXP_PARAMSET(2)%param(1)%scale_val_2nd_wvlen = 4.0
        EXP_PARAMSET(2)%param(1)%scale_fact_lower_limit = 5
        EXP_PARAMSET(2)%param(1)%scale_val_lower_limit = 5.0
        EXP_PARAMSET(2)%param(1)%scale_fact_upper_limit = 6
        EXP_PARAMSET(2)%param(1)%scale_val_upper_limit = 6.0
        allocate(EXP_PARAMSET(2)%param(1)%scale(1))
        EXP_PARAMSET(2)%param(1)%scale = 0.0
        EXP_PARAMSET(2)%param(1)%stat_miss_val = 7
        EXP_PARAMSET(2)%param(1)%leng_time_range_prev = 8
        EXP_PARAMSET(2)%param(1)%time_inc_betwn_succ_fld = 9
        EXP_PARAMSET(2)%param(1)%type_of_time_inc = "toi1"
        EXP_PARAMSET(2)%param(1)%stat_unit_time_key_succ = "su1"
        EXP_PARAMSET(2)%param(1)%bit_map_flag = "bm1"
        EXP_PARAMSET(2)%param(2)%post_avblfldidx = 102
        EXP_PARAMSET(2)%param(2)%shortname = "sn2"
        EXP_PARAMSET(2)%param(2)%longname = "longname2"
        EXP_PARAMSET(2)%param(2)%pname = "pn2"
        EXP_PARAMSET(2)%param(2)%table_info = "tb2"
        EXP_PARAMSET(2)%param(2)%stats_proc = "sp2"
        EXP_PARAMSET(2)%param(2)%fixed_sfc1_type = "fs1p"
        allocate(EXP_PARAMSET(2)%param(2)%scale_fact_fixed_sfc1(1))
        EXP_PARAMSET(2)%param(2)%scale_fact_fixed_sfc1 = 7
        allocate(EXP_PARAMSET(2)%param(2)%level(1))
        EXP_PARAMSET(2)%param(2)%level = 1000.0
        EXP_PARAMSET(2)%param(2)%fixed_sfc2_type = "fs2p"
        allocate(EXP_PARAMSET(2)%param(2)%scale_fact_fixed_sfc2(1))
        EXP_PARAMSET(2)%param(2)%scale_fact_fixed_sfc2 = 8
        allocate(EXP_PARAMSET(2)%param(2)%level2(1))
        EXP_PARAMSET(2)%param(2)%level2 = 2000.0
        EXP_PARAMSET(2)%param(2)%aerosol_type = "at2"
        EXP_PARAMSET(2)%param(2)%prob_type = "pt2"
        EXP_PARAMSET(2)%param(2)%typ_intvl_size = "tis2"
        EXP_PARAMSET(2)%param(2)%scale_fact_1st_size = 1
        EXP_PARAMSET(2)%param(2)%scale_val_1st_size = 1.5
        EXP_PARAMSET(2)%param(2)%scale_fact_2nd_size = 2
        EXP_PARAMSET(2)%param(2)%scale_val_2nd_size = 2.5
        EXP_PARAMSET(2)%param(2)%typ_intvl_wvlen = "tw2"
        EXP_PARAMSET(2)%param(2)%scale_fact_1st_wvlen = 3
        EXP_PARAMSET(2)%param(2)%scale_val_1st_wvlen = 3.5
        EXP_PARAMSET(2)%param(2)%scale_fact_2nd_wvlen = 4
        EXP_PARAMSET(2)%param(2)%scale_val_2nd_wvlen = 4.5
        EXP_PARAMSET(2)%param(2)%scale_fact_lower_limit = 5
        EXP_PARAMSET(2)%param(2)%scale_val_lower_limit = 5.5
        EXP_PARAMSET(2)%param(2)%scale_fact_upper_limit = 6
        EXP_PARAMSET(2)%param(2)%scale_val_upper_limit = 6.5
        allocate(EXP_PARAMSET(2)%param(2)%scale(1))
        EXP_PARAMSET(2)%param(2)%scale = 9.0
        EXP_PARAMSET(2)%param(2)%stat_miss_val = 10
        EXP_PARAMSET(2)%param(2)%leng_time_range_prev = 11
        EXP_PARAMSET(2)%param(2)%time_inc_betwn_succ_fld = 12
        EXP_PARAMSET(2)%param(2)%type_of_time_inc = "toi2"
        EXP_PARAMSET(2)%param(2)%stat_unit_time_key_succ = "su2"
        EXP_PARAMSET(2)%param(2)%bit_map_flag = "bm2"

        ! Test Case 3: gen_proc_type = 'ens_fcst'
        EXP_PARAMSET(3)%datset = "ps2"
        EXP_PARAMSET(3)%grid_num = 3
        EXP_PARAMSET(3)%sub_center = "sc2"
        EXP_PARAMSET(3)%version_no = "v2"
        EXP_PARAMSET(3)%local_table_vers_no = "lt2"
        EXP_PARAMSET(3)%sigreftime = "sr2"
        EXP_PARAMSET(3)%prod_status = "ps2"
        EXP_PARAMSET(3)%data_type = "dt2"
        EXP_PARAMSET(3)%gen_proc_type = "ens_fcst"
        EXP_PARAMSET(3)%time_range_unit = "hour"
        EXP_PARAMSET(3)%orig_center = "oc2"
        EXP_PARAMSET(3)%gen_proc = "gp2"
        EXP_PARAMSET(3)%packing_method = "pm2"
        EXP_PARAMSET(3)%field_datatype = "fd2"
        EXP_PARAMSET(3)%comprs_type = "ct2"
        EXP_PARAMSET(3)%type_ens_fcst = "ens1"
        allocate(EXP_PARAMSET(3)%param(1))
        EXP_PARAMSET(3)%param(1)%post_avblfldidx = 201
        EXP_PARAMSET(3)%param(1)%shortname = "sn2"
        EXP_PARAMSET(3)%param(1)%longname = "longname2"
        EXP_PARAMSET(3)%param(1)%pname = "pn2"
        EXP_PARAMSET(3)%param(1)%table_info = "tb2"
        EXP_PARAMSET(3)%param(1)%stats_proc = "sp2"
        EXP_PARAMSET(3)%param(1)%fixed_sfc1_type = "fs1p"
        allocate(EXP_PARAMSET(3)%param(1)%scale_fact_fixed_sfc1(1))
        EXP_PARAMSET(3)%param(1)%scale_fact_fixed_sfc1 = 7
        allocate(EXP_PARAMSET(3)%param(1)%level(1))
        EXP_PARAMSET(3)%param(1)%level = 1000.0
        EXP_PARAMSET(3)%param(1)%fixed_sfc2_type = "fs2p"
        allocate(EXP_PARAMSET(3)%param(1)%scale_fact_fixed_sfc2(1))
        EXP_PARAMSET(3)%param(1)%scale_fact_fixed_sfc2 = 8
        allocate(EXP_PARAMSET(3)%param(1)%level2(1))
        EXP_PARAMSET(3)%param(1)%level2 = 2000.0
        EXP_PARAMSET(3)%param(1)%aerosol_type = "at2"
        EXP_PARAMSET(3)%param(1)%prob_type = "pt2"
        EXP_PARAMSET(3)%param(1)%typ_intvl_size = "tis2"
        EXP_PARAMSET(3)%param(1)%scale_fact_1st_size = 1
        EXP_PARAMSET(3)%param(1)%scale_val_1st_size = 1.5
        EXP_PARAMSET(3)%param(1)%scale_fact_2nd_size = 2
        EXP_PARAMSET(3)%param(1)%scale_val_2nd_size = 2.5
        EXP_PARAMSET(3)%param(1)%typ_intvl_wvlen = "tw2"
        EXP_PARAMSET(3)%param(1)%scale_fact_1st_wvlen = 3
        EXP_PARAMSET(3)%param(1)%scale_val_1st_wvlen = 3.5
        EXP_PARAMSET(3)%param(1)%scale_fact_2nd_wvlen = 4
        EXP_PARAMSET(3)%param(1)%scale_val_2nd_wvlen = 4.5
        EXP_PARAMSET(3)%param(1)%scale_fact_lower_limit = 5
        EXP_PARAMSET(3)%param(1)%scale_val_lower_limit = 5.5
        EXP_PARAMSET(3)%param(1)%scale_fact_upper_limit = 6
        EXP_PARAMSET(3)%param(1)%scale_val_upper_limit = 6.5
        allocate(EXP_PARAMSET(3)%param(1)%scale(1))
        EXP_PARAMSET(3)%param(1)%scale = 9.0
        EXP_PARAMSET(3)%param(1)%stat_miss_val = 10
        EXP_PARAMSET(3)%param(1)%leng_time_range_prev = 11
        EXP_PARAMSET(3)%param(1)%time_inc_betwn_succ_fld = 12
        EXP_PARAMSET(3)%param(1)%type_of_time_inc = "toi2"
        EXP_PARAMSET(3)%param(1)%stat_unit_time_key_succ = "su2"
        EXP_PARAMSET(3)%param(1)%bit_map_flag = "bm2"

    end subroutine setup_expected_paramsets

end program test_xml_perl_data