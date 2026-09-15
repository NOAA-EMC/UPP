    ! This is a test program for UPP.
    !
    ! This program tests routines from the SeverityMaps module.
    !
    ! Alyson Stahl, 9/2026
    program test_severity_maps
    use SeverityMaps
    implicit none

    real, parameter :: tol = 1.0e-5
    integer :: res
    
    res = 0
    
    call test_twp_map(res)
    if (res .ne. 0) stop 10

    call test_t_map(res)
    if (res .ne. 0) stop 11

    call test_prcpCondensate_map(res)
    if (res .ne. 0) stop 12

    call test_deltaZ_map(res)
    if (res .ne. 0) stop 13

    call test_ctt_map(res)
    if (res .ne. 0) stop 14

    call test_vv_map(res)
    if (res .ne. 0) stop 15

    call test_cldTopDist_map(res)
    if (res .ne. 0) stop 16

    call test_cldBaseDist_map(res)
    if (res .ne. 0) stop 17

    call test_deltaQ_map(res)
    if (res .ne. 0) stop 18

    call test_rh_map(res)
    if (res .ne. 0) stop 19

    call test_condensate_map(res)
    if (res .ne. 0) stop 20

    call test_convect_t_map(res)
    if (res .ne. 0) stop 21

    call test_convect_qpf_map(res)
    if (res .ne. 0) stop 22

    call test_convect_cape_map(res)
    if (res .ne. 0) stop 23

    call test_convect_liftedIdx_map(res)
    if (res .ne. 0) stop 24

    call test_convect_kIdx_map(res)
    if (res .ne. 0) stop 25

    call test_convect_totals_map(res)
    if (res .ne. 0) stop 26
    
    call test_moisture_map_cond(res)
    if (res .ne. 0) stop 27

    call test_moisture_map_cwat(res)
    if (res .ne. 0) stop 28

    print *, "SUCCESS!"

contains

    subroutine test_twp_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 16
        real :: v(ntests)
        integer :: scenario(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: Precipitation above warmnose scenario w/ v < 0
        v(1) = -1.0
        scenario(1) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(1) = 0.0

        ! Test Case 2: Precipitation above warmnose scenario w/ v >= 0
        v(2) = 0.0
        scenario(2) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(2) = 0.0

        ! Test Case 3: Precipitation above warmnose scenario w/ 0 < v < 1000
        v(3) = 500.0
        scenario(3) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(3) = 0.5

        ! Test Case 4: Precipitation above warmnose scenario w/ v = 1000.
        v(4) = 1000.0
        scenario(4) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(4) = 1.0

        ! Test Case 5: Precipitation above warmnose scenario w/ v > 1000.
        v(5) = 1500.0
        scenario(5) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(5) = 1.0

        ! Test Case 6: All snow scenario w/ 0 < v < 1000
        v(6) = 250.0
        scenario(6) = SCENARIOS%ALL_SNOW
        expected(6) = 0.25

        ! Test Case 7: Cold rain scenario w/ 0 < v < 1000
        v(7) = 750.0
        scenario(7) = SCENARIOS%COLD_RAIN
        expected(7) = 0.75

        ! Test Case 8: Freezing precipitation scenario w/ 0 < v < 1000
        v(8) = 600.0
        scenario(8) = SCENARIOS%FREEZING_PRECIPITAION
        expected(8) = 0.6

        ! Test Case 9: Warm precipitation scenario w/ v < 0
        v(9) = -1.0
        scenario(9) = SCENARIOS%WARM_PRECIPITAION
        expected(9) = 0.0

        ! Test Case 10: Warm precipitation scenario w/ v = 0
        v(10) = 0.0
        scenario(10) = SCENARIOS%WARM_PRECIPITAION
        expected(10) = 0.0

        ! Test Case 11: Warm precipitation scenario w/ 0 < v < 500
        v(11) = 250.0
        scenario(11) = SCENARIOS%WARM_PRECIPITAION
        expected(11) = 0.5

        ! Test Case 12: Warm precipitation scenario w/ v = 500
        v(12) = 500.0
        scenario(12) = SCENARIOS%WARM_PRECIPITAION
        expected(12) = 1.0

        ! Test Case 13: Warm precipitation scenario w/ v > 500
        v(13) = 750.0
        scenario(13) = SCENARIOS%WARM_PRECIPITAION
        expected(13) = 1.0

        ! Test Case 14: No precipitation scenario (should return 0 for all v)
        v(14) = 500.0
        scenario(14) = SCENARIOS%NO_PRECIPITAION
        expected(14) = 0.0

        ! Test Case 15: Precipitation below warmnose scenario (should return 0 for all v)
        v(15) = 500.0
        scenario(15) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(15) = 0.0

        ! Test Case 16: Convection scenario (should return 0 for all v)`
        v(16) = 500.0
        scenario(16) = SCENARIOS%CONVECTION
        expected(16) = 0.0

        do i = 1, ntests
        calculated = twp_map(v(i), scenario(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "twp_map() failed for test case ", i, ": expected ", expected(i), &
                " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_twp_map

    subroutine test_t_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 30
        real :: v(ntests)
        integer :: scenario(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: Precipitation below warmnose scenario w/ v < 269.15
        v(1) = 268.0
        scenario(1) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(1) = 1.0

        ! Test Case 2: Precipitation below warmnose scenario w/ v = 269.15
        v(2) = 269.15
        scenario(2) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(2) = 1.0
        
        ! Test Case 3: Precipitation below warmnose scenario w/ 269.15 < v < 273.15
        v(3) = 270.0
        scenario(3) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(3) = 0.7875

        ! Test Case 4: Precipitation below warmnose scenario w/ v = 273.15
        v(4) = 273.15
        scenario(4) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(4) = 0.0

        ! Test Case 5: Precipitation below warmnose scenario w/ v > 273.15
        v(5) = 274.0
        scenario(5) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(5) = 0.0

        ! NOTE: Only precip below warmnose should have a different temperature 
        ! map. The following test cases cover all other scenarios for one temp map.
        ! Test Case 6: No precipitation scenario w/ v < 248.15
        v(6) = 247.0
        scenario(6) = SCENARIOS%NO_PRECIPITAION
        expected(6) = 1.0

        ! Test Case 7: No precipitation scenario w/ v = 248.15
        v(7) = 248.15
        scenario(7) = SCENARIOS%NO_PRECIPITAION
        expected(7) = 1.0

        ! Test Case 8: No precipitation scenario w/ 248.15 < v < 248.65
        v(8) = 248.5
        scenario(8) = SCENARIOS%NO_PRECIPITAION
        expected(8) = 0.86273

        ! Test Case 9: No precipitation scenario w/ v = 248.65
        v(9) = 248.65
        scenario(9) = SCENARIOS%NO_PRECIPITAION
        expected(9) = 0.8039

        ! Test Case 10: No precipitation scenario w/ 248.65 < v < 249.15
        v(10) = 249.0
        scenario(10) = SCENARIOS%NO_PRECIPITAION
        expected(10) = 0.74699

        ! Test Case 11: No precipitation scenario w/ v = 249.15
        v(11) = 249.15
        scenario(11) = SCENARIOS%NO_PRECIPITAION
        expected(11) = 0.7226

        ! Test Case 12: No precipitation scenario w/ 249.15 < v < 251.15
        v(12) = 250.0
        scenario(12) = SCENARIOS%NO_PRECIPITAION
        expected(12) = 0.636325

        ! Test Case 13: No precipitation scenario w/ v = 251.15
        v(13) = 251.15
        scenario(13) = SCENARIOS%NO_PRECIPITAION
        expected(13) = 0.5196

        ! Test Case 14: No precipitation scenario w/ 251.15 < v < 253.15
        v(14) = 252.0
        scenario(14) = SCENARIOS%NO_PRECIPITAION
        expected(14) = 0.460185

        ! Test Case 15: No precipitation scenario w/ v = 253.15
        v(15) = 253.15
        scenario(15) = SCENARIOS%NO_PRECIPITAION
        expected(15) = 0.3798

        ! Test Case 16: No precipitation scenario w/ 253.15 < v < 255.15
        v(16) = 254.0
        scenario(16) = SCENARIOS%NO_PRECIPITAION
        expected(16) = 0.33152

        ! Test Case 17: No precipitation scenario w/ v = 255.15
        v(17) = 255.15
        scenario(17) = SCENARIOS%NO_PRECIPITAION
        expected(17) = 0.2662

        ! Test Case 18: No precipitation scenario w/ 255.15 < v < 257.15
        v(18) = 256.0
        scenario(18) = SCENARIOS%NO_PRECIPITAION
        expected(18) = 0.2244225

        ! Test Case 19: No precipitation scenario w/ v = 257.15
        v(19) = 257.15
        scenario(19) = SCENARIOS%NO_PRECIPITAION
        expected(19) = 0.1679

        ! Test Case 20: No precipitation scenario w/ 257.15 < v < 259.15
        v(20) = 258.0
        scenario(20) = SCENARIOS%NO_PRECIPITAION
        expected(20) = 0.130585

        ! Test Case 21: No precipitation scenario w/ v = 259.15
        v(21) = 259.15
        scenario(21) = SCENARIOS%NO_PRECIPITAION
        expected(21) = 0.0801

        ! Test Case 22: No precipitation scenario w/ 259.15 < v < 261.15
        v(22) = 260.0
        scenario(22) = SCENARIOS%NO_PRECIPITAION
        expected(22) = 0.0460575

        ! Test Case 23: No precipitation scenario w/ v = 261.15
        v(23) = 261.15
        scenario(23) = SCENARIOS%NO_PRECIPITAION
        expected(23) = 0.0

        ! Test Case 24: No precipitation scenario w/ v > 261.15
        v(24) = 262.0
        scenario(24) = SCENARIOS%NO_PRECIPITAION
        expected(24) = 0.0
        
        ! Test Case 25: Precipitation above warmnose scenario w/ v < 248.15
        v(25) = 247.0
        scenario(25) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(25) = 1.0

        ! Test Case 26: All snow scenario w/ v < 248.15
        v(26) = 247.0
        scenario(26) = SCENARIOS%ALL_SNOW
        expected(26) = 1.0

        ! Test Case 27: Cold rain scenario w/ v < 248.15
        v(27) = 247.0
        scenario(27) = SCENARIOS%COLD_RAIN
        expected(27) = 1.0

        ! Test Case 28: Warm precipitation scenario w/ v < 248.15
        v(28) = 247.0
        scenario(28) = SCENARIOS%WARM_PRECIPITAION
        expected(28) = 1.0

        ! Test Case 29: Freezing precipitation scenario w/ v < 248.15
        v(29) = 247.0
        scenario(29) = SCENARIOS%FREEZING_PRECIPITAION
        expected(29) = 1.0

        ! Test Case 30: Convection scenario w/ v < 248.15
        v(30) = 247.0
        scenario(30) = SCENARIOS%CONVECTION
        expected(30) = 1.0

        do i = 1, ntests
        calculated = t_map(v(i), scenario(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "t_map() failed for test case ", i, ": expected ", expected(i), &
                " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_t_map

    subroutine test_prcpCondensate_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 34
        real :: v(ntests)
        integer :: scenario(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: Precipitation below warmnose scenario w/ v < 0.05
        v(1) = 0.04
        scenario(1) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(1) = 0.0

        ! Test Case 2: Precipitation below warmnose scenario w/ v = 0.05
        v(2) = 0.05
        scenario(2) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(2) = 0.0

        ! Test Case 3: Precipitation below warmnose scenario w/ 0.05 < v < 0.2
        v(3) = 0.11
        scenario(3) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(3) = 0.4

        ! Test Case 4: Precipitation below warmnose scenario w/ v = 0.2
        v(4) = 0.2
        scenario(4) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(4) = 1.0

        ! Test Case 5: Precipitation below warmnose scenario w/ v > 0.2
        v(5) = 0.25
        scenario(5) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(5) = 1.0

        ! Test Case 6: Precipitation above warmnose scenario w/ v < 0.05
        v(6) = 0.04
        scenario(6) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(6) = 0.0

        ! Test Case 7: Precipitation above warmnose scenario w/ v = 0.05
        v(7) = 0.05
        scenario(7) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(7) = 0.0

        ! Test Case 8: Precipitation above warmnose scenario w/ 0.05 < v < 0.2
        v(8) = 0.11
        scenario(8) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(8) = 0.4

        ! Test Case 9: Precipitation above warmnose scenario w/ v = 0.2
        v(9) = 0.2
        scenario(9) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(9) = 1.0

        ! Test Case 10: Precipitation above warmnose scenario w/ v > 0.2
        v(10) = 0.25
        scenario(10) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(10) = 1.0

        ! Test Case 11: Cold rain scenario w/ v < 0.05
        v(11) = 0.04
        scenario(11) = SCENARIOS%COLD_RAIN
        expected(11) = 0.0

        ! Test Case 12: Cold rain scenario w/ v = 0.05
        v(12) = 0.05
        scenario(12) = SCENARIOS%COLD_RAIN
        expected(12) = 0.0

        ! Test Case 13: Cold rain scenario w/ 0.05 < v < 0.2
        v(13) = 0.11
        scenario(13) = SCENARIOS%COLD_RAIN
        expected(13) = 0.4

        ! Test Case 14: Cold rain scenario w/ v = 0.2
        v(14) = 0.2
        scenario(14) = SCENARIOS%COLD_RAIN
        expected(14) = 1.0

        ! Test Case 15: Cold rain scenario w/ v > 0.2
        v(15) = 0.25
        scenario(15) = SCENARIOS%COLD_RAIN
        expected(15) = 1.0

        ! Test Case 16: All snow scenario w/ v < 0.05
        v(16) = 0.04
        scenario(16) = SCENARIOS%ALL_SNOW
        expected(16) = 0.0

        ! Test Case 17: All snow scenario w/ v = 0.05
        v(17) = 0.05
        scenario(17) = SCENARIOS%ALL_SNOW
        expected(17) = 0.0

        ! Test Case 18: All snow scenario w/ 0.05 < v < 0.25
        v(18) = 0.1
        scenario(18) = SCENARIOS%ALL_SNOW
        expected(18) = 0.25

        ! Test Case 19: All snow scenario w/ v = 0.25
        v(19) = 0.25
        scenario(19) = SCENARIOS%ALL_SNOW
        expected(19) = 1.0

        ! Test Case 20: All snow scenario w/ v > 0.25
        v(20) = 0.30
        scenario(20) = SCENARIOS%ALL_SNOW
        expected(20) = 1.0

        ! Test Case 21: Warm precipitation scenario w/ v < 0.05
        v(21) = 0.04
        scenario(21) = SCENARIOS%WARM_PRECIPITAION
        expected(21) = 0.0

        ! Test Case 22: Warm precipitation scenario scenario w/ v = 0.05
        v(22) = 0.05
        scenario(22) = SCENARIOS%WARM_PRECIPITAION
        expected(22) = 0.0

        ! Test Case 23: Warm precipitation scenario w/ 0.05 < v < 0.15
        v(23) = 0.1
        scenario(23) = SCENARIOS%WARM_PRECIPITAION
        expected(23) = 0.25

        ! Test Case 24: Warm precipitation scenario w/ v = 0.15
        v(24) = 0.15
        scenario(24) = SCENARIOS%WARM_PRECIPITAION
        expected(24) = 0.5

        ! Test Case 25: Warm precipitation scenario w/ v > 0.15
        v(25) = 0.2
        scenario(25) = SCENARIOS%WARM_PRECIPITAION
        expected(25) = 0.5

        ! Test Case 26: Freezing precipitation scenario w/ v < 0.05
        v(26) = 0.04
        scenario(26) = SCENARIOS%FREEZING_PRECIPITAION
        expected(26) = 0.0

        ! Test Case 27: Freezing precipitation scenario w/ v = 0.05
        v(27) = 0.05
        scenario(27) = SCENARIOS%FREEZING_PRECIPITAION
        expected(27) = 0.0

        ! Test Case 28: Freezing precipitation scenario w/ 0.05 < v < 0.15
        v(28) = 0.1
        scenario(28) = SCENARIOS%FREEZING_PRECIPITAION
        expected(28) = 0.25

        ! Test Case 29: Freezing precipitation scenario w/ v = 0.15
        v(29) = 0.15
        scenario(29) = SCENARIOS%FREEZING_PRECIPITAION
        expected(29) = 0.5

        ! Test Case 30: Freezing precipitation scenario w/ v > 0.15
        v(30) = 0.2
        scenario(30) = SCENARIOS%FREEZING_PRECIPITAION
        expected(30) = 0.5

        ! Test Case 31: Convection scenario w/ v < 0.05 (Should return 0 for all v)
        v(31) = 0.04
        scenario(31) = SCENARIOS%CONVECTION
        expected(31) = 0.0

        ! Test Case 32: Convection scenario w/ v > 0.25 (Should return 0 for all v)
        v(32) = 0.3
        scenario(32) = SCENARIOS%CONVECTION
        expected(32) = 0.0

        ! Test Case 33: No precipitation scenario w/ v < 0.05 (Should return 0 for all v)
        v(33) = 0.04
        scenario(33) = SCENARIOS%NO_PRECIPITAION
        expected(33) = 0.0

        ! Test Case 34: No precipitation scenario w/ v > 0.25 (Should return 0 for all v)
        v(34) = 0.3
        scenario(34) = SCENARIOS%NO_PRECIPITAION
        expected(34) = 0.0

        do i = 1, ntests
        calculated = prcpCondensate_map(v(i), scenario(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "prcpCondensate_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_prcpCondensate_map

    subroutine test_deltaZ_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 30
        real :: v(ntests)
        integer :: scenario(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: No precipitation scenario w/ v < 0
        v(1) = -0.1
        scenario(1) = SCENARIOS%NO_PRECIPITAION
        expected(1) = 0.0

        ! Test Case 2: No precipitation scenario w/ v = 0
        v(2) = 0.0
        scenario(2) = SCENARIOS%NO_PRECIPITAION
        expected(2) = 0.0

        ! Test Case 3: No precipitation scenario w/ 0 < v < 1828.8
        v(3) = 457.2
        scenario(3) = SCENARIOS%NO_PRECIPITAION
        expected(3) = 0.25

        ! Test Case 4: No precipitation scenario w/ v = 1828.8
        v(4) = 1828.8
        scenario(4) = SCENARIOS%NO_PRECIPITAION
        expected(4) = 1.0

        ! Test Case 5: No precipitation scenario w/ v > 1828.8
        v(5) = 2000.0
        scenario(5) = SCENARIOS%NO_PRECIPITAION
        expected(5) = 1.0

        ! Test Case 6: Warm precipitation scenario w/ v < 0
        v(6) = -0.2
        scenario(6) = SCENARIOS%WARM_PRECIPITAION
        expected(6) = 0.0

        ! Test Case 7: Warm precipitation scenario w/ v = 0
        v(7) = 0.0
        scenario(7) = SCENARIOS%WARM_PRECIPITAION
        expected(7) = 0.0

        ! Test Case 8: Warm precipitation scenario w/ 0 < v < 1828.8
        v(8) = 914.4
        scenario(8) = SCENARIOS%WARM_PRECIPITAION
        expected(8) = 0.5

        ! Test Case 9: Warm precipitation scenario w/ v = 1828.8
        v(9) = 1828.8
        scenario(9) = SCENARIOS%WARM_PRECIPITAION
        expected(9) = 1.0

        ! Test Case 10: Warm precipitation scenario w/ v > 1828.8
        v(10) = 2200.0
        scenario(10) = SCENARIOS%WARM_PRECIPITAION
        expected(10) = 1.0

        ! Test Case 11: Precipitation below warmnose scenario w/ v < 30.5
        v(11) = 20.0
        scenario(11) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(11) = 0.0

        ! Test Case 12: Precipitation below warmnose scenario w/ v = 30.5
        v(12) = 30.5
        scenario(12) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(12) = 1.0

        ! Test Case 13: Precipitation below warmnose scenario w/ v > 30.5
        v(13) = 45.0
        scenario(13) = SCENARIOS%PRECIPITAION_BELOW_WARMNOSE
        expected(13) = 1.0

        ! Test Case 14: All snow scenario w/ v < 914.4
        v(14) = 900.0
        scenario(14) = SCENARIOS%ALL_SNOW
        expected(14) = 0.0

        ! Test Case 15: All snow scenario w/ v = 914.4
        v(15) = 914.4
        scenario(15) = SCENARIOS%ALL_SNOW
        expected(15) = 0.0

        ! Test Case 16: All snow scenario w/ 914.4 < v < 2743.2
        v(16) = 1371.6
        scenario(16) = SCENARIOS%ALL_SNOW
        expected(16) = 0.25

        ! Test Case 17: All snow scenario w/ v = 2743.2
        v(17) = 2743.2
        scenario(17) = SCENARIOS%ALL_SNOW
        expected(17) = 1.0

        ! Test Case 18: All snow scenario w/ v > 2743.2
        v(18) = 3000.0
        scenario(18) = SCENARIOS%ALL_SNOW
        expected(18) = 1.0

        ! Test Case 19: Cold rain scenario w/ v < 1524
        v(19) = 1523.0
        scenario(19) = SCENARIOS%COLD_RAIN
        expected(19) = 0.0

        ! Test Case 20: Cold rain scenario w/ v = 1524
        v(20) = 1524.0
        scenario(20) = SCENARIOS%COLD_RAIN
        expected(20) = 0.0

        ! Test Case 21: Cold rain scenario w/ 1524 < v < 4267.2
        v(21) = 2346.96
        scenario(21) = SCENARIOS%COLD_RAIN
        expected(21) = 0.3

        ! Test Case 22: Cold rain scenario w/ v = 4267.2
        v(22) = 4267.2
        scenario(22) = SCENARIOS%COLD_RAIN
        expected(22) = 1.0

        ! Test Case 23: Cold rain scenario w/ v > 4267.2
        v(23) = 5000.0
        scenario(23) = SCENARIOS%COLD_RAIN
        expected(23) = 1.0

        ! Test Case 24: Freezing precipitation scenario w/ v < 1524
        v(24) = 1500.0
        scenario(24) = SCENARIOS%FREEZING_PRECIPITAION
        expected(24) = 0.0

        ! Test Case 25: Freezing precipitation scenario w/ v = 1524
        v(25) = 1524.0
        scenario(25) = SCENARIOS%FREEZING_PRECIPITAION
        expected(25) = 0.0

        ! Test Case 26: Freezing precipitation scenario w/ 1524 < v < 4267.2
        v(26) = 3444.24
        scenario(26) = SCENARIOS%FREEZING_PRECIPITAION
        expected(26) = 0.7

        ! Test Case 27: Freezing precipitation scenario w/ v = 4267.2
        v(27) = 4267.2
        scenario(27) = SCENARIOS%FREEZING_PRECIPITAION
        expected(27) = 1.0

        ! Test Case 28: Freezing precipitation scenario w/ v > 4267.2
        v(28) = 4500.0
        scenario(28) = SCENARIOS%FREEZING_PRECIPITAION
        expected(28) = 1.0

        ! Test Case 29: Convection scenario (should return 0 for all v)
        v(29) = 5000.0
        scenario(29) = SCENARIOS%CONVECTION
        expected(29) = 0.0

        ! Test Case 30: Precipitation above warmnose scenario (should return 0 for all v)
        v(30) = 5000.0
        scenario(30) = SCENARIOS%PRECIPITAION_ABOVE_WARMNOSE
        expected(30) = 0.0
        
        do i = 1, ntests
        calculated = deltaZ_map(v(i), scenario(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "deltaZ_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_deltaZ_map

    subroutine test_ctt_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 11
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 223.15
        v(1) = 200.0
        expected(1) = 0.8

        ! Test Case 2: v = 223.15
        v(2) = 223.15
        expected(2) = 0.8

        ! Test Case 3: 223.15 < v < 233.15
        v(3) = 230.0
        expected(3) = 0.762051

        ! Test Case 4: v = 233.15
        v(4) = 233.15
        expected(4) = 0.7446

        ! Test Case 5: 233.15 < v < 243.15
        v(5) = 240.0
        expected(5) = 0.630753

        ! Test Case 6: v = 243.15
        v(6) = 243.15
        expected(6) = 0.5784

        ! Test Case 7: 243.15 < v < 253.15
        v(7) = 250.0
        expected(7) = 0.388655

        ! Test Case 8: v = 253.15
        v(8) = 253.15
        expected(8) = 0.3014

        ! Test Case 9: 253.15 < v < 261.15
        v(9) = 255.0
        expected(9) = 0.23170125

        ! Test Case 10: v = 261.15
        v(10) = 261.15
        expected(10) = 0.0

        ! Test Case 11: v > 261.15
        v(11) = 270.0
        expected(11) = 0.0

        do i = 1, ntests
        calculated = ctt_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "ctt_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_ctt_map

    subroutine test_vv_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < -0.5
        v(1) = -1.0
        expected(1) = 1.0

        ! Test Case 2: v = -0.5
        v(2) = -0.5
        expected(2) = 1.0

        ! Test Case 3: -0.5 < v < 0.0
        v(3) = -0.25
        expected(3) = 0.5

        ! Test Case 4: v = 0.0
        v(4) = 0.0
        expected(4) = 0.0

        ! Test Case 5: v > 0
        v(5) = 0.25
        expected(5) = 0.0

        do i = 1, ntests
        calculated = vv_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "vv_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_vv_map

    subroutine test_cldTopDist_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 609.6
        v(1) = 500.0
        expected(1) = 1.0

        ! Test Case 2: v = 609.6
        v(2) = 609.6
        expected(2) = 1.0

        ! Test Case 3: 609.6 < v < 3048.0
        v(3) = 2072.64
        expected(3) = 0.4

        ! Test Case 4: v = 3048.0
        v(4) = 3048.0
        expected(4) = 0.0

        ! Test Case 5: v > 3048.0
        v(5) = 4000.0
        expected(5) = 0.0

        do i = 1, ntests
        calculated = cldTopDist_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "cldTopDist_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_cldTopDist_map

    subroutine test_cldBaseDist_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 304.8
        v(1) = 200.0
        expected(1) = 1.0

        ! Test Case 2: v = 304.8
        v(2) = 304.8
        expected(2) = 1.0

        ! Test Case 3: 304.8 < v < 1524.0
        v(3) = 853.44
        expected(3) = 0.55

        ! Test Case 4: v = 1524.0
        v(4) = 1524.0
        expected(4) = 0.0

        ! Test Case 5: v > 1524.0
        v(5) = 2000.0
        expected(5) = 0.0

        do i = 1, ntests
        calculated = cldBaseDist_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "cldBaseDist_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_cldBaseDist_map

    subroutine test_deltaQ_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 0
        v(1) = -1.0
        expected(1) = 0.0

        ! Test Case 2: v = 0
        v(2) = 0.0
        expected(2) = 0.0

        ! Test Case 3: 0 < v < 1
        v(3) = 0.5
        expected(3) = 0.5

        ! Test Case 4: v = 1
        v(4) = 1.0
        expected(4) = 1.0

        ! Test Case 5: v > 1
        v(5) = 2.0
        expected(5) = 1.0

        do i = 1, ntests
        calculated = deltaQ_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "deltaQ_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_deltaQ_map

    subroutine test_rh_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 70
        v(1) = 50.0
        expected(1) = 0.0

        ! Test Case 2: v = 70
        v(2) = 70.0
        expected(2) = 0.0

        ! Test Case 3: 70 < v < 100
        v(3) = 85.0
        expected(3) = 0.5

        ! Test Case 4: v = 100
        v(4) = 100.0
        expected(4) = 1.0

        ! Test Case 5: v > 100
        v(5) = 150.0
        expected(5) = 1.0

        do i = 1, ntests
        calculated = rh_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "rh_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_rh_map

    subroutine test_condensate_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 0.004
        v(1) = 0.002
        expected(1) = 0.0

        ! Test Case 2: v = 0.004
        v(2) = 0.004
        expected(2) = 0.0

        ! Test Case 3: 0.004 < v < 0.2
        v(3) = 0.0824
        expected(3) = 0.4

        ! Test Case 4: v = 0.2
        v(4) = 0.2
        expected(4) = 1.0

        ! Test Case 5: v > 0.2
        v(5) = 0.3
        expected(5) = 1.0

        do i = 1, ntests
        calculated = condensate_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "condensate_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_condensate_map

    subroutine test_convect_t_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 15
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 243.15
        v(1) = 200.0
        expected(1) = 0.0

        ! Test Case 2: v = 243.15
        v(2) = 243.15
        expected(2) = 0.0

        ! Test Case 3: 243.15 < v < 265.15
        v(3) = 254.15
        expected(3) = 0.5

        ! Test Case 4: v = 265.15
        v(4) = 265.15
        expected(4) = 1.0

        ! Test Case 5: 265.15 < v < 269.15
        v(5) = 267.15
        expected(5) = 1.0

        ! Test Case 6: v = 269.15
        v(6) = 269.15
        expected(6) = 1.0

        ! Test Case 7: 269.15 < v < 270.15
        v(7) = 269.5
        expected(7) = 0.9545

        ! Test Case 8: v = 270.15
        v(8) = 270.15
        expected(8) = 0.87

        ! Test Case 9: 270.15 < v < 271.15
        v(9) = 271.0
        expected(9) = 0.734

        ! Test Case 10: v = 271.15
        v(10) = 271.15
        expected(10) = 0.71

        ! Test Case 11: 271.15 < v < 272.15
        v(11) = 272.0
        expected(11) = 0.5315

        ! Test Case 12: v = 272.15
        v(12) = 272.15
        expected(12) = 0.5

        ! Test Case 13: 272.15 < v < 273.15
        v(13) = 272.5
        expected(13) = 0.325

        ! Test Case 14: v = 273.15
        v(14) = 273.15
        expected(14) = 0.0

        ! Test Case 15: v > 273.15
        v(15) = 274.0
        expected(15) = 0.0

        do i = 1, ntests
        calculated = convect_t_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "convect_t_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_convect_t_map

    subroutine test_convect_qpf_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 1
        v(1) = 0.5
        expected(1) = 0.0

        ! Test Case 2: v = 1
        v(2) = 1.0
        expected(2) = 0.0

        ! Test Case 3: 1 < v < 3
        v(3) = 1.5
        expected(3) = 0.25

        ! Test Case 4: v = 3
        v(4) = 3.0
        expected(4) = 1.0

        ! Test Case 5: v > 3
        v(5) = 4.0
        expected(5) = 1.0

        do i = 1, ntests
        calculated = convect_qpf_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "convect_qpf_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_convect_qpf_map

    subroutine test_convect_cape_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 1000
        v(1) = 500
        expected(1) = 0.0

        ! Test Case 2: v = 1000
        v(2) = 1000
        expected(2) = 0.0

        ! Test Case 3: 1000 < v < 2500
        v(3) = 2200
        expected(3) = 0.8

        ! Test Case 4: v = 2500
        v(4) = 2500
        expected(4) = 1.0

        ! Test Case 5: v > 2500
        v(5) = 3000
        expected(5) = 1.0

        do i = 1, ntests
        calculated = convect_cape_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "convect_cape_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_convect_cape_map

    subroutine test_convect_liftedIdx_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < -10
        v(1) = -15
        expected(1) = 1.0

        ! Test Case 2: v = -10
        v(2) = -10
        expected(2) = 1.0

        ! Test Case 3: -10 < v < 0
        v(3) = -5
        expected(3) = 0.5

        ! Test Case 4: v = 0
        v(4) = 0
        expected(4) = 0.0

        ! Test Case 5: v > 0
        v(5) = 1
        expected(5) = 0.0

        do i = 1, ntests
        calculated = convect_liftedIdx_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "convect_liftedIdx_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_convect_liftedIdx_map

    subroutine test_convect_kIdx_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 20
        v(1) = 10
        expected(1) = 0.0

        ! Test Case 2: v = 20
        v(2) = 20
        expected(2) = 0.0

        ! Test Case 3: 20 < v < 40
        v(3) = 30
        expected(3) = 0.5

        ! Test Case 4: v = 40
        v(4) = 40
        expected(4) = 1.0

        ! Test Case 5: v > 40
        v(5) = 50
        expected(5) = 1.0

        do i = 1, ntests
        calculated = convect_kIdx_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "convect_kIdx_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_convect_kIdx_map

    subroutine test_convect_totals_map(res)
        integer, intent(out) :: res
        integer, parameter :: ntests = 5
        real :: v(ntests)
        real :: expected(ntests), calculated
        integer :: i

        ! Test Case 1: v < 20
        v(1) = 10
        expected(1) = 0.0

        ! Test Case 2: v = 20
        v(2) = 20
        expected(2) = 0.0

        ! Test Case 3: 20 < v < 55
        v(3) = 34
        expected(3) = 0.4

        ! Test Case 4: v = 55
        v(4) = 55
        expected(4) = 1.0

        ! Test Case 5: v > 55
        v(5) = 60
        expected(5) = 1.0

        do i = 1, ntests
        calculated = convect_totals_map(v(i))
        if (abs(calculated - expected(i)) > tol) then
            print *, "convect_totals_map() failed for test case ", i, ": expected ", &
                expected(i), " but got ", calculated
            res = 1
        end if
        end do
    end subroutine test_convect_totals_map

    subroutine test_moisture_map_cond(res)
        integer, intent(out) :: res
        real :: rh, liqCond, iceCond, pres, t
        real :: expected, calculated
        rh = 85.0
        liqCond = 0.102
        iceCond = 0.102
        pres = 77503.5
        t = 270.0
        expected = 0.625

        calculated = moisture_map_cond(rh, liqCond, iceCond, pres, t)
        if (abs(calculated - expected) > tol) then
        print *, "moisture_map_cond() failed: expected ", expected, " but got ", calculated
        res = 1
        end if
    end subroutine test_moisture_map_cond

    subroutine test_moisture_map_cwat(res)
        integer, intent(out) :: res
        real :: rh, cwat, pres, t
        real :: expected, calculated
        rh = 85.0
        cwat = 0.102
        pres = 77503.5
        t = 270.0
        expected = 0.625

        calculated = moisture_map_cwat(rh, cwat, pres, t)
        if (abs(calculated - expected) > tol) then
        print *, "moisture_map_cwat() failed: expected ", expected, " but got ", calculated
        res = 1
        end if
    end subroutine test_moisture_map_cwat

end program test_severity_maps