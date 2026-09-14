! This is a test program for UPP.
!
! This program tests routines from the DerivedFields module.
!
! Alyson Stahl, 8/2026
program test_derived_fields
    use DerivedFields, only: derive_fields, PRECIPS
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: nz = 10
    integer :: i, j, k, res
    integer :: imp_physics, topoK
    real, dimension(nz) :: t, rh, pres, hgt, totalWater, totalCond
    real :: hprcp, hcprcp, cin, cape
    !
    real, dimension(nz) :: ept, wbt, twp
    real :: pc, kx, lx, tott
    integer :: prcpType
    !
    real, dimension(nz) :: exp_ept, exp_wbt, exp_twp
    real :: exp_pc, exp_kx, exp_lx, exp_tott
    integer :: exp_prcpType

    ! Test Case 1: Convection
    imp_physics = 99
    topoK = nz

    do k = 1, nz
        pres(k) = 5000.0 + 10000.0 * real(k)
        hgt(k) = 17000.0 - 1800.0 * real(k - 1)
        t(k) = 220.0 + 9.0 * real(k - 1)
        rh(k) = 40.0 + 6.0 * real(k - 1)
    end do

    totalWater = 0.0
    totalWater(3:5) = 0.002
    totalWater(8:10) = 0.002
    totalCond = 0.0

    hprcp = 0.10
    hcprcp = 1.50
    cin = -50.0
    cape = 1800.0

    exp_pc = 0.0
    exp_kx = 34.2969360
    exp_lx = -43.6572571
    exp_tott = 60.0737305
    exp_prcpType = PRECIPS%CONVECTION

    exp_ept(1) = 3.78657196E+02
    exp_ept(2) = 3.40894562E+02
    exp_ept(3) = 3.22361389E+02
    exp_ept(4) = 3.12390656E+02
    exp_ept(5) = 3.07577484E+02
    exp_ept(6) = 3.06790283E+02
    exp_ept(7) = 3.10039032E+02
    exp_ept(8) = 3.18258545E+02
    exp_ept(9) = 3.33539337E+02
    exp_ept(10) = 3.59922333E+02

    exp_wbt(1) = 2.19768265E+02
    exp_wbt(2) = 2.28567307E+02
    exp_wbt(3) = 2.37400467E+02
    exp_wbt(4) = 2.46090515E+02
    exp_wbt(5) = 2.54864166E+02
    exp_wbt(6) = 2.63451843E+02
    exp_wbt(7) = 2.72300842E+02
    exp_wbt(8) = 2.81353973E+02
    exp_wbt(9) = 2.90605042E+02
    exp_wbt(10) = 3.00208923E+02

    exp_twp = 0.0
    exp_twp(3:5) = 6.82362318E+00

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 1, res)

    if (res .ne. 0) stop 10

    ! Test Case 2: Precipitation threshold fails (Precipitation type is None)
    imp_physics = 11
    hcprcp = 0.0
    exp_prcpType = PRECIPS%NONE

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 2, res)

    if (res .ne. 0) stop 20

    ! Test Case 3: Eligible cloud not found (Precipitation type is None)
    totalCond = 0.0
    totalCond(8:10) = 0.004

    exp_pc = 1.20000001E-02
    exp_prcpType = PRECIPS%NONE

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 3, res)

    if (res .ne. 0) stop 30

    ! Test Case 4: Snow 
    rh(9) = 95.0
    t(9) = 260.0
    t(10) = 271.0

    exp_pc = 1.20000001E-02
    exp_kx = 34.2969360
    exp_lx = 25.7447357
    exp_tott = 60.0737305
    exp_prcpType = PRECIPS%SNOW

    exp_ept(9) =    267.728638    
    exp_ept(10) =    275.201019    
    exp_wbt(9) =    259.842285    
    exp_wbt(10) =    270.686951    
    exp_twp(8) =    9.44163132    
    exp_twp(9) =    9.44163132    
    exp_twp(10) =    9.44163132 

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 4, res)

    if (res .ne. 0) stop 40

    ! Test Case 5: Rain with tColdArea < 350.0
    rh(10) = 40.0
    t(10) = 274.0

    exp_pc = 1.20000001E-02
    exp_kx = 34.2969360
    exp_lx = 26.2238159
    exp_tott = 60.0737305
    exp_prcpType = PRECIPS%RAIN

    exp_ept(10) =     274.63641357
    exp_wbt(10) =     270.41192627
    exp_twp(8) =       4.58242607
    exp_twp(9) =       4.58242607
    exp_twp(10) =       4.58242607

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 5, res)

    if (res .ne. 0) stop 50

    ! Test Case 6: Rain with tColdArea >= 350.0
    rh(5) = 95.0
    rh(9) = 88.0
    rh(10) = 94.0
    t(9) = 292.0
    t(10) = 301.0
    
    exp_pc = 1.20000001E-02
    exp_kx = 34.2969360
    exp_lx = -43.6572571 
    exp_tott = 60.0737305
    exp_prcpType = PRECIPS%RAIN

    exp_ept(5) =     309.33511353
    exp_ept(9) =     333.53933716
    exp_ept(10) =     359.92233276
    exp_wbt(5) =     255.84771729
    exp_wbt(9) =     290.60504150
    exp_wbt(10) =     300.20892334
    exp_twp(8) =       0.00000000
    exp_twp(9) =       0.00000000
    exp_twp(10) =       0.00000000

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 6, res)

    if (res .ne. 0) stop 60

    ! Test Case 7: Other with coldTemp > 265.15
    rh(5) = 64.0
    rh(7) = 95.0

    exp_pc = 1.20000001E-02
    exp_kx = 35.8169861
    exp_lx = -43.6572571 
    exp_tott = 60.0737305
    exp_prcpType = PRECIPS%OTHER

    exp_ept(5) =     307.57748413
    exp_ept(7) =     313.02865601
    exp_wbt(5) =     254.86416626
    exp_wbt(7) =     273.64514160

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 7, res)

    if (res .ne. 0) stop 70

    ! Test Case 8: Other with coldTemp <= 265.15 & tColdArea >= 350.0 & wetBuldArea <= -250.0
    rh(5) = 95.0
    rh(7) = 76.0
    rh(10) = 40.0
    t(9) = 270.0
    t(10) = 268.0

    exp_pc = 1.20000001E-02
    exp_kx = 34.2969360
    exp_lx = 32.3780365 
    exp_tott = 60.0737305
    exp_prcpType = PRECIPS%OTHER

    exp_ept(5) =     309.33511353
    exp_wbt(5) =     255.84771729
    exp_ept(7) =     310.03903198
    exp_wbt(7) =     272.30084229
    exp_twp(8) =       9.32630539
    exp_ept(9) =     281.87905884
    exp_wbt(9) =     269.36096191
    exp_twp(9) =       9.32630539
    exp_ept(10) =     267.13070679
    exp_wbt(10) =     265.39135742
    exp_twp(10) =       9.32630539
    
    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 8, res)

    if (res .ne. 0) stop 80

    ! Test Case 9: Other with coldTemp <= 265.15 & tColdArea >= 350.0 & wetBuldArea > -250.0 & t(k) <= 273.15
    t(10) = 273.1
    rh(10) = 100.0

    exp_pc = 1.20000001E-02
    exp_kx = 34.2969360
    exp_lx = 22.4877319
    exp_tott = 60.0737305
    exp_prcpType = PRECIPS%OTHER

    exp_ept(10) = 279.16961670
    exp_wbt(10) = 273.10000610
    exp_twp(8) = 9.23454666
    exp_twp(9) = 9.23454666
    exp_twp(10) = 9.23454666

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 9, res)

    if (res .ne. 0) stop 90

    ! Test Case 10: abs(pres(k)- 50000.0) <= 0.1 & abs(pres(k)- 70000.0) <= 0.1) branches 
    ! in calc_indice()
    do k = 1, nz
        pres(k) = 10000.0 * real(k)
    end do

    exp_pc = 1.20000001E-02
    exp_kx = 17.78497314
    exp_lx = 23.35687256
    exp_tott = 38.68481445
    exp_prcpType = PRECIPS%OTHER

    exp_ept(1) =     425.37011719
    exp_ept(2) =     363.49542236
    exp_ept(3) =     337.06890869
    exp_ept(4) =     323.35174561
    exp_ept(5) =     318.45797729
    exp_ept(6) =     314.49588013
    exp_ept(7) =     317.14675903
    exp_ept(8) =     325.26940918
    exp_ept(9) =     286.71722412
    exp_ept(10) =     283.59588623
    exp_wbt(1) =     219.53652954
    exp_wbt(2) =     228.56730652
    exp_wbt(3) =     237.40046692
    exp_wbt(4) =     246.09051514
    exp_wbt(5) =     255.84771729
    exp_wbt(6) =     263.45184326
    exp_wbt(7) =     272.24224854
    exp_wbt(8) =     281.26251221
    exp_wbt(9) =     269.36096191
    exp_wbt(10) =     273.10000610
    exp_twp = 0.0
    exp_twp(3:5) = 6.06132603
    exp_twp(8:10) = 8.77268791

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 10, res)

    if (res .ne. 0) stop 100

    ! Test Case 11: abs(pres(k-1)- 50000.0) <= 0.1) & abs(pres(k-1)- 70000.0) <= 0.1) & 
    ! (abs(pres(k-1)- 85000.0) <= 0.1) branches in calc_indice()
    pres(5) = 49999.95
    pres(7) = 69999.95
    pres(8) = 84999.95
    pres(10) = 100000.0

    exp_pc = 1.20000001E-02
    exp_kx = 30.17388916
    exp_lx = 23.35687256
    exp_tott = 51.07373047
    exp_prcpType = PRECIPS%OTHER

    exp_ept(1) =     425.37011719
    exp_ept(2) =     363.49542236
    exp_ept(3) =     337.06890869
    exp_ept(4) =     323.35174561
    exp_ept(5) =     318.45806885
    exp_ept(6) =     314.49588013
    exp_ept(7) =     317.14682007
    exp_ept(8) =     318.25857544
    exp_ept(9) =     286.71722412
    exp_ept(10) =     283.59588623
    exp_wbt(1) =     219.53652954
    exp_wbt(2) =     228.56730652
    exp_wbt(3) =     237.40046692
    exp_wbt(4) =     246.09051514
    exp_wbt(5) =     255.84771729
    exp_wbt(6) =     263.45184326
    exp_wbt(7) =     272.24224854
    exp_wbt(8) =     281.35397339
    exp_wbt(9) =     269.36096191
    exp_wbt(10) =     273.10000610
    exp_twp = 0.0
    exp_twp(3:5) = 6.06132317
    exp_twp(8:10) = 8.77268791

    res = 0
    call derive_fields(imp_physics,t, rh, pres, hgt, totalWater, totalCond,&
                           nz, topoK, hprcp, hcprcp, cin, cape, &
                           ept, wbt, twp, pc, kx, lx, tott, prcpType)

    
    call check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                               exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                exp_ept, exp_wbt, exp_twp, 11, res)
                                
    if (res .ne. 0) stop 110

    print *, "SUCCESS!"

    contains

    subroutine check_expected_values(pc, kx, lx, tott, prcpType, ept, wbt, twp, &
                                        exp_pc, exp_kx, exp_lx, exp_tott, exp_prcpType, &
                                        exp_ept, exp_wbt, exp_twp, test_case, res)
        implicit none
        real, intent(in) :: pc, kx, lx, tott
        integer, intent(in) :: prcpType
        real, intent(in) :: ept(nz), wbt(nz), twp(:)
        real, intent(in) :: exp_pc, exp_kx, exp_lx, exp_tott
        integer, intent(in) :: exp_prcpType
        real, intent(in) :: exp_ept(nz), exp_wbt(nz), exp_twp(:)
        integer, intent(in) :: test_case
        integer, intent(inout) :: res
        integer :: i

    if (prcpType .ne. exp_prcpType) then
        print *, "Test Case ", test_case, " returned wrong precipitation type: ", prcpType
        res = 1
    end if

    if (abs(pc - exp_pc) > tol) then
        print *, "Test Case ", test_case, " expected pc = ", exp_pc, ", but got ", pc
        res = 1
    end if
    if (abs(kx - exp_kx) > tol) then
        print *, "Test Case ", test_case, " expected kx = ", exp_kx, ", but got ", kx
        res = 1
    end if
    if (abs(lx - exp_lx) > tol) then
        print *, "Test Case ", test_case, " expected lx = ", exp_lx, ", but got ", lx
        res = 1
    end if
    if (abs(tott - exp_tott) > tol) then
        print *, "Test Case ", test_case, " expected tott = ", exp_tott, ", but got ", tott
        res = 1
    end if

    do i = 1, nz
        if (abs(ept(i) - exp_ept(i)) > tol) then
            print '(A,I0,A,I0,A,F16.8,A,F16.8)', "Test Case ", test_case, " expected ept(", i, ") = ", &
                exp_ept(i), ", but got ", ept(i)
            res = 1
        end if
        if (abs(wbt(i) - exp_wbt(i)) > tol) then
            print '(A,I0,A,I0,A,F16.8,A,F16.8)', "Test Case ", test_case, " expected wbt(", i, ") = ", & 
                exp_wbt(i), ", but got ", wbt(i)
            res = 1
        end if
        if (abs(twp(i) - exp_twp(i)) > tol) then
            print '(A,I0,A,I0,A,F16.8,A,F16.8)', "Test Case ", test_case, " expected twp(", i, ") = ", & 
                exp_twp(i), ", but got ", twp(i)
            res = 1
        end if
    end do
    end subroutine check_expected_values

end program test_derived_fields