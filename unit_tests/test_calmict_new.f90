! This is a test program for UPP.
!
! This program tests the CALMICT_new() subroutine.
!
! Alyson Stahl, 10/2026
program test_calmict_new
    use params_mod, only: dbzmin, epsq, tfrz, eps, rd, d608, oneps, nlimin, pi
    use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, ista, iend, ista_2l, iend_2u
    use rhgrd_mod, only: rhgrd
    use cmassi_mod, only: t_ice, rqr_drmin, n0rmin, cn0r_dmrmin, mdrmin, &
            rqr_drmax, cn0r_dmrmax, mdrmax, n0r0, xmrmin, xmrmax, &
            massi, cn0r0, mdimin, xmimax, mdimax, NLImax 
    implicit none

    real, parameter :: tol = 1.0e-6
    integer, parameter :: ntests = 15
    !
    integer :: i, res
    real, dimension(1, ntests) :: P1D, T1D, Q1D, C1D, FI1D, FR1D, FS1D, CUREFL
    real, dimension(1, ntests) :: QW1, QI1, QR1, QS1, DBZ1, DBZR1, DBZI1, DBZC1, NLICE1, NRAIN1
    real, dimension(1, ntests) :: EXP_QW1, EXP_QI1, EXP_QR1, EXP_QS1, EXP_DBZ1, EXP_DBZR1, & 
                                    EXP_DBZI1, EXP_DBZC1, EXP_NLICE1, EXP_NRAIN1

    interface
        subroutine CALMICT_new(P1D,T1D,Q1D,C1D,FI1D,FR1D,FS1D,CUREFL,   &
                                QW1,QI1,QR1,QS1,DBZ1,DBZR1,DBZI1,DBZC1, &
                                NLICE1,NRAIN1)
            use ctlblk_mod, only: jsta, jend, jsta_2l, jend_2u, ista, iend, &
                                    ista_2l, iend_2u
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(in)  ::  &
                                P1D, T1D, Q1D, C1D, FI1D, FR1D, FS1D, CUREFL
            real, dimension(ista_2l:iend_2u,jsta_2l:jend_2u), intent(inout) :: &
                                QW1, QI1, QR1, QS1, DBZ1, DBZR1, DBZI1, DBZC1, NLICE1, NRAIN1
        end subroutine CALMICT_new
    end interface

    ista = 1
    ista_2l = 1
    iend = 1
    iend_2u = 1
    jsta = 1
    jsta_2l = 1
    jend = ntests
    jend_2u = ntests

    rhgrd = 1.0
    NLImax = 5000.0
    xmrmax = 1000.0
    mdrmax = xmrmax

    t_ice = -40.0
    rqr_drmin = 1.54192e-7
    rqr_drmax = 2.48731e-2
    cn0r0 = 1.0e6 / sqrt(sqrt(pi * 1000.0 * n0r0))
    cn0r_dmrmin = 1.0 / (1000.0 * pi * (50.0e-6)**4)
    cn0r_dmrmax = 1.0 / (1000.0 * pi * (1000.0e-6)**4)
    
    do i = mdimin, mdimax
        massi(i) = 0.523599 * 100.0 * (real(i) * 1.0e-6)**3
    end do
    
    ! Test Case 1: 
    ! Default case
    ! Exercises following branches:
    ! - T1D - TFRZ > T_ICE and 0 < Fice < 1 
    ! - QW1 > 0 and 0 < Frain < 1
    ! - RQR_DRmin < RQR < RQR_DRmax
    ! - T1D - TFRZ < 0 and WVQW >= QSIgrd (sets FLARGE = 2)
    ! - XLI <= MASSI(MDImax)
    P1D(1,1) = 80000.0
    T1D(1,1) = 263.15
    Q1D(1,1) = 0.002
    C1D(1,1) = 1.0e-4
    FI1D(1,1) = 0.5
    FR1D(1,1) = 0.25
    FS1D(1,1) = 1.5
    CUREFL(1,1) = 0.0

    EXP_QW1(1,1) = 3.74999981E-05
    EXP_QI1(1,1) = 2.49608424E-11
    EXP_QR1(1,1) = 1.24999997E-05
    EXP_QS1(1,1) = 4.99999733E-05
    EXP_DBZ1(1,1) = 13.7838554    
    EXP_DBZR1(1,1) = 10.1326132    
    EXP_DBZI1(1,1) = 11.3319588    
    EXP_DBZC1(1,1) = -20.0  
    EXP_NLICE1(1,1) = 3363.78760    
    EXP_NRAIN1(1,1) = 1208.0 
    
    ! Test Case 2: QW1 > 0 and Frain = 1
    P1D(1,2) = 80000.0
    T1D(1,2) = 263.15
    Q1D(1,2) = 0.002
    C1D(1,2) = 1.0e-4
    FI1D(1,2) = 0.5
    FR1D(1,2) = 1.0
    FS1D(1,2) = 1.5
    CUREFL(1,2) = 0.0

    EXP_QW1(1,2) = 0.0    
    EXP_QI1(1,2) = 2.49608424E-11
    EXP_QR1(1,2) = 4.99999987E-05
    EXP_QS1(1,2) = 4.99999733E-05
    EXP_DBZ1(1,2) = 21.2050934    
    EXP_DBZR1(1,2) = 20.7331924    
    EXP_DBZI1(1,2) = 11.3319588    
    EXP_DBZC1(1,2) = -20.0    
    EXP_NLICE1(1,2) = 3363.78760    
    EXP_NRAIN1(1,2) = 1712.0

    ! Test Case 3: T1D - TFRZ > T_ICE and Fice = 1
    P1D(1,3) = 80000.0
    T1D(1,3) = 263.15
    Q1D(1,3) = 0.002
    C1D(1,3) = 1.0e-4
    FI1D(1,3) = 1.0
    FR1D(1,3) = 0.25
    FS1D(1,3) = 1.5
    CUREFL(1,3) = 0.0

    EXP_QW1(1,3) = 0.0    
    EXP_QI1(1,3) = 2.49608424E-11
    EXP_QR1(1,3) = 0.0   
    EXP_QS1(1,3) = 9.99999756E-05
    EXP_DBZ1(1,3) = 14.3422604    
    EXP_DBZR1(1,3) = -20.0   
    EXP_DBZI1(1,3) = 14.3422604    
    EXP_DBZC1(1,3) = -20.0    
    EXP_NLICE1(1,3) = 6727.57666    
    EXP_NRAIN1(1,3) = 0.0

    ! Test Case 4: T1D - TFRZ > T_ICE and Fice = 0
    P1D(1,4) = 80000.0
    T1D(1,4) = 263.15
    Q1D(1,4) = 0.002
    C1D(1,4) = 1.0e-4
    FI1D(1,4) = 0.0
    FR1D(1,4) = 0.25
    FS1D(1,4) = 1.5
    CUREFL(1,4) = 0.0

    EXP_QW1(1,4) = 7.49999963E-05
    EXP_QI1(1,4) = 0.0    
    EXP_QR1(1,4) = 2.49999994E-05
    EXP_QS1(1,4) = 0.0   
    EXP_DBZ1(1,4) = 15.4733009    
    EXP_DBZR1(1,4) = 15.4733009    
    EXP_DBZI1(1,4) = -20.0    
    EXP_DBZC1(1,4) = -20.0   
    EXP_NLICE1(1,4) = 0.0   
    EXP_NRAIN1(1,4) = 1440.0 

    ! Test Case 5: T1D - TFRZ >= 0 
    P1D(1,5) = 80000.0
    T1D(1,5) = 293.15
    Q1D(1,5) = 0.002
    C1D(1,5) = 1.0e-4
    FI1D(1,5) = 0.5
    FR1D(1,5) = 0.25
    FS1D(1,5) = 1.5
    CUREFL(1,5) = 0.0

    EXP_QW1(1,5) = 3.74999981E-05
    EXP_QI1(1,5) = 0.0 
    EXP_QR1(1,5) = 1.24999997E-05
    EXP_QS1(1,5) = 4.99999987E-05
    EXP_DBZ1(1,5) = 21.1522789    
    EXP_DBZR1(1,5) = 9.31643867    
    EXP_DBZI1(1,5) = 20.8579483    
    EXP_DBZC1(1,5) = -20.0    
    EXP_NLICE1(1,5) = 1599.53943    
    EXP_NRAIN1(1,5) = 1176.0 

    ! Test Case 6: RQLICE < NLI_min*MASSI(MDImin)  (min large mixing ratio)
    P1D(1,6) = 80000.0
    T1D(1,6) = 293.15
    Q1D(1,6) = 0.002
    C1D(1,6) = 5.0e-9
    FI1D(1,6) = 1.0
    FR1D(1,6) = 0.0
    FS1D(1,6) = 1.0
    CUREFL(1,6) = 0.0

    EXP_QW1(1,6) = 0.0    
    EXP_QI1(1,6) = 0.0   
    EXP_QR1(1,6) = 0.0  
    EXP_QS1(1,6) = 4.99999997E-09
    EXP_DBZ1(1,6) = -20.0   
    EXP_DBZR1(1,6) = -20.0    
    EXP_DBZI1(1,6) = -20.0    
    EXP_DBZC1(1,6) = -20.0    
    EXP_NLICE1(1,6) = 725.422668    
    EXP_NRAIN1(1,6) = 0.0

    ! Test Case 7: XLI <= MASSI(MDImin)
    P1D(1,7) = 80000.0
    T1D(1,7) = 293.15
    Q1D(1,7) = 0.002
    C1D(1,7) = 8.0e-9
    FI1D(1,7) = 1.0
    FR1D(1,7) = 0.0
    FS1D(1,7) = 1.5
    CUREFL(1,7) = 0.0

    EXP_QW1(1,7) = 0.0   
    EXP_QI1(1,7) = 0.0    
    EXP_QR1(1,7) = 0.0    
    EXP_QS1(1,7) = 7.99999977E-09
    EXP_DBZ1(1,7) = -20.0   
    EXP_DBZR1(1,7) = -20.0    
    EXP_DBZI1(1,7) = -20.0    
    EXP_DBZC1(1,7) = -20.0    
    EXP_NLICE1(1,7) = 773.784119    
    EXP_NRAIN1(1,7) = 0.0    

    ! Test Case 8: MASSI(MDIMin) < XLI <= MASSI(450)
    P1D(1,8) = 80000.0
    T1D(1,8) = 253.15
    Q1D(1,8) = 0.0001
    C1D(1,8) = 7.0e-5
    FI1D(1,8) = 1.0
    FR1D(1,8) = 0.0
    FS1D(1,8) = 1.5
    CUREFL(1,8) = 0.0
    
    EXP_QW1(1,8) = 0.0    
    EXP_QI1(1,8) = 0.0    
    EXP_QR1(1,8) = 0.0    
    EXP_QS1(1,8) = 6.99903248E-05
    EXP_DBZ1(1,8) = 2.89214897  
    EXP_DBZR1(1,8) = -20.0    
    EXP_DBZI1(1,8) = 2.89214897  
    EXP_DBZC1(1,8) = -20.0    
    EXP_NLICE1(1,8) = 49842.5117  
    EXP_NRAIN1(1,8) = 0.0  
 
    ! Test Case 9: RQR > RQmix and RQLICE > RQmix and RQR <= RQLICE
    P1D(1,9) = 80000.0
    T1D(1,9) = 263.15
    Q1D(1,9) = 0.002
    C1D(1,9) = 1.0e-2
    FI1D(1,9) = 0.5
    FR1D(1,9) = 0.5
    FS1D(1,9) = 1.5
    CUREFL(1,9) = 0.0

    EXP_QW1(1,9) = 2.49999994E-03
    EXP_QI1(1,9) = 2.49608424E-11
    EXP_QR1(1,9) = 2.49999994E-03
    EXP_QS1(1,9) = 4.99999989E-03
    EXP_DBZ1(1,9) = 50.7193451    
    EXP_DBZR1(1,9) = 50.4620819    
    EXP_DBZI1(1,9) = 38.3172874    
    EXP_DBZC1(1,9) = -20.0000000    
    EXP_NLICE1(1,9) = 67343.5859    
    EXP_NRAIN1(1,9) = 4552.00000    

    ! Test Case 10: RQR > RQmix and RQLICE > RQmix and RQR <= RQLICE
    P1D(1,10) = 80000.0
    T1D(1,10) = 263.15
    Q1D(1,10) = 0.002
    C1D(1,10) = 1.0e-1
    FI1D(1,10) = 0.95
    FR1D(1,10) = 0.5
    FS1D(1,10) = 1.5
    CUREFL(1,10) = 0.0

    EXP_QW1(1,10) = 2.50000134E-03
    EXP_QI1(1,10) = 2.49608424E-11
    EXP_QR1(1,10) = 2.50000134E-03
    EXP_QS1(1,10) = 9.49999988E-02
    EXP_DBZ1(1,10) = 57.8279572    
    EXP_DBZR1(1,10) = 50.4620819    
    EXP_DBZI1(1,10) = 56.9480286    
    EXP_DBZC1(1,10) = -20.0000000    
    EXP_NLICE1(1,10) = 1279528.12    
    EXP_NRAIN1(1,10) = 4552.00000  

    ! Test Case 11: RQR <= RQR_DRmin
    P1D(1,11) = 80000.0
    T1D(1,11) = 263.15
    Q1D(1,11) = 0.002
    C1D(1,11) = 1.0e-4
    FI1D(1,11) = 0.5
    FR1D(1,11) = 0.001
    FS1D(1,11) = 1.5
    CUREFL(1,11) = 0.0

    EXP_QW1(1,11) = 4.99499984E-05
    EXP_QI1(1,11) = 2.49608424E-11
    EXP_QR1(1,11) = 5.00000006E-08
    EXP_QS1(1,11) = 4.99999733E-05
    EXP_DBZ1(1,11) = 11.3324432    
    EXP_DBZR1(1,11) = -20.0000000    
    EXP_DBZI1(1,11) = 11.3319588    
    EXP_DBZC1(1,11) = -20.0000000    
    EXP_NLICE1(1,11) = 3363.78760    
    EXP_NRAIN1(1,11) = 134.687271  

    ! Test Case 12: RQR >= RQR_DRmax
    P1D(1,12) = 80000.0
    T1D(1,12) = 263.15
    Q1D(1,12) = 0.002
    C1D(1,12) = 1.0e-1
    FI1D(1,12) = 0.0
    FR1D(1,12) = 1.0
    FS1D(1,12) = 1.5
    CUREFL(1,12) = 0.0

    EXP_QW1(1,12) = 0.0   
    EXP_QI1(1,12) = 0.0   
    EXP_QR1(1,12) = 0.100000001    
    EXP_QS1(1,12) = 0.0   
    EXP_DBZ1(1,12) = 73.8459854    
    EXP_DBZR1(1,12) = 73.8459854    
    EXP_DBZI1(1,12) = -20.0   
    EXP_DBZC1(1,12) = -20.0    
    EXP_NLICE1(1,12) = 0.0    
    EXP_NRAIN1(1,12) = 33671.8008 

    ! Test Case 13: RimeF > 10
    P1D(1,13) = 80000.0
    T1D(1,13) = 263.15
    Q1D(1,13) = 0.002
    C1D(1,13) = 1.0e-4
    FI1D(1,13) = 1.0
    FR1D(1,13) = 0.0
    FS1D(1,13) = 11.0
    CUREFL(1,13) = 0.0

    EXP_QW1(1,13) = 0.0  
    EXP_QI1(1,13) = 2.49608424E-11
    EXP_QR1(1,13) = 0.0    
    EXP_QS1(1,13) = 9.99999756E-05
    EXP_DBZ1(1,13) = 18.0097847    
    EXP_DBZR1(1,13) = -20.0   
    EXP_DBZI1(1,13) = 18.0097847    
    EXP_DBZC1(1,13) = -20.0    
    EXP_NLICE1(1,13) = 2891.38599    
    EXP_NRAIN1(1,13) = 0.0  

    ! Test Case 14: QI1 <= EPSQ
    P1D(1,14) = 80000.0
    T1D(1,14) = 263.15
    Q1D(1,14) = 0.002
    C1D(1,14) = 1.0e-4
    FI1D(1,14) = 0.0
    FR1D(1,14) = 0.25
    FS1D(1,14) = 1.5
    CUREFL(1,14) = 0.0

    EXP_QW1(1,14) = 7.49999963E-05
    EXP_QI1(1,14) = 0.0    
    EXP_QR1(1,14) = 2.49999994E-05
    EXP_QS1(1,14) = 0.0    
    EXP_DBZ1(1,14) = 15.4733009    
    EXP_DBZR1(1,14) = 15.4733009    
    EXP_DBZI1(1,14) = -20.0000000    
    EXP_DBZC1(1,14) = -20.0000000    
    EXP_NLICE1(1,14) = 0.0    
    EXP_NRAIN1(1,14) = 1440.0 

    ! Test Case 15: C1D <= EPSQ
    P1D(1,15) = 80000.0
    T1D(1,15) = 263.15
    Q1D(1,15) = 0.002
    C1D(1,15) = 0.0
    FI1D(1,15) = 0.5
    FR1D(1,15) = 0.25
    FS1D(1,15) = 1.5
    CUREFL(1,15) = 0.0

    EXP_QW1(1,15) = 0.0    
    EXP_QI1(1,15) = 0.0    
    EXP_QR1(1,15) = 0.0    
    EXP_QS1(1,15) = 0.0    
    EXP_DBZ1(1,15) = -20.0    
    EXP_DBZR1(1,15) = -20.0    
    EXP_DBZI1(1,15) = -20.0    
    EXP_DBZC1(1,15) = -20.0    
    EXP_NLICE1(1,15) = 0.0    
    EXP_NRAIN1(1,15) = 0.0  

    call CALMICT_new(P1D,T1D,Q1D,C1D,FI1D,FR1D,FS1D,CUREFL,   &
                     QW1,QI1,QR1,QS1,DBZ1,DBZR1,DBZI1,DBZC1, &
                     NLICE1,NRAIN1)

    res = 0
    do i = 1, ntests
        if (abs(QW1(1,i) - EXP_QW1(1,i)) > tol) then
            print *, "Test case ", i, " failed for QW1: Expected ", EXP_QW1(1,i), &
                     " but got ", QW1(1,i)
            res = 1
        end if
        if (abs(QI1(1,i) - EXP_QI1(1,i)) > tol) then
            print *, "Test case ", i, " failed for QI1: Expected ", EXP_QI1(1,i), &
                     " but got ", QI1(1,i)
            res = 1
        end if
        if (abs(QR1(1,i) - EXP_QR1(1,i)) > tol) then
            print *, "Test case ", i, " failed for QR1: Expected ", EXP_QR1(1,i), &
                     " but got ", QR1(1,i)
            res = 1
        end if
        if (abs(QS1(1,i) - EXP_QS1(1,i)) > tol) then
            print *, "Test case ", i, " failed for QS1: Expected ", EXP_QS1(1,i), &
                     " but got ", QS1(1,i)
            res = 1
        end if
        if (abs(DBZ1(1,i) - EXP_DBZ1(1,i)) > tol) then
            print *, "Test case ", i, " failed for DBZ1: Expected ", EXP_DBZ1(1,i), &
                     " but got ", DBZ1(1,i)
            res = 1
        end if
        if (abs(DBZR1(1,i) - EXP_DBZR1(1,i)) > tol) then
            print *, "Test case ", i, " failed for DBZR1: Expected ", EXP_DBZR1(1,i), &
                     " but got ", DBZR1(1,i)
            res = 1
        end if
        if (abs(DBZI1(1,i) - EXP_DBZI1(1,i)) > tol) then
            print *, "Test case ", i, " failed for DBZI1: Expected ", EXP_DBZI1(1,i), &
                     " but got ", DBZI1(1,i)
            res = 1
        end if
        if (abs(DBZC1(1,i) - EXP_DBZC1(1,i)) > tol) then
            print *, "Test case ", i, " failed for DBZC1: Expected ", EXP_DBZC1(1,i), &
                     " but got ", DBZC1(1,i)
            res = 1
        end if
        ! Using relative tolerance for NLICE1
        if (abs(NLICE1(1,i) - EXP_NLICE1(1,i))/EXP_NLICE1(1,i) > tol) then
            print *, "Test case ", i, " failed for NLICE1: Expected ", EXP_NLICE1(1,i), &
                     " but got ", NLICE1(1,i)
            res = 1
        end if
        if (abs(NRAIN1(1,i) - EXP_NRAIN1(1,i)) > tol) then
            print *, "Test case ", i, " failed for NRAIN1: Expected ", EXP_NRAIN1(1,i), &
                     " but got ", NRAIN1(1,i)
            res = 1
        end if
    end do

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"
end program test_calmict_new