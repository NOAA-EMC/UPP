! This is a test program for UPP.
!
! This program tests the SELECT_CHANNELS_L() subroutine in SELECT_CHANNELS.f.
!
! The subroutine SELECT_CHANNELS() is intentionally skipped because it's
! no longer being used.
!
! Alyson Stahl, 3/2026
program test_select_channels
    use crtm_channelinfo_define, only: crtm_channelinfo_type
    implicit none

    integer, parameter :: NCHANNELS = 5
    integer :: CHANNELS(1:NCHANNELS), L(1:NCHANNELS), IGOT
    type(crtm_channelinfo_type) :: CHANNELINFO
    integer :: i, res

    interface
        subroutine SELECT_CHANNELS_L(CHANNELINFO, NCHANNELS, CHANNELS, L, IGOT)
            use crtm_channelinfo_define, only: crtm_channelinfo_type
            type(crtm_channelinfo_type), intent(inout) :: CHANNELINFO
            integer, intent(in) :: NCHANNELS, CHANNELS(NCHANNELS)
            integer :: IGOT, L(NCHANNELS)
        end subroutine SELECT_CHANNELS_L
    end interface

    CHANNELINFO%Is_Allocated = .true.
    CHANNELINFO%n_channels = NCHANNELS
    allocate(CHANNELINFO%Channel_Index(NCHANNELS))
    allocate(CHANNELINFO%Process_Channel(NCHANNELS))
    allocate(CHANNELINFO%Sensor_Channel(NCHANNELS))

    do i = 1, NCHANNELS
        CHANNELINFO%Channel_Index(i) = 10 * i
        CHANNELINFO%Process_Channel(i) = .true.
        CHANNELINFO%Sensor_Channel(i) = i
    end do

    CHANNELS = (/ 1, 2, 3, 4, 5 /)

    res = 0

    ! Test Case 1: All L = 0
    L = 0
    IGOT = 99

    call SELECT_CHANNELS_L(CHANNELINFO, NCHANNELS, CHANNELS, L, IGOT)

    if (IGOT .ne. 0) then
        print *, 'Test Case 1 Failed: Expected IGOT = 0, got ', IGOT
        res = 1
    end if

    do i = 1, NCHANNELS
        if (CHANNELINFO%Process_Channel(i)) then
            print *, 'Test 1 failed: Process_Channel(', i, ') should be .FALSE.'
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(CHANNELINFO%Channel_Index)
        deallocate(CHANNELINFO%Process_Channel)
        deallocate(CHANNELINFO%Sensor_Channel)
        stop 10
    end if

    do i = 1, NCHANNELS
        CHANNELINFO%Process_Channel(i) = .true.
    end do

    ! Test case 2: Mixed L values
    L    = (/ 1, 0, 1, 0, 1 /)
    IGOT = 5

    call SELECT_CHANNELS_L(CHANNELINFO, NCHANNELS, CHANNELS, L, IGOT)
    
    if (IGOT .ne. 5) then
        print *, 'Test Case 2 Failed: Expected IGOT = 5, got ', IGOT
        res = 1
    end if

    do i = 1, NCHANNELS
        if (L(i) .eq. 0 .and. CHANNELINFO%Process_Channel(i)) then
            print *, 'Test 2 failed: Process_Channel(', i, ') should be .FALSE.'
            res = 1
        else if (L(i) .eq. 1 .and. .not. CHANNELINFO%Process_Channel(i)) then
            print *, 'Test 2 failed: Process_Channel(', i, ') should be .TRUE.'
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(CHANNELINFO%Channel_Index)
        deallocate(CHANNELINFO%Process_Channel)
        deallocate(CHANNELINFO%Sensor_Channel)
        stop 20
    end if

    do i = 1, NCHANNELS
        CHANNELINFO%Process_Channel(i) = .true.
    end do
    
    ! Test case 3: L contains a value other than 0 or 1
    L = (/ 2, 1, 0, 2, 1 /)
    IGOT = 7

    call SELECT_CHANNELS_L(CHANNELINFO, NCHANNELS, CHANNELS, L, IGOT)

    if (IGOT .ne. 7) then
        print *, 'Test Case 3 Failed: Expected IGOT = 7, got ', IGOT
        res = 1
    end if

    do i = 1, NCHANNELS
        ! Only i = 3 should be turned off
        if (i .eq. 3 .and. CHANNELINFO%Process_Channel(i)) then
            print *, 'Test 3 failed: Process_Channel(', i, ') should be .FALSE.'
            res = 1
        else if (i .ne. 3 .and. .not. CHANNELINFO%Process_Channel(i)) then
            print *, 'Test 3 failed: Process_Channel(', i, ') should be .TRUE.'
            res = 1
        end if
    end do

    if (res .ne. 0) then
        deallocate(CHANNELINFO%Channel_Index)
        deallocate(CHANNELINFO%Process_Channel)
        deallocate(CHANNELINFO%Sensor_Channel)
        stop 30
    end if

    deallocate(CHANNELINFO%Channel_Index)
    deallocate(CHANNELINFO%Process_Channel)
    deallocate(CHANNELINFO%Sensor_Channel)
    print *, 'SUCCESS!'
end program test_select_channels