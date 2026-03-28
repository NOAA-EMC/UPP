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

    real, parameter :: tol = 1.0e-8

    interface
        subroutine SELECT_CHANNELS_L(CHANNELINFO, NCHANNELS, CHANNELS, L, IGOT)
            use crtm_channelinfo_define, only: crtm_channelinfo_type
            type(crtm_channelinfo_type), intent(inout) :: CHANNELINFO
            integer, intent(in) :: NCHANNELS
            integer, intent(in) :: CHANNELS(NCHANNELS), L(NCHANNELS)
            integer, intent(out) :: IGOT
        end subroutine SELECT_CHANNELS_L
    end interface


    
    print *, 'SUCCESS!'
end program test_select_channels