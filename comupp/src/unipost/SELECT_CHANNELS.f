      subroutine SELECT_CHANNELS(channelinfo,nchannels,channels)
        use crtm_channelinfo_define, only: crtm_channelinfo_type
        implicit none
        type(crtm_channelinfo_type),intent(inout) :: channelinfo
        integer, intent(in) :: nchannels,channels(nchannels)
        integer :: i,j
        integer :: temp(nchannels)

        if(nchannels>channelinfo%n_channels) then
           write(6,*) 'ERROR*** tried to use more channels than sensor has'
           write(6,*) '  ',nchannels,' > ',channelinfo%n_channels
           stop 18
        endif
        check: do i=1,nchannels
           if(channels(i)<1 .or. channels(i)>channelinfo%n_channels) then
              write(6,*) 'ERROR*** invalid channel id: ',channels(i)
              write(6,*) '  in SELECT_CHANNELS at index ',i
              stop 19
           endif
           temp(i)=channelinfo%Channel_Index(channels(i))
        enddo check

        channelinfo%n_channels=nchannels
        channelinfo%Channel_Index(1:nchannels)=temp
      end subroutine SELECT_CHANNELS
