!       function written early Dec. 1999 by M. Pyle to support  workstation
!       Eta for users with etime but not timef functionality (like  certain
!mp     HPs)  Designed to duplicate timef (elapsed time in milliseconds)
!
        function timef()
        use mpi
        implicit none
        real *8  et(2),rtc
        data et/0.0,0.0/
        real*8 timef, etime
        if(et(1) .eq. 0) et(1)=mpi_wtime()
        et(2)=mpi_wtime()
        timef=(et(2)-et(1))
!        timef=(et(2)-et(1))*1.e3
!        timef=mpi_wtime() *1.e3 -ti
        end

        function rtcfake()
        real*8 rtc, etime
        rtcfake=mpi_wtime() *1.e3
!        rtcfake=rtc*1.e3
        end
