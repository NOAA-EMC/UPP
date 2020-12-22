!       function written early Dec. 1999 by M. Pyle to support  workstation
!       Eta for users with etime but not timef functionality (like  certain
!mp     HPs)  Designed to duplicate timef (elapsed time in milliseconds)
!
        function timef()
        implicit none
        real et(2)
        real*8 timef, etime
        timef=etime(et)
        timef=timef*1.e3
        end

        function rtc()
        implicit none
        real et(2)
        real*8 rtc, etime
        rtc=etime(et)
        rtc=rtc*1.e3
        end
