program satgrib2
  implicit none
  integer :: ipack, iresult, iunit, ounit
  character(len=:), pointer :: infile, outfile
  character(len=5), parameter :: text_version = '0.0.1'

  integer,parameter :: mingrib=500
  integer, parameter :: maxpts=40000000,msk1=32000

  iresult = scan_arg_list(ipack,infile,outfile)
  if(iresult==0) goto 999 ! satgrib2 -h with no other args, return 0
  if(iresult<0) stop 1 ! error, exit non-zero

1 format('satgrib2 -p',I0,' "',A,'" "',A,'"')
  print 1,ipack,infile,outfile

  call open_files(infile,outfile,iunit,ounit)
  call cnv12(infile,outfile,iunit,ounit)
  call close_files(infile,outfile,iunit,ounit)

999 continue ! avoid printing "0" when running satgrib2 -h

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer function channelmap(nchannels,channels,level)
    integer, intent(inout) :: level
    integer, intent(in) :: nchannels, channels(nchannels)
    if(level<1 .or. level>nchannels) then
       channelmap=-1
       write(0,35) level,nchannels
       return
    endif
    level=channels(level)
    channelmap=0
35  format("Bad level ",I0," not in [1,",I0,"].  Ignoring record.")
  end function channelmap

  integer function hlev2table432(kpds,listsec1,ipdsnum,ipdstmpl)
    integer, intent(in) :: kpds(200)
    integer, intent(inout) :: listsec1(13), ipdsnum, ipdstmpl(15)


!  Read_SpcCoeff_Binary(INFORMATION) : FILE: lib/crtm2/src/fix/SpcCoeff/Big_Endian/ahi_himawari8.SpcCoeff.bin; 
!  SpcCoeff RELEASE.VERSION:  7.01  N_CHANNELS=10
! lib/crtm2/src/fix/SpcCoeff/Big_Endian/ahi_himawari8.SpcCoeff.bin: 
!    Has 10 channels, 0 FOVs, release=7 version=1
!    WMO satellite id=173 WMO sensor id=297
!    Sensor "ahi_himawari8"
!    RCS id "$Id: SpcCoeff_Define.f90 22707 2012-11-21 21:09:10Z paul.vandelst@noaa.gov $"
!    Band   1: wavenumber=257593.86633 m-1, wavelength=3.882 um, freq=77224.69 Ghz, Unpolarized
!    Band   2: wavenumber=160990.92907 m-1, wavelength=6.212 um, freq=48263.86 Ghz, Unpolarized
!    Band   3: wavenumber=143921.76108 m-1, wavelength=6.948 um, freq=43146.66 Ghz, Unpolarized
!    Band   4: wavenumber=136199.66035 m-1, wavelength=7.342 um, freq=40831.63 Ghz, Unpolarized
!    Band   5: wavenumber=116446.90808 m-1, wavelength=8.588 um, freq=34909.90 Ghz, Unpolarized
!    Band   6: wavenumber=103906.23067 m-1, wavelength=9.624 um, freq=31150.30 Ghz, Unpolarized
!    Band   7: wavenumber=96145.03574 m-1, wavelength=10.401 um, freq=28823.55 Ghz, Unpolarized
!    Band   8: wavenumber=89155.43992 m-1, wavelength=11.216 um, freq=26728.13 Ghz, Unpolarized
!    Band   9: wavenumber=80976.33020 m-1, wavelength=12.349 um, freq=24276.09 Ghz, Unpolarized
!    Band  10: wavenumber=75413.09962 m-1, wavelength=13.260 um, freq=22608.28 Ghz, Unpolarized

    ! Central wavenumber in m^-1 * 100
    integer, parameter :: h8ahi(10) = (/ 25759386, 16099092, 14392176, 13619966, &
         11644690, 10390623, 9614503, 8915543, 8097633, 7541309 /)
    integer, parameter :: seviri(8) = (/  &
         25658851,15945984,13598205,11478975,10346736,9287793,8379711,7496626 /)
    integer, parameter :: goes12(4) = (/ 25648222,15423775,9334097,7511192 /)
    integer, parameter :: goes13(4) = (/ 25639572,15282016,9372177,7526775 /)
    integer, parameter :: goes15(4) = (/ 25639344,15261432,9359576,7531833 /)
    integer, parameter :: insat3d(4) = (/ 25484201,14544097,9250980,8373859 /)
    integer, parameter :: mtsat1r(4) = (/ 9269174,8331543,14823088,26536647 /)
    integer, parameter :: mtsat2(4) = (/ 9264610,8356672,14766938,26841158 /)
    integer, parameter :: ssmif13(7) = (/ 6454,6454,7416,12341,12341,28519,28519 /)
    integer, parameter :: ssmif14(7) = (/ 6454,6454,7416,12341,12341,28519,28519 /)
    integer, parameter :: ssmif15(7) = (/ 6454,6454,7416,12341,12341,28519,28519 /)
    integer, parameter :: ssmisf16(24) = (/ &
         16778,17612,17877,18145,18512,19109,19813,50034,61145,61145,61145,6454,&
         6454,7416,12341,12341,30572,30572,21109,20278,20278,20278,20278,20278 /)
    integer, parameter :: ssmisf17(24) = (/ &
         16778,17612,17877,18145,18512,19109,19813,50034,61145,61145,61145,6454,&
         6454,7416,12341,12341,30572,30572,21109,20278,20278,20278,20278,20278 /)
    integer, parameter :: ssmisf18(24) = (/ &
         16778,17612,17877,18145,18512,19109,19813,50034,61145,61145,61145,6454,&
         6454,7416,12341,12341,30572,30572,21109,20278,20278,20278,20278,20278 /)
    integer, parameter :: ssmisf19(24) = (/ &
         16778,17612,17877,18145,18512,19109,19813,50034,61145,61145,61145,6454,&
         6454,7416,12341,12341,30572,30572,21109,20278,20278,20278,20278,20278 /)
    integer, parameter :: ssmisf20(24) = (/ &
         16778,17612,17877,18145,18512,19109,19813,50034,61145,61145,61145,6454,&
         6454,7416,12341,12341,30572,30572,21109,20278,20278,20278,20278,20278 /)
    
    integer :: level, satellite, channel

    integer, parameter :: Unpolarized = 1*8192, H_Lin_Pol = 2*8192
    integer, parameter :: V_Lin_Pol=3*8192, R_Circ_Pol=4*8192

    ! Polarization information for F13-F20 SSMI and SSMIS
    integer, parameter :: polf13(7) = (/ &
         V_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,V_Lin_Pol,H_Lin_Pol /)
    integer, parameter :: polf14(7) = (/ &
         V_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,V_Lin_Pol,H_Lin_Pol /)
    integer, parameter :: polf15(7) = (/ &
         V_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,V_Lin_Pol,H_Lin_Pol /)
    integer, parameter :: polf16(24) = (/ &
         V_Lin_Pol,V_Lin_Pol,V_Lin_Pol,V_Lin_Pol,V_Lin_Pol,R_Circ_Pol,R_Circ_Pol,&
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,&
         H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,R_Circ_Pol,&
         R_Circ_Pol,R_Circ_Pol,R_Circ_Pol /)
    integer, parameter :: polf17(24) = (/ &
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,&
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,&
         H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,R_Circ_Pol,&
         R_Circ_Pol,R_Circ_Pol,R_Circ_Pol /)
    integer, parameter :: polf18(24) = (/ &
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,&
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,&
         H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,R_Circ_Pol,&
         R_Circ_Pol,R_Circ_Pol,R_Circ_Pol /)
    integer, parameter :: polf19(24) = (/ &
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,&
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,&
         H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,R_Circ_Pol,&
         R_Circ_Pol,R_Circ_Pol,R_Circ_Pol /)
    integer, parameter :: polf20(24) = (/ &
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,&
         H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,&
         H_Lin_Pol,V_Lin_Pol,V_Lin_Pol,H_Lin_Pol,R_Circ_Pol,R_Circ_Pol,R_Circ_Pol,&
         R_Circ_Pol,R_Circ_Pol,R_Circ_Pol /)

    hlev2table432=5

    !7, 89, 193, 192, 118, 109,
    if(kpds(2)==129 .and. kpds(6)>=213 .and. kpds(6)<=216) then
       ! GOES 12, nadir view.  Trick code below into working.
       satellite=12
       channel=kpds(6)-212
       level=satellite*100+channel
   elseif(kpds(1)==7 .and. kpds(2)==89 .and. kpds(3)==193 .and. kpds(4)==192 .and. kpds(5)==118 .and. kpds(6)==109) then
      level=kpds(7)
      satellite=level/100
      channel=level-satellite*100
   elseif(kpds(1)==7 .and. kpds(2)==89 .and. kpds(3)==255 .and. kpds(4)==128 .and. kpds(5)==118 .and. kpds(6)==109) then
      level=kpds(7)
      satellite=level/100
      channel=level-satellite*100
    elseif(kpds(1)/=7 .or. kpds(2)/=89 .or. (kpds(3)/=193 .and. kpds(3)/=255) .or. kpds(4)/=192 .or. &
         kpds(5)/=118 .or. kpds(6)/=109) then
       ! Unknown satellite or not a satellite product.
38     format('Ignoring non-satellite record with kpds(1:6) = (/ ',5(I0,", "),I0,' /)')
       write(0,38) kpds(1:6)
       listsec1=0
       ipdsnum=0
       ipdstmpl=0
       hlev2table432=1
       return
    else
       level=kpds(7)
       satellite=level/100
       channel=level-satellite*100
    end if

    listsec1=0
    listsec1(1)=kpds(1)   ! Id of orginating centre (Common Code Table C-1)
    listsec1(2)=kpds(23)  ! Id of orginating sub-centre (local table)/Table C of ON388
    listsec1(3)=8         ! GRIB Master Tables Version Number (Code Table 1.0)
    listsec1(4)=1         ! GRIB Local Tables Version Number (Code Table 1.1)
    listsec1(5)=1         ! Significance of Reference Time (Code Table 1.2)
    listsec1(6)=(kpds(21)-1)*100+kpds(8)   ! Reference Time - Year (4 digits)
    listsec1(7)=kpds(9)   ! Reference Time - Month
    listsec1(8)=kpds(10)  ! Reference Time - Day
    listsec1(9)=kpds(11)  ! Reference Time - Hour
    listsec1(10)=kpds(12) ! Reference Time - Minute
    listsec1(11)=0       ! Reference Time - Second
    listsec1(12)=0       ! Production status of data (Code Table 1.3)
    listsec1(13)=1       ! Type of processed data (Code Table 1.4)

    ipdsnum=32               ! Satellite
    ipdstmpl=0
    ipdstmpl(1)=4            ! Short-wave radiation
    ipdstmpl(2)=4            ! Short-wave brightness temperature
    ipdstmpl(3)=2            ! Forecast
    ipdstmpl(4:5)=0          ! Generating process identifier
    ipdstmpl(6)=0            ! hours after reference time
    ipdstmpl(7)=0            ! minutes after reference time
    ipdstmpl(8)=kpds(13)     ! forecast lead time units: hours
    ipdstmpl(9)=kpds(14)     ! forecast lead time
    ipdstmpl(10)=1           ! Number of spectral bands
    
    ! Product description: first and only spectral band
    ipdstmpl(11)=511         ! satellite constellation (BUFR 0 02 020)
    ipdstmpl(12)=854         ! satellite (BUFR 0 01 007)
    ipdstmpl(13)=2047        ! instrument (BUFR 0 02 019) and polarization
    ipdstmpl(14)=2           ! Scale factor for wavenumber
    ipdstmpl(15)=0           ! Center wavenumber is 0 m-1

    select case(satellite)
    case( 10) ! SEVIRI M10 
       hlev2table432=channelmap(7,(/ 2,3,4,5,6,7,8 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 333,  57, 207+unpolarized,   2, seviri(channel)   /)
    case( 13) ! GOES 13 imager
       hlev2table432=channelmap(4,(/ 1,2,3,4 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 241, 257, 615+unpolarized,   2, goes13(channel)   /)
    case( 12) ! GOES 13 imager
       hlev2table432=channelmap(4,(/ 1,2,3,4 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 241, 257, 615+unpolarized,   2, goes12(channel)   /)
    case( 15) ! GOES 15 imager
       hlev2table432=channelmap(4,(/ 1,2,3,4 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 241, 259, 615+unpolarized,   2, goes15(channel)   /)
    case( 30) ! INSAT-3D 
       hlev2table432=channelmap(4,(/ 1,2,3,4 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 301, 471, 270+unpolarized,   2, insat3d(channel)  /)
    case(100) ! MTSAT-1r JAMI imager
       hlev2table432=channelmap(4,(/ 1,2,3,4 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 272, 172, 294+unpolarized,   2, mtsat1r(channel)  /)
    case(200) ! MTSAT-2 imager
       hlev2table432=channelmap(4,(/ 1,2,3,4 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 272, 172, 295+unpolarized,   2, mtsat2(channel)   /)
    case(280) ! Himawari 8 AHI
       hlev2table432=channelmap(10,(/ 1,2,3,4,5,6,7,8,9,10 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/ 273, 173, 297+unpolarized,   2, h8ahi(channel) /)
    case(113) ! F13 SSM/I
       hlev2table432=channelmap(6,(/ 1,2,4,5,6,7 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 246, 905+polf13(channel), 2, ssmif13(channel)  /)
    case(114) ! F14 SSM/I
       hlev2table432=channelmap(6,(/ 1,2,4,5,6,7 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 247, 905+polf14(channel), 2, ssmif14(channel)  /)
    case(115) ! F15 SSM/I
       hlev2table432=channelmap(6,(/ 1,2,4,5,6,7 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 248, 905+polf15(channel), 2, ssmif15(channel)  /)
    case(116) ! F16 SSMIS
       hlev2table432=channelmap(7,(/ 9,12,13,15,16,17,18 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 249, 908+polf16(channel), 2, ssmisf16(channel) /)
    case(117) ! F17 SSMIS
       hlev2table432=channelmap(7,(/ 9,12,13,15,16,17,18 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 285, 908+polf17(channel), 2, ssmisf17(channel) /)
    case(118) ! F18 SSMIS
       hlev2table432=channelmap(7,(/ 9,12,13,15,16,17,18 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 286, 908+polf18(channel), 2, ssmisf18(channel) /)
    case(119) ! F19 SSMIS
       hlev2table432=channelmap(7,(/ 9,12,13,15,16,17,18 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 287, 908+polf19(channel), 2, ssmisf19(channel) /)
    case(120) ! F20 SSMIS not in 0 01 007 yet
       hlev2table432=channelmap(7,(/ 9,12,13,15,16,17,18 /),channel)
       if(hlev2table432/=0) return
       ipdstmpl(11:15)=(/  31, 853, 908+polf20(channel), 2, ssmisf20(channel) /)
    end select

    if(ipdstmpl(15)==0) then
       ! Could not find sensor.
       hlev2table432=1
       return
    endif

21  format('Satellite ',I0,' ',I0,' instrument ',I0,' polarization ',I0,&
           ' wavenumber ',I0,'/100')
    print 21,ipdstmpl(11),ipdstmpl(12),iand(ipdstmpl(13),2047),&
             ipdstmpl(13)/8192,ipdstmpl(15)

    if(ipdstmpl(15)<3e7) then
       ! Switch to longwave
       ipdstmpl(1:2)=(/ 5, 7 /)
    endif

    hlev2table432=0
  end function hlev2table432

  subroutine cnv12(infile,outfile,iunit,ounit)
    implicit none

    character(len=:), pointer, intent(out) :: infile, outfile
    integer, intent(in) :: iunit,ounit
    
    integer :: iseek, ifli1, lskip, lgrib, currlen, lcgrib, ifield

    integer :: listsec0(2),listsec1(13)
    logical*1,allocatable,dimension(:) :: bmp,bmpv
    real,allocatable,dimension(:) :: fld
    real,allocatable,dimension(:) :: fldv
    real,allocatable,dimension(:) :: coordlist
    character(len=1),allocatable,dimension(:) :: cgrib,cgribin
    integer :: idrstmpl(200),idrstmplv(200)
    integer ::  KPDS(200),KGDS(200),KPTR(200), i
    integer :: ideflist(MAXPTS),idefnum,idrsnum,numcoord
    integer :: igds(5)=(/0,0,0,0,0/),igdstmpl(200), ibmap
    integer :: icnd, lengrib , iret, numpts, ipdstmpl(200), ipdsnum, ierr
    logical :: usemiss
    real :: rmiss

    usemiss=.false.

    icnd=0
    ifli1=0
    iseek=0
    currlen=0
    listsec0=(/0,2/)
    ifield=0
    lgrib=0
    lskip=0
    allocate(fld(maxpts))
    allocate(coordlist(maxpts))
    coordlist=0
    allocate(bmp(maxpts))

    gribloop: do
       lgrib=-1
       call skgb(iunit,iseek,msk1,lskip,lgrib)
       if(lgrib==0) then
          write(0,*) 'Exit loop: lgrib=',lgrib
          exit gribloop
       endif
       need_more_space: if(lgrib>currlen) then
          if(allocated(cgribin)) deallocate(cgribin)
          if(allocated(cgrib)) deallocate(cgrib)

          allocate(cgribin(lgrib))
          currlen=lgrib
          lcgrib=max(mingrib,lgrib*2)
          allocate(cgrib(lcgrib))
       end if need_more_space

       call baread(iunit,lskip,lgrib,lengrib,cgribin)
       if(lgrib==lengrib) then
          iret=0
          call w3fi63(cgribin,KPDS,KGDS,BMP,FLD,KPTR,IRET)
          numpts=KPTR(10)
          if (iret.ne.0) then
             call warning('error unpacking GRIB field',iret,infile)
             iseek=lskip+lgrib
             cycle gribloop
          endif
        else
           call warning('I/O error on input grib file',lgrib-lengrib,infile)
           cycle gribloop
       endif
       iseek=lskip+lgrib
       ifield=ifield+1

40     format("Field ",I4," kpds(1..23) = ",23(I0,", "))

       ! Get listsec1 and GRIB2 PDS for this satellite
       iret=hlev2table432(kpds,listsec1,ipdsnum,ipdstmpl)

       if(iret/=0) cycle gribloop ! skip non-satellite fields and bad records

       print 40,ifield,kpds(1:23)

       call gribcreate(cgrib,lcgrib,listsec0,listsec1,ierr)
       if (ierr.ne.0) then
          write(6,*) ' ERROR creating new GRIB2 field = ',ierr
          cycle gribloop
       endif

       ! Convert GDS from grib1 to grib2
        call gds2gdt(kgds,igds,igdstmpl,idefnum,ideflist,ierr)
        if (ierr.ne.0) then
           write(0,*) 'ERROR converting GDS from GRIB1 to GRIB2',ierr
           cycle gribloop
        endif
        if (listsec1(1) .eq. 7 ) igdstmpl(1)=6    ! FOR NWS/NCEP
        if ((listsec1(1) .eq. 7 .and. igds(5).eq.20 &   ! For Snow Cover Analysis 
             .and. kpds(2).eq.25 ) .and. &              ! Polar Stereographic Grid
             (kpds(5).eq.91 .or. kpds(5).eq.238)) then 
           igdstmpl(1)=2
        end if 
        call addgrid(cgrib,lcgrib,igds,igdstmpl,200,ideflist,idefnum,ierr)
        if (ierr.ne.0) then
          write(6,*) ' ERROR adding GRIB2 grid = ',ierr
          cycle
        endif

        ! Handle bitmap:
        idrstmpl=0
        if (btest(kpds(4),6)) then
          ibmap=0
          !fld=pack(fld,mask=bmp(1:numpts))
          !itemp=count(bmp(1:numpts))
          !numpts=itemp
          !
          !   convert bitmap to "missing" values, if requested.
          !
          if ( (usemiss) .AND. (ipack.eq.2 .OR. ipack.eq.31 .OR. &
                                ipack.eq.32) ) then
             ibmap=255
             rmiss=minval(fld(1:numpts))
             if ( rmiss .lt. -9999.0 ) then
                rmiss=rmiss*10.0
             else
                rmiss=-9999.0
             endif
             do i=1,numpts
                if ( .NOT. bmp(i) ) then
                   fld(i)=rmiss
                   bmp(i)=.true.
                endif
             enddo
             idrstmpl(7)=1                   ! Missing value management 
             call mkieee(rmiss,idrstmpl(8),1)
          endif
        else
          ibmap=255
          idrstmpl(7)=0                   ! No missing values
        endif
        
        ! Set packing info:
        if ( ipack.eq.0 ) then
           idrsnum=0
        elseif ( ipack.eq.2 ) then
           idrsnum=2
           idrstmpl(6)=1                   ! general group split
        elseif ( ipack.eq.31.OR.ipack.eq.32 ) then
           idrsnum=ipack/10
           idrstmpl(6)=1                   ! general group split
           idrstmpl(17)=mod(ipack,10)      ! order of s.d.
        elseif ( ipack.eq.40 .OR. ipack.eq.41 .OR. &
                 ipack.eq.40000 .OR. ipack.eq.40010 ) then
           idrsnum=ipack
           idrstmpl(6)=0
           idrstmpl(7)=255
        else
           idrsnum=3
           idrstmpl(17)=1                  ! order of s.d.
           idrstmpl(6)=1                   ! general group split
           if (kpds(5).eq.61) idrsnum=2
        endif
        idrstmpl(2)=KPTR(19)       ! binary scale
        idrstmpl(3)=kpds(22)       ! decimal scale

        numcoord=0
        call addfield(cgrib,lcgrib,ipdsnum,ipdstmpl,15,&
                      coordlist,numcoord,idrsnum,idrstmpl,200,&
                      fld,numpts,ibmap,bmp,ierr)
        call gribend(cgrib,lcgrib,lengrib,ierr)
        if (ierr.ne.0) then
          write(6,*) ' ERROR ending new GRIB2 message = ',ierr
          cycle
        endif
!        print *,' writing ',lengrib,' bytes...'
        call wryte(ounit,lengrib,cgrib)

    end do gribloop

    if (allocated(cgribin)) deallocate(cgribin)
    if (allocated(cgrib)) deallocate(cgrib)
    if (allocated(fld)) deallocate(fld)
    if (allocated(fldv)) deallocate(fldv)
    if (allocated(coordlist)) deallocate(coordlist)
    if (allocated(bmp)) deallocate(bmp)
    
    return
  end subroutine cnv12

  subroutine abort(why,ierr,filename)
    character(len=*), intent(in) :: why
    character(len=*), intent(in),optional :: filename
    integer, intent(in) :: ierr
    call warning(why,ierr,filename)
    stop 1
  end subroutine abort

  subroutine warning(why,ierr,filename)
    character(len=*), intent(in) :: why
    character(len=*), intent(in),optional :: filename
    integer, intent(in) :: ierr

10  format(A,': ',A,' (ierr=',I0,')')
20  format('FATAL ERROR: ',A,' (ierr=',I0,')')

    if(present(filename)) then
       write(0,10) filename,why,ierr
    else
       write(0,20) why,ierr
    end if
  end subroutine warning

  subroutine open_files(infile,outfile,iunit,ounit)
    character(len=:), pointer, intent(out) :: infile, outfile
    integer, intent(inout) :: iunit,ounit
    integer :: ios

    iunit=10
    call BAOPENR(iunit,infile,ios)
    if(ios/=0) call abort('cannot open for read',ios,filename=infile)

    ounit=50
    call BAOPENW(ounit,outfile,ios)
    if(ios/=0) call abort('cannot open for write',ios,filename=outfile)
  end subroutine open_files

  subroutine close_files(infile,outfile,iunit,ounit)
    character(len=:), pointer, intent(out) :: infile, outfile
    integer, intent(inout) :: iunit,ounit
    integer :: ios

    ios=0
    CALL BACLOSE(iunit,IOS)
    ios=0
    CALL BACLOSE(ounit,IOS)
    if(ios/=0) call warning('error closing input file',ios,filename=infile)
    if(ios/=0) call warning('error closing output file',ios,filename=outfile)
  end subroutine close_files

  integer function scan_arg_list(ipack,infile,outfile)
    implicit none
    character(len=:), pointer, intent(out) :: infile, outfile
    integer, intent(out) :: ipack

    character(len=:), pointer :: arg, message
    integer :: arglen, nargs, iarg, filearg
    nargs=command_argument_count()
    if(nargs<1) then
       call usage("Please specify the input and output file.")
       scan_arg_list=-1
       return
    endif

    nullify(arg,message)
    ipack=32
    
20  format('satgrib2 version ',A)
    scanargs: do iarg=1,nargs
       filearg=iarg
       if(associated(arg)) deallocate(arg)
       call get_command_argument(iarg,length=arglen)
       allocate(character(len=arglen) :: arg)
       call get_command_argument(iarg,arg)

       if(arg(1:1) /= '-') exit scanargs

       select case(arg)
       case('--')
          filearg=iarg+1
          exit scanargs
       case('--version')
          print 20, text_version
          scan_arg_list=0
          return
       case('-h','--help')
          call usage('-h')
          if(nargs==1) then
             scan_arg_list=0
          else
             scan_arg_list=-1
          end if
          return
       case('-g12')   ; cycle scanargs ! ignore -g12
       case('-p0')    ; ipack=0
       case('-p2')    ; ipack=2
       case('-p31')   ; ipack=31
       case('-p32')   ; ipack=32
       case('-p40')   ; ipack=40
       case('-p41')   ; ipack=41
       case default
10        format('Invalid argument: "',A,'"')
          allocate(character(len=len(arg)+20) :: message)
          write(message,10) arg
          call usage(message)
          scan_arg_list=-1
          return
       end select
    end do scanargs

    if(associated(arg)) deallocate(arg)

    if(filearg/=nargs-1) then
       call usage('Exactly two filenames are required: ingrib1 and outgrib2.')
       scan_arg_list=-1
       return
    endif

    call get_command_argument(filearg,length=arglen)
    allocate(character(len=arglen) :: infile)
    call get_command_argument(filearg,infile)
    
    call get_command_argument(filearg+1,length=arglen)
    allocate(character(len=arglen) :: outfile)
    call get_command_argument(filearg+1,outfile)

    scan_arg_list=1
  end function scan_arg_list

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine usage(why)
    implicit none
    logical :: full
    character(len=*), intent(in), optional :: why
    integer :: u

10  format(A)
20  format('INVALID ARGUMENTS: ',A)

    if(.not.present(why)) then
       ! run without arguments, so print short message to stdout
       full=.false.
       u=6
    else if(why == "-h") then
       ! write to stdout with full message if -h given
       u=6
       full=.true.
    else ! error message present, so print to stderr
       u=0 ! errors go to stdout
       full=.false.
    end if

    write(u,10) 'Program: satgrib2'
    write(u,10) 'Synopsis: Converts HWRF post satellite GRIB1 output to'
    write(u,10) '          GRIB2 synthetic satellite template 4.32.'

    if_full: if(full) then
       write(u,10) 'http://www.nco.ncep.noaa.gov/pmb/docs/grib2/grib2_temp4-32.shtml'
       write(u,10) ' '
       write(u,10) 'Usage: satgrib2 [options] ingrib1 outgrib2'
       write(u,10) 'Usage: satgrib2 [ -h | --help ]   (print this message)'
       write(u,10) 'Usage: satgrib2 --version         (print version message)'
       write(u,10) ' '
       write(u,10) 'Arguments:'
       write(u,10) '   ingrib1 = input GRIB1 hwrfsat file from HWRF satellite post'
       write(u,10) '   outgrib2 = output GRIB2 file that will contain template 4.32'
       write(u,10) '              ("synthetic satellite") data'
       write(u,10) ' '
       write(u,10) 'Options:'
       write(u,10) '   -g12 = ignored'
       write(u,10) '   -p0  = simple packing'
       write(u,10) '   -p2  = complex packing'
       write(u,10) '   -p31 = complex packing with 1st order diffs'
       write(u,10) '   -p32 = complex packing with 1st order diffs'
       write(u,10) '   -p40 = JPEG2000 encoding'
       write(u,10) '   -p41 = PNG encoding'
       write(u,10) '   --   = terminate option parsing (use to handle filenames'
       write(u,10) '          that begin with a hyphen)'
    else
       write(u,10) 'Usage: satgrib2 [-g12] [{-p0|-p2|-p31|-p32|-p40|-p41}]'
       write(u,10) '                ingrib1 outgrib2'
       write(u,10) 'Usage: satgrib2 [ -h | --help ]   print detailed usage message'
       write(u,10) 'Usage: satgrib2 --version         print version message'
    end if if_full

    if(present(why)) then
       if(why /= '-h') then
          write(u,10) ' '
          write(u,20) why
       endif
    endif
  end subroutine usage
end program satgrib2
