module fancy_module
    USE SpcCoeff_Binary_IO
    USE Type_Kinds         , ONLY: Long
    USE File_Utility       , ONLY: File_Open, File_Exists
    USE Message_Handler    , ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, Display_Message
    USE Binary_File_Utility, ONLY: Open_Binary_File
    USE SpcCoeff_Define    , ONLY: SpcCoeff_type, &
                                   Associated_SpcCoeff, &
                                   Allocate_SpcCoeff, &
                                   Destroy_SpcCoeff, &
                                   CheckRelease_SpcCoeff, &
                                   Info_SpcCoeff
    USE AntCorr_Binary_IO  , ONLY: Read_AntCorr_Binary, Write_AntCorr_Binary
    USE SpcCoeff_Define, ONLY: Destroy_SpcCoeff
  CHARACTER(*), PRIVATE, PARAMETER :: MODULE_RCS_ID = &
    '$Id: SpcCoeff_Binary_IO.f90 22707 2012-11-21 21:09:10Z paul.vandelst@noaa.gov $'
  ! Keyword set value
  INTEGER, PRIVATE, PARAMETER :: SET = 1
  ! Message character length
  INTEGER, PARAMETER :: ML = 512

  CONTAINS
    FUNCTION Fancy_Inquire_SpcCoeff_Binary( Filename        , &  ! Input
                                    n_Channels      , &  ! Optional Output
                                    n_FOVs          , &  ! Optional Output
                                    Release         , &  ! Optional Output
                                    Version         , &  ! Optional Output
                                    Sensor_Id       , &  ! Optional Output
                                    WMO_Satellite_Id, &  ! Optional Output
                                    WMO_Sensor_Id   , &  ! Optional Output
                                    RCS_Id          , &  ! Revision control
                                    Message_Log     , &  ! Error messaging
                                    Wavenumber      , &  ! Optional Output
                                    Polarization    ) &  ! Optional Output
                                  RESULT ( Error_Status )

      ! This function is similar to Inquire_SpcCoeff_Binary but also
      ! provides wavenumber and polarization information.

    ! Arguments
    USE Type_Kinds           , ONLY: Long, Double
    implicit none
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Channels
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_FOVs
    INTEGER     , OPTIONAL, INTENT(OUT) :: Release
    INTEGER     , OPTIONAL, INTENT(OUT) :: Version
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Sensor_Id       
    INTEGER     , OPTIONAL, INTENT(OUT) :: WMO_Satellite_Id
    INTEGER     , OPTIONAL, INTENT(OUT) :: WMO_Sensor_Id   
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: RCS_Id
    REAL(Double), OPTIONAL, POINTER, INTENT(OUT) :: Wavenumber(:)
    INTEGER(Long), OPTIONAL, POINTER, INTENT(OUT) :: Polarization(:)
    CHARACTER(*), OPTIONAL, INTENT(IN)  :: Message_Log
    ! Function result
    INTEGER :: Error_Status, Ignore
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Inquire_SpcCoeff_Binary'
    ! Local variables
    CHARACTER(ML) :: Message
    CHARACTER(ML) :: Close_Message
    INTEGER :: IO_Status
    INTEGER :: FileID
    TYPE(SpcCoeff_type) :: SpcCoeff

    Error_Status=Read_SpcCoeff_Binary(Filename=Filename,SpcCoeff=SpcCoeff,RCS_Id=RCS_Id,Message_Log=Message_Log)

    ! Assign the return arguments
    IF ( PRESENT(n_Channels      ) ) n_Channels      = SpcCoeff%n_Channels      
    IF ( PRESENT(n_FOVs          ) ) n_FOVs          = SpcCoeff%AC%n_FOVs          
    IF ( PRESENT(Release         ) ) Release         = SpcCoeff%Release         
    IF ( PRESENT(Version         ) ) Version         = SpcCoeff%Version         
    IF ( PRESENT(Sensor_Id       ) ) Sensor_Id       = SpcCoeff%Sensor_Id       
    IF ( PRESENT(WMO_Satellite_Id) ) WMO_Satellite_Id= SpcCoeff%WMO_Satellite_Id
    IF ( PRESENT(WMO_Sensor_Id   ) ) WMO_Sensor_Id   = SpcCoeff%WMO_Sensor_Id   

    IF ( PRESENT(Wavenumber) ) THEN
       allocate(Wavenumber(size(SpcCoeff%Wavenumber)))
       Wavenumber = SpcCoeff%Wavenumber
    else
       write(0,*) 'Wavenumber is not present'
       stop 2
    END IF

    IF ( PRESENT(Polarization) ) THEN
       allocate(Polarization(size(SpcCoeff%Polarization)))
       Polarization = SpcCoeff%Polarization
    else
       write(0,*) 'Polarization is not present'
       stop 2
    END IF

    Error_Status=Destroy_SpcCoeff(SpcCoeff,RCS_Id=RCS_Id,Message_Log=Message_Log)

    return
    
666 continue  ! Error handling
    
    IF ( File_Open( Filename ) ) THEN
       CLOSE( FileID, IOSTAT=IO_Status )
       IF ( IO_Status /= 0 ) THEN
          WRITE( Close_Message,'("; Error closing input file during error cleanup. IOSTAT=",i0)') &
               IO_Status
          Message = TRIM(Message)//TRIM(Close_Message)
       END IF
    END IF
    ! Set error status and print error message
    Error_Status = FAILURE
    CALL Display_Message( ROUTINE_NAME,TRIM(Message), &
         Error_Status,Message_Log=Message_Log )

    Ignore=Destroy_SpcCoeff(SpcCoeff,RCS_Id=RCS_Id,Message_Log=Message_Log)
    
    END FUNCTION Fancy_Inquire_SpcCoeff_Binary
end module fancy_module

PROGRAM SpcCoeff_Inspect
  USE Type_Kinds           , ONLY: Long, Double
  USE fancy_module
  USE SensorInfo_Parameters
  implicit none
  integer :: nargs,arglen,iarg,ierr,i
  character(len=:), allocatable :: filename

  integer :: n_Channels,n_FOVs, Release, Version, WMO_Satellite_Id, WMO_Sensor_Id
  character(1000) :: Sensor_Id, RCS_Id
  REAL(Double) , POINTER :: Wavenumber(:)
  INTEGER(Long), POINTER :: Polarization(:)
  REAL(Double) :: wn,c

  c=299792458.
  
10 format('error: ',A,': Inquire_SpcCoeff_Binary failed with ierr=',I0)
20 format(A,': ')
30 format('   Has ',I0,' channels, ',I0,' FOVs, release=',I0,' version=',I0)
40 format('   WMO satellite id=',I0,' WMO sensor id=',I0)
50 format('   Sensor "',A,'"')
60 format('   RCS id "',A,'"')
70 format('   Band ',I3,': wavenumber   = ',F10.2,"  m-1")
71 format('   Band ',I3,': wavelength   = ',F10.4,'  um')
72 format('   Band ',I3,': frequency    = ',F10.4,'  GHz')
80 format('   Band ',I3,': polarization = ',A)
81 format('   Band ',I3,': polarization = unknown value 'I0)
90 format('   Band ',I3,': wavenumber=',F0.5,' m-1, wavelength=',F0.3,' um, freq=',F0.2,' Ghz, ',A)
91 format('   Band ',I3,': wavenumber=',F0.5,' m-1, wavelength=',F0.3,' um, freq=',F0.2,' Ghz, unknown polarization ',I0)
99 format(A)  
  nargs=command_argument_count()
  argloop: do iarg=1,nargs
     call get_command_argument(iarg,length=arglen)
     if(allocated(filename)) deallocate(filename)
     allocate(character(len=arglen) :: filename)
     call get_command_argument(iarg,filename)

     if(associated(Wavenumber)) then
        deallocate(Wavenumber)
        nullify(Wavenumber)
     end if
     
     if(associated(Polarization)) then
        deallocate(Polarization)
        nullify(Polarization)
     end if
     
     ierr=Fancy_Inquire_SpcCoeff_Binary(&
          filename,n_Channels,n_FOVs,Release,Version,&
          Sensor_Id,WMO_Satellite_Id,WMO_Sensor_Id,RCS_Id,&
          Wavenumber=Wavenumber,Polarization=Polarization)
     if(ierr/=0) then
        write(0,10) filename,ierr
        cycle argloop
     endif
     
     print 20,filename
     print 30,n_Channels,n_FOVs,Release, Version
     print 40,WMO_Satellite_Id,WMO_Sensor_Id
     print 50,trim(Sensor_Id)
     print 60,trim(RCS_Id)

     do i=1,size(Wavenumber)
        wn=Wavenumber(i)*100.0
        !print 99,' '
        !print 70,i,wn
        !print 71,i,1.e6/wn
        !print 72,i,c*1e-9*wn
        select case(Polarization(i))
        case(INVALID_POLARIZATION)    ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'(invalid)'
        case(UNPOLARIZED)             ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'Unpolarized'
        !case(INTENSITY)               ; print 90,i,'Intensity' = unpolarized ! duplicate of Unpolarized
        !case(FIRST_STOKES_COMPONENT)  ; print 90,i,'1st Stokes' = unpolarized ! duplicate of Unpolarized
        case(SECOND_STOKES_COMPONENT) ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'2nd Stokes'
        case(THIRD_STOKES_COMPONENT)  ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'3rd Stokes'
        case(FOURTH_STOKES_COMPONENT) ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'4th Stokes'
        case(VL_POLARIZATION)         ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'V. Lin Pol'
        case(HL_POLARIZATION)         ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'H. Lin Pol'
        case(plus45L_POLARIZATION)    ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'+45 L Pol'
        case(minus45L_POLARIZATION)   ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'-45 L Pol'
        case(VL_MIXED_POLARIZATION)   ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'V. Mixed Pol'
        case(HL_MIXED_POLARIZATION)   ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'H. Mixed Pol'
        case(RC_POLARIZATION)         ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'R. Circ. Pol'
        case(LC_POLARIZATION)         ; print 90,i,wn,1.e6/wn,c*1e-9*wn,'L. Circ. Pol'
        case default                  ; print 91,i,wn,1.e6/wn,c*1e-9*wn,Polarization(i)
        end select
     end do

  end do argloop

CONTAINS

end PROGRAM SpcCoeff_Inspect
