! Test program for the module "XMLPARSE"
! Read a set of sample XML-files and report the
! possible failure to the screen
!
! $Id: tstparse.f90,v 1.3 2006/03/26 19:05:48 arjenmarkus Exp $
!
! Per test case read:
! - The title
! - The XML file to be
! - The expected parse data:
!   - which tags and the attributes that belong to them
!   - which data lines
!
! If the XML file was read but the actually reported data did
! not match the expected ones, then the title of the test
! case is printed on screen.
!
! So, a successful test does not produce output.
!
! The test cases are defined in a file "tstpars.inp"
! which has the following format:
! - Title line
! - A line which starts with "-----"
! - The contents of the XML file
! - Another line with "-----"
! - Lines with keywords (in order) that define what
!   tags, attributes and data to expect:
!   tag: string
!   endtag: T/F
!   attribute: A=B
!   data: string
! - Another line with "-----" to separate it from the
!   next case
! - Lines with "#" in the first column are considered
!   comment and are therefore ignored
!
program tstparse

   use xmlparse

   implicit none

   character(len=60) :: title
   character(len=80) :: string
   logical           :: has_title
   integer           :: ierr


   !
   ! Start the loop
   !
   open( 10, file = 'tstparse.inp', status = 'old' )
   open( 20, file = 'tstparse.out' )

   has_title = .false.

   do while ( ierr .eq. 0 )
      read( 10, '(a)', iostat = ierr ) string
      !
      ! End of file or some error condition
      !
      if ( ierr .ne. 0 ) then
         exit
      endif
      !
      ! Skip comments
      !
      if ( string(1:1) .eq. '#' ) then
         cycle
      endif
      !
      ! Store the title and start reading the XML data
      !
      if ( .not. has_title ) then
         title = string
         has_title = .true.
      endif
      if ( has_title .and. string(1:5) .eq. '-----' ) then
         open( 30, file = 'tstparse.xml' )
         do
            read(  10, '(a)', iostat = ierr ) string
            if ( ierr .ne. 0 ) then
               exit
            endif
            if ( string(1:5) .eq. '-----' ) then
               exit
            endif
            write( 30, '(a)' ) string
         enddo

         !
         ! Now we are ready to read the XML file
         ! and analyse the report
         !
         close( 30 )
         call check_xml_file
      endif
   enddo

contains

subroutine check_xml_file
   logical           :: mustread
   type(XML_PARSE)   :: info

   character(len=80)                      :: tag
   logical                                :: endtag
   character(len=80), dimension(1:2,1:20) :: attribs
   integer                                :: no_attribs
   character(len=200), dimension(1:100)   :: data
   integer                                :: no_data
   integer                                :: i

   call xml_open( info, 'tstparse.xml', .true. )
   call xml_options( info, report_lun = 20, report_details = .true. )

   do
      call xml_get( info, tag, endtag, attribs, no_attribs, data, no_data )
      if ( xml_error(info) ) then
         write(*,*)  'Error ', title
         write(20,*) 'Error ', title
      endif

      write( 20,* ) tag, endtag
      do i = 1,no_attribs
         write( 20,'(i3,1x,3a)') &
            i, trim(attribs(1,i)), '=', trim(attribs(2,i))
      enddo
      write( 20,'(i3,1x,3a)') (i, '>',trim(data(i)), '<', i=1,no_data)
      if ( .not. xml_ok(info) ) exit
   enddo

   do
      read(10,'(a)' ) string
      if ( string(1:5) .eq. '-----' ) then
         exit
      endif
   enddo

   call xml_close( info )

end subroutine check_xml_file
end program tstparse
