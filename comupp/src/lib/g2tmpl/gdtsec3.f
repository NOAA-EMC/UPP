      module gdtsec3
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! MODULE:    gdtsec3 
!   PRGMMR: Vuong         ORG: W/NP11    DATE: 2010-02-05
!
! ABSTRACT: This Fortran Module allows access to predefined GRIB2 Grid 
!   Definition Templates stored in a file.  The GDTs are represented by 
!   a predefined number.
!
!   At the first request, all the grid GDT entries in the file associated
!   with input Fortran file unit number, lunit, are read into a linked list
!   named gridlist.  This list is searched for the requested entry.
!
!   Users of this Fortran module should only call routines getgdtnum.
!
!   The format of the file scanned by routines in this module is as follows.
!   Each line contains one Grid entry containing eight fields, each separated
!   by space.  The fields are:
!      1) - The first character is a space
!      2) - predefined grid ID number
!      3) - Grid Definition Template number
!      4) - Number of entries in the Grid Definition Template
!      5) - Contains information read from the appropriate GRIB Grid
!      6) - A list of values for each entry in the Grid Definition Template.
!      7) - The number of entries in array ideflist i.e., number of rows
!           (or columns) for which optional grid points are defined. eg. WAFS grid
!      8) - Optional integer array containing the number of grid points contained
!           in each row (or columns).
!
!   As an example, this is the entry for the 1x1 GFS global grid 
!    3 0 19 0 41760 0 0 0 6 0 0 0 0 0 0 360 181 0 0 90000000 0 48 -90000000 359000000 1000000 1000000 0 0 0
!
!   Comments can be included in the file by specifying the symbol "#" as the
!   first character on the line.  These lines are ignored.
!
! PROGRAM HISTORY LOG:
! 2010-02-05  Vuong
!
! USAGE:    use gdtsec3
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
      implicit none
      private
!
      integer,parameter :: MAXTEMP=200

      type,private :: gdts3
          integer :: grid_num
          integer :: gdt_num
          integer :: gdt_len
          integer :: gds_num(5)
          integer,dimension(MAXTEMP) :: gridtmpl
          integer :: def_num
          integer,dimension(100) :: def_list
          type(gdts3),pointer :: next
      end type gdts3

      type(gdts3),pointer,private :: gridlist
      integer :: num_grids=0, iret

      public getgdtnum
!
      contains

      integer function readgrids(lunit)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    readgrids
!   PRGMMR: Gilbert         ORG: W/NP11    DATE: 2001-02-05
!
! ABSTRACT: This function reads the list of GDT entries in the file 
!   associated with fortran unit, lunit.  All the entries are stored in a
!   linked list called gdts3.
!
! PROGRAM HISTORY LOG:
! 2010-02-05  Vuong
!
! USAGE:    number=readgrids(lunit)
!   INPUT ARGUMENT LIST:
!     lunit   - Fortran unit number associated the the GDT file.
!
! RETURNS:  The number of Grid Definition Templates read in.
!
! REMARKS: None
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
           integer,intent(in) :: lunit

           integer,parameter :: linelen=1280
           character(len=8) :: desc
           character(len=linelen) :: cline
           integer  ient,igdtn,igdtmpl(200),igdtlen,idefnum,i,j
           integer  pos(10), igds(5),ixlen,ideflist(100)

           type(gdts3),pointer :: gtemp
           type(gdts3),pointer :: prev
           integer count,iret

           count=0
           iret=0

           !   For each line in the file....
           DO
             !  Read line into buffer
             !
             pos=0
             igds=0
             cline(1:linelen)=' '
             read(lunit,end=999,fmt='(a)') cline
             print *,'in read,count=',count,'cline=',cline
             !
             !  Skip line if commented out
             !
             if (cline(1:1).eq.'#') cycle
             !
             !  find positions of delimiters, " " (space)
             !
             do i=1,8
                pos(i)=index(cline,' ')
                cline(pos(i):pos(i))=';'
             end do
             if ( pos(1).eq.0 .or. pos(2).eq.0 .or. pos(3).eq.0 .or.    &
     &            pos(4).eq.0) cycle
             !  Read each of the four fields.
             !
             If (pos(1).eq.1) then
                ixlen=0
                pos(9)=index(cline,' ')
                cline(pos(9):pos(9))=';'
                read(cline(pos(1)+1:pos(2)-1),*) ient
                read(cline(pos(2)+1:pos(3)-1),*) igdtn
                read(cline(pos(3)+1:pos(4)-1),*) igdtlen
                read(cline(pos(4)+1:pos(5)-1),*) igds(1)
                read(cline(pos(5)+1:pos(6)-1),*) igds(2)
                read(cline(pos(6)+1:pos(7)-1),*) igds(3)
                read(cline(pos(7)+1:pos(8)-1),*) igds(4)
                read(cline(pos(8)+1:pos(9)-1),*) igds(5)
                if (ient.ge.21.and.ient.le.24) ixlen=37
                if (ient.eq.25.or.ient.eq.26)  ixlen=19
                if (ient.ge.37.and.ient.le.44) ixlen=73
                if (ient.ge.61.and.ient.le.64) ixlen=46
                read(cline(pos(9)+1:linelen),*)                       &
     &           (igdtmpl(j),j=1,igdtlen+1+ixlen)

          print *,'22 igds = ', igds(1:5)
          print *,'22 igdtmpl = ', (igdtmpl(j),j=1,igdtlen+1+ixlen)

                idefnum=igdtmpl(igdtlen+1)

               print *,' idefnum = ', idefnum
               print *,' before ixlen = ', ixlen,'ient = ', ient
               print *,' ient=',ient,'igdtn=',igdtn,'igdtlen=',igdtlen

                ideflist(1:ixlen)=igdtmpl(igdtlen+2:igdtlen+2+ixlen)

               print *,' ideflist = ', ideflist(1:ixlen)
             else
                print *,' Error Reading File : gdtsec3.template' 
                iret = 0
             end if
            !
            !  Allocate new type(gdts3) variable to store the GDT
            !
!!!             allocate(gtemp,stat=iom)
             count=count+1
             gtemp%grid_num=ient
             gtemp%gdt_num=igdtn
             gtemp%gds_num=igds
             gtemp%gdt_len=igdtlen
             gtemp%gridtmpl=igdtmpl
             gtemp%def_num=idefnum
             gtemp%def_list=ideflist
             nullify(gtemp%next)              ! defines end of linked list.
             if ( count .eq. 1 ) then
                gridlist => gtemp
             else                       ! make sure previous entry in list
                prev%next => gtemp      ! points to the new entry,
             endif
             prev => gtemp
             print *,'count=',count

           enddo

 999       readgrids=count
           return

         end function

         subroutine getgdtnum(lunit,number,igdtn,igdtlen,igds,           &
     &                        igdtmpl,idefnum,ideflist,iret)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .                                       .
! SUBPROGRAM:    getgdtnum 
!   PRGMMR: Vuong           ORG: W/NP11    DATE: 2010-02-05
!
! ABSTRACT: This subroutine searches a file referenced by fortran unit lunit
!   for a Grid Definition Template assigned to the requested number. 
!   The input file format is described at the top of this module.
!
! PROGRAM HISTORY LOG:
! 2010-01-06  Vuong   
!
! USAGE:    CALL getgdtnum(lunit,number,igdtn,igdtlen,igds,
!         &                igdtmpl,idefnum,ideflist,iret)
!   INPUT ARGUMENT LIST:
!     lunit    - Unit number of file containing Grid definitions 
!     number   - Grid number of the requested Grid definition 
!
!   OUTPUT ARGUMENT LIST:      
!     igdtn     - NN, indicating the number of the Grid Definition 
!                 Template 3.NN
!     igdtlen   - Number of entries in the Grid Definition Template
!     igds()    - Contains information read from the appropriate GRIB Grid
!                 igds(1)=Source of grid definition (see Code Table 3.0)
!                 igds(2)=Number of grid points in the defined grid.
!                 igds(3)=Number of octets needed for each
!                        additional grid points definition.
!                        Used to define number of points in 
!                        each row ( or column ) for non-regular grids.
!                        = 0, if using regular grid.
!                 igds(4)=Interpretation of list for optional points
!                            definition.  (Code Table 3.11)
!                 igds(5)=Grid Definition Template Number (Code Table 3.1)
!     igdtmpl() - An array containing the values of each entry in
!                 the Grid Definition Template.
!     idefnum   - The number of entries in array ideflist.
!                 i.e. number of rows ( or columns )
!                 for which optional grid points are defined.
!     ideflist()- Optional integer array containing the number of 
!                 grid points contained in each row (or column).
!     iret      - Error return code.
!                 0 = no error
!                 -1 = Undefined Grid number.
!                 3 = Could not read any grids from file.
!
! REMARKS: None
!
! ATTRIBUTES:
!   LANGUAGE: Fortran 90
!   MACHINE:  IBM SP
!
!$$$
           integer,intent(in) :: lunit,number
           integer,intent(out):: igdtn,igdtlen,igds(:),igdtmpl(:)
           integer,intent(out):: idefnum,ideflist(100),iret

           type(gdts3),pointer :: tempgrid

           iret=0
           igdtn=-1
           !igdtmpl=0

           !
           !  If no grids in list, try reading them from the file.
           !
           if ( num_grids.eq.0 ) then
              print *,'in gdtnum,before readgrids,lunit=',lunit
              num_grids=readgrids(lunit)
           endif

           if ( num_grids.eq.0 ) then
              iret=-1                         ! Undefined Grid number
              return
           endif

           if ( num_grids.eq.999 ) then
              iret=3                         ! problem reading file
              return
           endif

           tempgrid => gridlist

           !
           !  Search through list
           !
           do while ( associated(tempgrid) )
               if ( number .eq. tempgrid%grid_num ) then
                  igdtn=tempgrid%gdt_num
                  igdtlen=tempgrid%gdt_len
                  igds(1:5)=tempgrid%gds_num(1:5)
                  igdtmpl(1:tempgrid%gdt_len)=                              &
     &                        tempgrid%gridtmpl(1:tempgrid%gdt_len)
                  idefnum = tempgrid%def_num  ! number of entries in array ideflist.
!             Optional integer array containing number of 
!             grid points contained in each row (or column). 
                  ideflist= tempgrid%def_list 
                  print *,'find numer=',number
                  return
               else
                  tempgrid => tempgrid%next
               endif
           enddo 
 
           iret=-1
           return
 
         end subroutine

      end
