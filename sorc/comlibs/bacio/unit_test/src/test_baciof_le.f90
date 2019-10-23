         program test_baciof_le
!
    implicit none
!
    integer,parameter :: kindr=4,kindr8=8,kindi=4
    integer,parameter :: im=20,jm=10,lm=5
    character(255) cin
    real(kindr) r1(im),r2(im),r3(im)
    real(kindr8) s1(jm),s2(jm),s3(jm)
    integer(kindi)   i1(lm),i2(lm),i3(lm)
    character(32) astr1,astr2,astr3
    integer flunit,iret,i,j
    integer lenr1,lens1,leni1,lenc1
    integer iskip,iwrite,nwrite,iread,nread
    integer(8) iskip8,iwrite8,nwrite8,iread8,nread8
!
!--- test baread/bawrite/wytre subroutines in baciof.f  
     r2=0.;s2=0;i2=0
     r3=0.;s3=0;i3=0
     do  i=1,im
       r1(i)=i*10.
     enddo

     do j=1,jm
       s1(j)=2.*j
     enddo
!
     do i=1,lm
       i1(i)=5*i
     enddo
!
     astr1='this is a test string'
!
!open file
     flunit=1990
     cin='data_baciof_le'
     call baopenwt(flunit,trim(cin),iret)
!
!bawrite
! write out 4 byte real
     iskip=0
     iwrite=4
     lenr1=im*kindr
     call bawrite(flunit,iskip,iwrite,nwrite,lenr1)
     iskip=iskip+nwrite
     iwrite=im*kindr
     call bawrite(flunit,iskip,iwrite,nwrite,r1)
     iskip=iskip+nwrite
     iwrite=4
     call bawrite(flunit,iskip,iwrite,nwrite,lenr1)
! write out 8 byte real
     iskip=iskip+nwrite
     iwrite=4
     lens1=jm*kindr8
     call bawrite(flunit,iskip,iwrite,nwrite,lens1)
     iskip=iskip+nwrite
     iwrite=jm*kindr8
     call bawrite(flunit,iskip,iwrite,nwrite,s1)
     iskip=iskip+nwrite
     iwrite=4
     call bawrite(flunit,iskip,iwrite,nwrite,lens1)
! write out interger
     iskip=iskip+nwrite
     iwrite=4
     leni1=lm*kindi
     call bawrite(flunit,iskip,iwrite,nwrite,leni1)
     iskip=iskip+nwrite
     iwrite=lm*kindi
     call bawrite(flunit,iskip,iwrite,nwrite,i1)
!     print *,'aft bawrite,r1,nwrite=',nwrite ,'iwrite=',iwrite,  &
!       'i1=',i1(1:5)
     iskip=iskip+nwrite
     iwrite=4
     call bawrite(flunit,iskip,iwrite,nwrite,leni1)
! write out character string
     iskip=iskip+nwrite
     iwrite=4
     lenc1=len(astr1)
     call bawrite(flunit,iskip,iwrite,nwrite,lenc1)
     iskip=iskip+nwrite
     iwrite=len(astr1)
     call bawrite(flunit,iskip,iwrite,nwrite,astr1)
!     print *,'aft bawrite,astring,nwrite=',nwrite ,'iwrite=',iwrite,  &
!       'astr1=',astr1,'len=',len(astr1)
     iskip=iskip+nwrite
     iwrite=4
     call bawrite(flunit,iskip,iwrite,nwrite,lenc1)
!
!
!close file
!
     call baclose(flunit,iret)
!
! test reading
     call baopenr(flunit,trim(cin),iret)
!read out r4
     iskip=0
     iread=4
     call baread(flunit,iskip,iread,nread,lenr1)
     iskip=iskip+nread
     iread=lenr1
     call baread(flunit,iskip,iread,nread,r2)
     iskip=iskip+nread
     iread=4
     call baread(flunit,iskip,iread,nread,lenr1)
     if(maxval(r1-r2)/=0.or.minval(r1-r2)/=0) then
       print *,'Wrong in test_baciof_le, r4'
       stop
     endif
!
!read out r8
     iskip=iskip+nread
     iread=4
     call baread(flunit,iskip,iread,nread,lens1)
     iskip=iskip+nread
     iread=lens1
     call baread(flunit,iskip,iread,nread,s2)
     iskip=iskip+nread
     iread=4
     call baread(flunit,iskip,iread,nread,lens1)
     if(maxval(s1-s2)/=0.or.minval(s1-s2)/=0) then
       print *,'Wrong in test_baciof_le, r8'
       stop
     endif
!
!read int
     iskip=iskip+nread
     iread=4
     call baread(flunit,iskip,iread,nread,leni1)
     iskip=iskip+nread
     iread=leni1
     call baread(flunit,iskip,iread,nread,i2)
!     print *,'aft baread,i1,nread=',nread ,'iread=',iread,     &
!      'leni1=',leni1,'i2=',i2(1:5)
     iskip=iskip+nread
     iread=4
     call baread(flunit,iskip,iread,nread,leni1)
     if(maxval(i1-i2)/=0.or.minval(i1-i2)/=0) then
       print *,'Wrong in test_baciof_le, i4'
       stop
     endif
!
!read string
     iskip=iskip+nread
     iread=4
     call baread(flunit,iskip,iread,nread,lenc1)
     iskip=iskip+nread
     iread=lenc1
     call baread(flunit,iskip,iread,nread,astr2)
!     print *,'aft baread,str,nread=',nread ,'iread=',iread,     &
!      'lenc1=',lenc1,'astr2=',astr2,'len(astr2)=',len(astr2)
     iskip=iskip+nread
     iread=4
     call baread(flunit,iskip,iread,nread,lenc1)
     if(astr2/=astr1 ) then
       print *,'Wrong in test_baciof, string'
       stop
     endif
!
      call baclose(flunit,iret)
!
!test with Fortran read
!
     open (flunit,file=trim(cin),form='unformatted')
     read(flunit)r3
     read(flunit)s3
     read(flunit)i3
     read(flunit)astr3
     close(flunit)
!
     if(maxval(r1-r3)/=0.or.minval(r1-r3)/=0) then
       print *,'Wrong in test_baciof, Fortran r4'
       stop
     endif
     if(maxval(s1-s3)/=0.or.minval(s1-s3)/=0) then
       print *,'Wrong in test_baciof, Fortran r8'
       stop
     endif
     if(maxval(i1-i3)/=0.or.minval(i1-i3)/=0) then
       print *,'Wrong in test_baciof, Fortran i4'
       stop
     endif
     if(astr3/=astr1) then
       print *,'Wrong in test_baciof, Fortran string'
       stop
     endif
!
     print *,'BACIO: unit test test_bacio_le ends normally'

!
     end program  test_baciof_le
