    program test_bafrio
!
    implicit none
!
    integer,parameter :: im=20,jm=10,lm=5,kindr=4,kindr8=8
    character(255) cin
    real(kindr) r1(im),r2(im),r3(im)
    real(kindr8) s1(jm),s2(jm),s3(jm)
    integer i1(lm),i2(lm),i3(lm)
    character(32) astr1,astr2,astr3
    integer flunit,iret,i,j
    integer*4 iskip,iwrite,nwrite,iread,nread
    integer*8 iskip8,iwrite8,nwrite8,iread8,nread8
    character(16) machine_endian
    logical do_byteswap

! set local vars     
     r2=0.;s2=0.;r3=0.;s3=0.
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
!test write
     flunit=990
     cin='data_bafrio_be'
!
! get machine endian
     do_byteswap=.false.
     call chk_endianc(machine_endian)
! decide if need to do byteswap
     if(trim(machine_endian)=='big_endian') then
       if(flunit>=1000.and.flunit<=1999) do_byteswap=.true.
     elseif(trim(machine_endian)=='little_endian') then
       if(flunit<=999) then
          do_byteswap=.true.
       endif
     endif
     print *,'in test_baciof,machine_endian=',trim(machine_endian),  &
       'do_byteswap=',do_byteswap

     call baopenwt(flunit,trim(cin),iret)
!
     iskip=0
     iwrite=im*kindr
     if(do_byteswap) call byteswap(r1,kindr,im)
     call bafrwrite(flunit,iskip,iwrite,nwrite,r1)
     if(nwrite.lt.iwrite) then
       print *,'aft bawrite,r1,nwrite=',nwrite ,'iwrite=',iwrite
       stop
     endif
     if(do_byteswap) call byteswap(r1,kindr,im)

     iskip=iskip+nwrite
     iwrite=jm*kindr8
     if(do_byteswap) call byteswap(s1,kindr8,jm)
     call bafrwrite(flunit,iskip,iwrite,nwrite,s1)
     if(nwrite.lt.iwrite) then
       print *,'aft bawrite,s1,nwrite=',nwrite ,'iwrite=',iwrite
       stop
     endif
     if(do_byteswap) call byteswap(s1,kindr8,jm)
!
     iskip=iskip+nwrite
     iwrite=lm*kind(i1)
     if(do_byteswap) call byteswap(i1,kind(i1),lm)
     call bafrwrite(flunit,iskip,iwrite,nwrite,i1)
     if(nwrite.lt.iwrite) then
       print *,'aft bawrite,i1,nwrite=',nwrite ,'iwrite=',iwrite
       stop
     endif
     if(do_byteswap) call byteswap(i1,kind(i1),lm)

     iskip=iskip+nwrite
     iwrite=len(astr1)
     call bafrwrite(flunit,iskip,iwrite,nwrite,astr1)
     if(nwrite.lt.iwrite) then
       print *,'aft bawrite,astr1,nwrite=',nwrite ,'iwrite=',iwrite
       stop
     endif
!
      call baclose(flunit,iret)
!
!test read
     call baopenr(flunit,trim(cin),iret)
     iskip=0
     iread=im*kindr
     call bafrread(flunit,iskip,iread,nread,r2)
     if(do_byteswap) call byteswap(r2,kindr,im)
     if(maxval(r1-r2)/=0.or.minval(r1-r2)/=0) then
       print *,'Wrong in test_bafriof, bafrread r4'
       stop
     endif
!
     iskip=iskip+nread
     iread=jm*kindr8
     call bafrread(flunit,iskip,iread,nread,s2)
     if(do_byteswap) call byteswap(s2,kindr8,jm)
     if(maxval(s1-s2)/=0.or.minval(s1-s2)/=0) then
       print *,'Wrong in test_bafriof, bafrread, r8'
       stop
     endif
!
     iskip=iskip+nread
     iread=lm*kind(i2)
     call bafrread(flunit,iskip,iread,nread,i2)
     if(do_byteswap) call byteswap(i2,kind(i2),lm)
     if(maxval(i1-i2)/=0.or.minval(i1-i2)/=0) then
       print *,'Wrong in test_bafriof, bafrread, i4'
       stop
     endif
!
     iskip=iskip+nread
     iread=len(astr2)
     call bafrread(flunit,iskip,iread,nread,astr2)
     if(astr2/=astr1) then
       print *,'Wrong in test_bafriof, bafrread, character'
       stop
     endif

      call baclose(flunit,iret)
!
!test bareadl/bawritel
!
     flunit=991
     cin='data_bafrio_be1'

     call baopenwt(flunit,trim(cin),iret)
!
     iskip8=0
     iwrite8=im*kindr
     if(do_byteswap) call byteswap(r1,kindr,im)
     call bafrwritel(flunit,iskip8,iwrite8,nwrite8,r1)
     if(nwrite8.lt.iwrite8) then
       print *,'aft bawritel,r1,nwrite8=',nwrite8,'iwrite8=',iwrite8
       stop
     endif
     if(do_byteswap) call byteswap(r1,kindr,im)

     iskip8=iskip8+nwrite8
     iwrite8=jm*kindr8
     if(do_byteswap) call byteswap(s1,kindr8,jm)
     call bafrwritel(flunit,iskip8,iwrite8,nwrite8,s1)
     if(nwrite8.lt.iwrite8) then
       print *,'aft bawritel,r2,nwrite8=',nwrite8 ,'iwrite8=',iwrite8
       stop
     endif
     if(do_byteswap) call byteswap(s1,kindr8,jm)

     iskip8=iskip8+nwrite8
     iwrite8=lm*kind(i1)
     if(do_byteswap) call byteswap(i1,kind(i1),lm)
     call bafrwritel(flunit,iskip8,iwrite8,nwrite8,i1)
     if(nwrite8.lt.iwrite8) then
       print *,'aft bawritel,i1,nwrite8=',nwrite8 ,'iwrite8=',iwrite8
       stop
     endif
     if(do_byteswap) call byteswap(i1,kind(i1),lm)

     iskip8=iskip8+nwrite8
     iwrite8=len(astr1)
     call bafrwritel(flunit,iskip8,iwrite8,nwrite8,astr1)
     if(nwrite8.lt.iwrite8) then
       print *,'aft bawritel,astr2,nwrite8=',nwrite8 ,'iwrite8=',iwrite8
       stop
     endif
!
      call baclose(flunit,iret)
!
!test read
     r2=0.;s2=0.;i2=0;astr2=''
     call baopenr(flunit,trim(cin),iret)
     iskip8=0
     iread8=im*kindr
     call bafrreadl(flunit,iskip8,iread8,nread8,r2)
     if(do_byteswap) call byteswap(r2,kindr,im)
     if(maxval(r1-r2)/=0.or.minval(r1-r2)/=0) then
       print *,'Wrong in test_bafriof, bafrreadl r4'
       stop
     endif
!
     iskip8=iskip8+nread8
     iread8=jm*kindr8
     call bafrreadl(flunit,iskip8,iread8,nread8,s2)
     if(do_byteswap) call byteswap(s2,kindr8,jm)
     if(maxval(s1-s2)/=0.or.minval(s1-s2)/=0) then
       print *,'Wrong in test_bafriof, bafrreadl,r8'
       stop
     endif
!
     iskip8=iskip8+nread8
     iread8=lm*kind(i2)
     call bafrreadl(flunit,iskip8,iread8,nread8,i2)
     if(do_byteswap) call byteswap(i2,kind(i2),lm)
     if(maxval(i1-i2)/=0.or.minval(i1-i2)/=0) then
       print *,'Wrong in test_bafriof, bafrreadl,i4'
       stop
     endif
!
     iskip8=iskip8+nread8
     iread8=len(astr2)
     call bafrreadl(flunit,iskip8,iread8,nread8,astr2)
     if(astr2/=astr1) then
       print *,'Wrong in test_bafriof, bafrreadl,character'
       stop
     endif

!
      call baclose(flunit,iret)
!
!test fortran read
     open(flunit,file=trim(cin),form='unformatted')
     read(flunit)r3
     read(flunit)s3
     read(flunit)i3
     read(flunit)astr3
     close(flunit)
     if(maxval(r1-r3)/=0.or.minval(r1-r3)/=0) then
       print *,'Wrong in test_bafriof, Fortran r4'
       stop
     endif
     if(maxval(s1-s3)/=0.or.minval(s1-s3)/=0) then
       print *,'Wrong in test_bafriof, Fortran r8'
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
     print *,'BACIO: unit test test_bafrio ends normally'
!
     end program  test_bafrio
