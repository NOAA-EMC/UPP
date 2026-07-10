!> @file
!> @brief Subroutine that processes fields for AI models
!
!> This routine processes fields for AI models.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2026-06-26 | Wen Meng | initial  
!>
!> @author Wen Meng @date 2026-06-26
      subroutine process_ai
    
      use mpi, only: mpi_wtime
      use ctlblk_mod, only: ista, iend, jsta, jend, datapd, fld_info, &
              spval, grib, me, im, jm, lm, cfld, ntlfld, spl, tprec, ifmin, &
              ifhr, ista_2l, iend_2u, jsta_2l, jend_2u
      use vrbls3d, only: t, zmid, q, uh, vh, q, omga
      use vrbls2d, only: tshltr, u10, v10, slp, acprec
      use rqstfld_mod, only: iget, lvls, id, iavblfld, lvlsxml
      use upp_physics, only: calvor

!----------------------------------------------------------------------------------------------
      implicit none
      integer i, j, l, itprec, ifincr
      logical log1
      real, allocatable :: grid1(:,:),grid2(:,:),grid3(:,:)

      allocate(grid1(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(grid2(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(grid3(ista_2l:iend_2u,jsta_2l:jend_2u))

      cfld=0
      grid1=spval
      grid2=spval
      grid3=spval

!-----SHELTER LEVEL TEMPERATURE
      if(iget(106)>0) then
        do j=jsta,jend
          do i=ista,iend
            grid1(i,j)=tshltr(i,j)
          enddo
        enddo
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(106))
          datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
        endif
      endif

!-----ANEMOMETER LEVEL U WIND AND/OR V WIND.
      if((iget(064)>0).or.(iget(065)>0)) then
        do j=jsta,jend
          do i=ista,iend
            grid1(i,j)=u10(i,j)
            grid2(i,j)=v10(i,j)
          enddo
        enddo
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(064))
          datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
          cfld=cfld+1
          fld_info(cfld)%ifld=IAVBLFLD(IGET(065))
          datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid2(ista:iend,jsta:jend)
        endif
      endif

!-----SEA LEVEL PRESSURE
      if(iget(105)>0) then
        do j=jsta,jend
          do i=ista,iend
            grid1(i,j)=slp(i,j)
          enddo
        enddo
        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=iavblfld(iget(105))
          datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
        endif
      endif      

!-----ACCUMULATED TOTAL PRECIPITATION.
      if(iget(087)>0) then
        do j=jsta,jend
          do i=ista,iend
            grid1(i,j)=acprec(i,j)
          enddo
        enddo
        id(1:25)=0
        itprec=nint(tprec)
        if(itprec/=0)then
          ifincr=mod(ifhr,itprec)
          if(ifmin>=1)ifincr=mod(ifhr*60+ifmin,itprec*60)
        else
          ifincr=0
        endif
        id(18)=0
        id(19)=0
        if(ifmin>=1)id(19)=ifhr*60+ifmin-ifincr
        id(20)=4
        if(ifincr==0) then
          id(18)=ifhr-itprec
        else
          id(18)=ifhr-ifincr
          if(ifmin>=1)id(18)=ifhr*60+ifmin-ifincr
        endif

        if(grib=='grib2') then
          cfld=cfld+1
          fld_info(cfld)%ifld=iavblfld(iget(087))
          fld_info(cfld)%ntrange=1
          fld_info(cfld)%tinvstat=ifhr-id(18)
          datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
        endif
      endif

      do l=1,lm

!-----TEMPERATURE on PRESSURE LEVEL
      if(iget(013)>0) then
        if(lvls(l,iget(013))>0) then
          do j=jsta,jend
            do i=ista,iend
              grid1(i,j)=t(i,j,l)
            enddo
          enddo
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=iavblfld(iget(013))
            fld_info(cfld)%lvl=lvlsxml(l,iget(013))
            datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
          endif
        endif
      endif

!-----GEOPOTENTIAL on PRESSURE LEVEL
      if(iget(012)>0) then
        if(lvls(l,iget(012))>0) then
          do j=jsta,jend
            do i=ista,iend
              grid1(i,j)=zmid(i,j,l)
            enddo
          enddo
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=iavblfld(iget(012))
            fld_info(cfld)%lvl=lvlsxml(l,iget(012))
            datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
          endif
        endif
      endif

!-----SPECIFIC HUMIDITY ON PRESSURE LEVEL
      if(iget(016)>0) then
        if(lvls(l,iget(016))>0) then
          do j=jsta,jend
            do i=ista,iend
              grid1(i,j)=q(i,j,l)
            enddo
          enddo
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=iavblfld(iget(016))
            fld_info(cfld)%lvl=lvlsxml(l,iget(016))
            datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
          endif
        endif
      endif

!-----OMEGA ON PRESSURE LEVEL
      if(iget(020)>0) then
        if(lvls(l,iget(020))>0) then
          do j=jsta,jend
            do i=ista,iend
              grid1(i,j)=omga(i,j,l)
            enddo
          enddo
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=iavblfld(iget(020))
            fld_info(cfld)%lvl=lvlsxml(l,iget(020))
            datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
          endif
        endif
      endif

!-----U AND/OR V WIND ON PRESSURE LEVEL
      if(iget(018)>0.or.iget(019)>0) then
        log1=.false.
        if(lvls(l,iget(018))>0.or.lvls(l,iget(019))>0) log1=.true.
        if(log1) then
          do j=jsta,jend
            do i=ista,iend
              grid1(i,j)=uh(i,j,l)
              grid2(i,j)=vh(i,j,l)
            enddo
          enddo
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=iavblfld(iget(018))
            fld_info(cfld)%lvl=lvlsxml(l,iget(018))
            datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid1(ista:iend,jsta:jend)
            cfld=cfld+1
            fld_info(cfld)%ifld=iavblfld(iget(019))
            fld_info(cfld)%lvl=lvlsxml(l,iget(019))
            datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid2(ista:iend,jsta:jend)
          endif
        endif
      endif

!-----ABSOLUTE VORTICITY
      if(iget(021)>0) then
        if(lvls(l,iget(021))>0) then
          do j=jsta,jend
            do i=ista,iend
              grid1(i,j)=uh(i,j,l)
              grid2(i,j)=vh(i,j,l)
              grid3(i,j)=spval
            enddo
          enddo
          call calvor(grid1,grid2,grid3)
          if(grib=='grib2') then
            cfld=cfld+1
            fld_info(cfld)%ifld=iavblfld(iget(021))
            fld_info(cfld)%lvl=lvlsxml(l,iget(021))
            datapd(1:iend-ista+1,1:jend-jsta+1,cfld)=grid3(ista:iend,jsta:jend)
          endif
        endif
      endif


      enddo !end l=1,lm

      ntlfld=cfld

      deallocate(grid1)
      deallocate(grid2)
      deallocate(grid3)
      

      return
      end
