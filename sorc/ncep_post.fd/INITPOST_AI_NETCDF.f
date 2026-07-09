!> @file
!> @brief initpost_ai_netcdf() initializes post for AIGFS model.
!>
!> @author Wen Meng @date 2026-06-26
!>
!> This routine initializes constants and read variables from AI models
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 2026-06-26 | Wen Meng | Initial.
!>
!> @author Wen Meng @date 2026-06-26
!----------------------------------------------------------------------
!> @brief initializes constants and read variables from AI models
!> @param[in] ncid integer netCDF ID
!> @param[in] idate integer inilitial date 
!----------------------------------------------------------------------
      subroutine initpost_ai_netcdf(ncid,idate)

      use netcdf
      use iso_fortran_env, only: int64
      use vrbls3d, only: t,zmid,q,uh,vh,q,omga
      use vrbls2d, only: tshltr,u10,v10,slp,acprec,f
      use ctlblk_mod, only: me,im,jm,lm,ista,iend,jsta,jend,ista_2l,iend_2u,jsta_2l, jend_2u, &
              spval, gdsdegr, idat, sdat, ihrst, imin, ifhr, ifmin, &
              tprec, tclod, trdlw, trdsw, tsrfc, tmaxmin, global, &
              ista_2l, iend_2u, jsta_2l, jend_2u, iend_m
      use gridspec_mod, only: maptype, gridtype, latstart, latlast, lonstart, lonlast, &
              dxval, dyval
      use masks, only: gdlat, gdlon
      use params_mod, only: dtr
!----------------------------------------------------------------------------------------------
      implicit none
      INCLUDE "mpif.h"

      integer,intent(in) :: idate(8)
      integer ncid, status, varid, numdims
      integer i, j, l
      character(len=30)  :: varname
      real, allocatable :: glat1d(:), glon1d(:)
      integer, allocatable :: plevel(:)
      logical, parameter :: debugprint = .false.
      integer jdate(8)
      integer start(1), count(1)
      real rinc(5)

      allocate(glat1d(jm),glon1d(im),plevel(lm))

      if(me==0)then
      WRITE(6,*)'INITPOST:  ENTER INITPOST_NETCDF'
      WRITE(6,*)'me=',me,  &
           'jsta_2l=',jsta_2l,'jend_2u=', &
           jend_2u,'im=',im, &
           'ista_2l=',ista_2l,'iend_2u=',iend_2u, &
           'ista=',ista,'iend=',iend, &
           'iend_m=',iend_m
      endif

!-----Set default for AI models-------------------
      maptype=0 !set as regular latlon grid
      dxval=0.25*gdsdegr !dlon as 0.25 degree
      dyval=0.25*gdsdegr !dlat as 0.25 degree
      gridtype='A'
      global=.true.


!-----read lat/lon------------
      status=nf90_inq_varid(ncid,'lon',varid)
      status=nf90_inquire_variable(ncid,varid,ndims=numdims)
      if(numdims==1)then
        status=nf90_get_var(ncid,varid,glon1d)
        do j=jsta,jend
          do i=ista,iend
            !gdlon(i,j) = real(glon1d(i),kind=4)
            gdlon(i,j)=glon1d(i)
          end do
        end do
      end if
      lonstart=nint(glon1d(1)*gdsdegr)
      lonlast=nint(glon1d(im)*gdsdegr)
      if(maptype == 0) then
        if(lonstart<0.)then
         lonstart=lonstart+360.*gdsdegr
        end if
        if(lonlast<0.)then
         lonlast=lonlast+360.*gdsdegr
        end if
      end if

      status=nf90_inq_varid(ncid,'lat',varid)
      status=nf90_inquire_variable(ncid,varid,ndims=numdims)
      if(numdims==1)then
        status=nf90_get_var(ncid,varid,glat1d)
        do j=jsta,jend
          do i=ista,iend
            !gdlat(i,j) = real(glat1d(j),kind=4)
            gdlat(i,j)=glat1d(j)
          end do
        end do
      end if
      latstart = nint(glat1d(1)*gdsdegr)
      latlast  = nint(glat1d(jm)*gdsdegr)


      status=nf90_inq_varid(ncid,'level',varid)
      status=nf90_get_var(ncid,varid,plevel)

!-----set date---------------
      sdat(1)=idate(2)
      sdat(2)=idate(3)
      sdat(3)=idate(1)
      ihrst=idate(5)
      imin=idate(6)
      rinc=0.
      ifmin=0
      !idate=(/2019,12,30,0,18,0,0,0/)
      rinc(2)=real(ifhr)
      call w3movdat(rinc,idate,jdate)
      idat(1)=jdate(2)
      idat(2)=jdate(3)
      idat(3)=jdate(1)
      idat(4)=jdate(5)
      idat(5)=jdate(6)
      if(debugprint)then
      if(me==0) print*, "Initial Date: ", idate(1:7)
      if(me==0) print*, "forecast hour: ", rinc(2) 
      if(me==0) print*, "Valid Date: ", jdate(1:7)
      endif

!-----Set buckets -------------
      tprec=6
      tclod=tprec
      trdlw=tprec
      trdsw=tprec
      tsrfc=tprec
      tmaxmin=tprec

!-----read 2D fields--------------
      varname='2m_temperature'
      call read_netcdf_2d(ncid,varname,tshltr(ista:iend,jsta:jend))

      varname='10m_u_component_of_wind'
      call read_netcdf_2d(ncid,varname,u10(ista:iend,jsta:jend))

      varname='10m_v_component_of_wind'
      call read_netcdf_2d(ncid,varname,v10(ista:iend,jsta:jend))

      varname='mean_sea_level_pressure'
      call read_netcdf_2d(ncid,varname,slp(ista:iend,jsta:jend))

      varname='total_precipitation_6hr'
      call read_netcdf_2d(ncid,varname,acprec(ista:iend,jsta:jend))

!-----read 3D fields--------------
      varname='temperature'
      call read_netcdf_3d(ncid,varname,t(ista:iend,jsta:jend,1:lm),lm)

      varname='geopotential'
      call read_netcdf_3d(ncid,varname,zmid(ista:iend,jsta:jend,1:lm),lm)

      varname='specific_humidity'
      call read_netcdf_3d(ncid,varname,q(ista:iend,jsta:jend,1:lm),lm)

      varname='u_component_of_wind'
      call read_netcdf_3d(ncid,varname,uh(ista:iend,jsta:jend,1:lm),lm)

      varname='v_component_of_wind'
      call read_netcdf_3d(ncid,varname,vh(ista:iend,jsta:jend,1:lm),lm)

      varname='vertical_velocity'
      call read_netcdf_3d(ncid,varname,omga(ista:iend,jsta:jend,1:lm),lm)

!$omp parallel do private(i,j)
      do j=jsta,jend
        do i=ista,iend
          f(i,j)=1.454441e-4*sin(gdlat(i,j)*dtr) ! 2*omeg*sin(phi)
        enddo
      enddo

      deallocate(glat1d,glon1d,plevel)

      end

!----------------------------------------------------------------------
!> @brief read_netcdf_3d() reads 3D variables from AI models
!> 
!> @param[in] ncid integer netCDF ID.
!> @param[in] varname character Variable name in netCDF file.
!> @param[in] lm integer Model levels.
!> @param[out] buf real Variable values.
!----------------------------------------------------------------------
      subroutine read_netcdf_3d(ncid,varname,buf,lm)

      use netcdf
      use ctlblk_mod, only : me,ista,iend,jsta,jend,spval
      use params_mod, only : small
      use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
      implicit none
      include "mpif.h"

      character(len=*),intent(in) :: varname
      integer,intent(in) :: ncid,lm
      real,intent(out)   :: buf(ista:iend,jsta:jend,lm)
      integer            :: varid,iret,ii,jj,i,j,l,kk
      integer            :: start(3), count(3), stride(3)
      real,parameter     :: spval_netcdf=9.99e+20
      real               :: fill_value

      iret = nf90_inq_varid(ncid,trim(varname),varid)
      if (iret /= 0) then
        if (me == 0) print*,VarName," not found -Assigned missing values"
!$omp parallel do private(i,j,l)
          do l=1,lm
            do j=jsta,jend
              do i=ista,iend
                buf(i,j,l)=spval
              enddo
            enddo
          enddo
      else
        iret = nf90_get_att(ncid,varid,"_FillValue",fill_value)
        if (iret /= 0) fill_value = spval_netcdf
        start = (/ista,jsta,1/)
        ii=iend-ista+1
        jj=jend-jsta+1
        count = (/ii,jj,lm/)
        iret = nf90_get_var(ncid,varid,buf(ista:iend,jsta:jend,1:lm),start=start,count=count)
        if (iret /= 0) then
          print*," iret /=0, Error in reading varid "
        endif
        do l=1,lm
          do j=jsta,jend
            do i=ista,iend
              if (ieee_is_nan(buf(i,j,l)) .or. ieee_is_nan(fill_value)) then
                if (ieee_is_nan(buf(i,j,l))) buf(i,j,l) = spval
              else if (abs(buf(i,j,l) - fill_value) < small) then
                buf(i,j,l) = spval
              endif
            end do
          end do
        end do
      endif


      end subroutine read_netcdf_3d

!----------------------------------------------------------------------
!> @brief read_netcdf_2d() reads 2D variables from AI models
!> 
!> @param[in] ncid integer netCDF ID.
!> @param[in] varname character Variable name in netCDF file.
!> @param[out] buf real Variable values.
!----------------------------------------------------------------------
      subroutine read_netcdf_2d(ncid,varname,buf)

      use netcdf
      use ctlblk_mod, only : me,ista,iend,jsta,jend,spval
      use params_mod, only : small
      use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
      implicit none
      include "mpif.h"

      character(len=*),intent(in) :: varname
      integer,intent(in) :: ncid
      real,intent(out)   :: buf(ista:iend,jsta:jend)
      integer            :: varid,iret,ii,jj,i,j,l,kk
      integer            :: start(2), count(2), stride(3)
      real,parameter     :: spval_netcdf=9.99e+20
      real               :: fill_value

      iret = nf90_inq_varid(ncid,trim(varname),varid)
      if (iret /= 0) then
        if (me == 0) print*,VarName," not found -Assigned missing values"
!$omp parallel do private(i,j,l)
            do j=jsta,jend
              do i=ista,iend
                buf(i,j)=spval
              enddo
            enddo
      else
        iret = nf90_get_att(ncid,varid,"_FillValue",fill_value)
        if (iret /= 0) fill_value = spval_netcdf
        start = (/ista,jsta/)
        ii=iend-ista+1
        jj=jend-jsta+1
        count = (/ii,jj/)
        iret = nf90_get_var(ncid,varid,buf(ista:iend,jsta:jend),start=start,count=count)
        if (iret /= 0) then
          print*," iret /=0, Error in reading varid "
        endif
          do j=jsta,jend
            do i=ista,iend
              if (ieee_is_nan(buf(i,j)) .or. ieee_is_nan(fill_value)) then
                if (ieee_is_nan(buf(i,j))) buf(i,j) = spval
              else if(abs(buf(i,j)-fill_value)<small) then
                buf(i,j)=spval
              endif
            end do
          end do
      endif


      end subroutine read_netcdf_2d

