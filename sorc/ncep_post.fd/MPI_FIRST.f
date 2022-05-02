!> @file
!
!> SUBPROGRAM:    MPI_FIRST   SET UP MESSGAE PASSING INFO
!!   PRGRMMR: TUCCILLO        ORG: IBM
!!
!! ABSTRACT:
!!     SETS UP MESSAGE PASSING INFO
!!
!! PROGRAM HISTORY LOG:
!!   14-12-01   WM LEWIS: ADDED ADDNL VARIABLES FOR SAT OUTPUT
!!   00-01-06  TUCCILLO - ORIGINAL
!!   01-10-25  H CHUANG - MODIFIED TO PROCESS HYBRID MODEL OUTPUT
!!   02-06-19  MIKE BALDWIN - WRF VERSION
!!   11-12-16  SARAH LU - MODIFIED TO INITIALIZE AEROSOL FIELDS
!!   12-01-07  SARAH LU - MODIFIED TO INITIALIZE AIR DENSITY/LAYER THICKNESS
!!   21-07-07  JESSE MENG - 2D DECOMPOSITION
!!
!! USAGE:    CALL MPI_FIRST
!!   INPUT ARGUMENT LIST:
!!
!!   OUTPUT ARGUMENT LIST:
!!
!!   OUTPUT FILES:
!!     STDOUT  - RUN TIME STANDARD OUT.
!!
!!   SUBPROGRAMS CALLED:
!!       PARA_RANGE
!!     UTILITIES:
!!       NONE
!!     LIBRARY:
!!       COMMON - CTLBLK.comm
!!
!!   ATTRIBUTES:
!!     LANGUAGE: FORTRAN
!!     MACHINE : IBM RS/6000 SP
!!
!@PROCESS NOEXTCHK
      SUBROUTINE MPI_FIRST()

!
      use vrbls4d, only: dust, salt, soot, waso, suso, pp25, pp10
      use vrbls3d, only: u, v, t, q, uh, vh, wh, pmid, pmidv, pint, alpint, zmid,      &
              zint, q2, omga, t_adj, ttnd, rswtt, rlwtt, exch_h, train, tcucn,         &
              el_pbl, cwm, f_ice, f_rain, f_rimef, qqw, qqi, qqr, qqs,qqg, qqni, qqnr, &
              extcof55, cfr, dbz, dbzr, dbzi, dbzc, mcvg, nlice, nrain, o3, vdifftt,   &
              tcucns, vdiffmois, dconvmois, sconvmois, nradtt, o3vdiff, o3prod,        &
              o3tndy, mwpv, unknown, vdiffzacce, zgdrag, cnvctummixing, vdiffmacce,    &
              mgdrag, cnvctvmmixing, ncnvctcfrac, cnvctumflx, cnvctdmflx, cnvctdetmflx,&
              cnvctzgdrag, cnvctmgdrag, icing_gfip, asy, ssa, duem, dusd, dudp,        &
              duwt, suem, susd, sudp, suwt, ocem, ocsd, ocdp, ocwt, bcem, bcsd,        &
              bcdp, bcwt, ssem, sssd, ssdp, sswt, ext, dpres, rhomid
      use vrbls2d, only: wspd10max, w_up_max, w_dn_max, w_mean, refd_max, up_heli_max, &
              prate_max, fprate_max, swupt,                                            &
              up_heli_max16, grpl_max, up_heli, up_heli16, ltg1_max, ltg2_max,         &
              up_heli_min, up_heli_min16, up_heli_max02, up_heli_min02, up_heli_max03, &
              up_heli_min03, rel_vort_max, rel_vort_max01, wspd10umax, wspd10vmax,     &
              refdm10c_max, hail_max2d, hail_maxk1, ltg3_max,rel_vort_maxhy1,          &
              nci_ltg, nca_ltg, nci_wq, nca_wq, nci_refd,                              &
              u10, v10, tshltr, qshltr, mrshltr, smstav, ssroff, bgroff,               &
              nca_refd, vegfrc, acsnow, acsnom, cmc, sst, qz0, thz0, uz0, vz0, qs, ths,&
              sno, snonc, snoavg, psfcavg, t10m, t10avg, akmsavg, akhsavg, u10max,     &
              v10max, u10h, v10h, akms, akhs, cuprec, acprec, ancprc, cuppt,           &
              rainc_bucket, rainnc_bucket, pcp_bucket, snow_bucket, qrmax, tmax,       &
              snownc, graupelnc, tsnow, qvg, qv2m, rswin, rlwin, rlwtoa, tg, sfcshx,   &
              fis, t500, cfracl, cfracm, cfrach, acfrst, acfrcv, hbot, potevp,         &
              sfclhx, htop, aswin, alwin, aswout, alwout, aswtoa, alwtoa, czen, czmean,&
              sigt4, rswout, radot, ncfrst, ncfrcv, smstot, pctsno, pshltr, th10,      &
              q10, sr, prec, subshx, snopcx, sfcuvx, sfcevp, z0, ustar, pblh, mixht,   &
              twbs, qwbs, sfcexc, grnflx, soiltb, z1000, slp, pslp, f, albedo, albase, &
              cldfra, cprate, cnvcfr, ivgtyp, hbotd, htopd, hbots, isltyp, htops,      &
              cldefi, islope, si, lspa, rswinc, vis, pd, mxsnal, epsr, sfcux,          &
              sfcvx, sfcuxi, sfcvxi, avgalbedo, avgcprate, avgprec, ptop, pbot, avgcfrach, avgcfracm,  &
              avgcfracl, avgtcdc, auvbin, auvbinc, ptopl, pbotl, ttopl, ptopm,         &
              pbotm, ttopm, ptoph, pboth, ttoph, sfcugs, sfcvgs, pblcfr, cldwork,      &
              gtaux, gtauy, mdltaux, mdltauy, runoff, maxtshltr, mintshltr,            &
              maxrhshltr, minrhshltr, dzice, alwinc, alwoutc, alwtoac, aswinc,         &
              aswoutc,aswtoac, aswintoa, smcwlt, suntime, fieldcapa, avisbeamswin,     &
              avisdiffswin, airbeamswin, airdiffswin, snowfall, dusmass, ducmass,      &
              dusmass25, susmass, sucmass, susmass25, sucmass25, ocsmass, occmass,     &
              ocsmass25, occmass25, bcsmass, bccmass, bcsmass25, bccmass25,            &
              sssmass, sscmass, sssmass25, sscmass25, ducmass25,                       &
              dustcb, sscb, bccb, occb, sulfcb, dustallcb, ssallcb,dustpm,sspm, pp25cb,&
              pp10cb, ti
      use soil, only:  smc, stc, sh2o, sldpth, rtdpth, sllevel
      use masks, only: htm, vtm, hbm2, sm, sice, lmh, gdlat, gdlon, dx, dy, lmv
      use ctlblk_mod, only: me, num_procs, jm, jsta, jend, jsta_m, jsta_m2,ista,iend , &
              jend_m, jend_m2, iup, idn, icnt, im, idsp, jsta_2l, jend_2u,idsp2,icnt2, &
              jvend_2u, lm, lp1, jsta_2l, jend_2u, nsoil, nbin_du, nbin_ss,            &
              nbin_bc, nbin_oc, nbin_su,                                               &
              ISTA_M,IEND_M,ISTA_M2,IEND_M2, iSTA_M,IEND_M,ISTA_M2,IEND_M2,            &
              ileft,iright,ileftb,irightb,ibsize,ibsum, isxa,iexa,jsxa,jexa,           &
              icoords,ibcoords,bufs,ibufs, rbufs, rcoords,rbcoords,                    &  
              ISTA_2L, IEND_2U,IVEND_2U,numx,MODELNAME   

!
!     use params_mod
!- - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - - - - 
      implicit none
     
      include 'mpif.h'
!
      integer ierr,i,jsx,jex,isx,iex,j
      integer size,ubound,lbound
      integer isumm,isum ,ii,jj,isumm2
      integer , allocatable :: ibuff(:)
      real    , allocatable :: rbuff(:)
      integer, allocatable :: ipole(:),ipoles(:,:)
      real   , allocatable :: rpole(:),rpoles(:,:)
     
      isumm=0
      isumm2=0

      if ( me == 0 ) then
        write(0,*) ' NUM_PROCS,NUMX,NUMY = ',num_procs,numx,num_procs/numx
      end if

      if ( num_procs > 1024 ) then
         print *, ' too many MPI tasks, max is 1024, stopping'
         call mpi_abort(MPI_COMM_WORLD,1,ierr)
         stop
      end if
 
!     error check
 
      if ( num_procs > JM/2 ) then
         print *, ' too many MPI tasks, max is ',jm/2,' stopping'
         call mpi_abort(MPI_COMM_WORLD,1,ierr)
         stop
      end if
 
!     global loop ranges
!
!  para_range2 supports a 2D decomposition.  The rest of the post
!  supports 1D still and the call here is the special case where each
!  processor gets all of the longitudes in the latitude 1D subdomain
!  jsta:jend.  The X decomposition will be specified by the third
!  argument (currently 1) and the Y decoposition will be specified by
!  the fourth argument (currently all of the ranks)   When X is
!  subdivided the third and fourth arguments will have to be integral
!  factors of num_procs 

      call para_range2(im,jm,numx,num_procs/numx,me,ista,iend,jsta,jend)

      jsta_m  = jsta
      jsta_m2 = jsta
      jend_m  = jend
      jend_m2 = jend
      ista_m  = ista
      ista_m2 = ista
      iend_m  = iend
      iend_m2 = iend

      if (me<numx)then
        jsta_m=2
        jsta_m2=3
      end if

      if(mod(me,numx)==0)then
        ista_m=2
        ista_m2=3
      end if

      if (me>=(num_procs-numx))then
        jend_m=jm-1
        jend_m2=jm-2
      end if

      if(mod(me+1,numx)==0)then
        iend_m=im-1
        iend_m2=im-2
      end if

 102  format(6i10,a20)

! 
      if ( me == 0 ) then
        idn = MPI_PROC_NULL
      end if
      if ( me == num_procs - 1 ) then
        iup = MPI_PROC_NULL
      end if
!
! GWV.  Array of i/j coordinates for bookkeeping tests.  Not used in
! calculations but to check if scatter,gather, and exchanges are doing as
! expected. Both real and integer arrays are sent.  Integer will be needed
! for very large domains because real mantissas overflow  and both coordinates'
! information can't be packed into a real mantisa.  Real is easier to use
! because the datatype is the same as for actual data
  
      allocate(icoords(im,jm))
      allocate(rcoords(im,jm))
      allocate(ibuff(im*jm))
      allocate(rbuff(im*jm))
      do j=1,jm
        do i=1,im
          icoords(i,j)=10000*I+j  ! both I and J information is in each element
          rcoords(i,j)=4000*i+j   ! both I and J information is in each element but it overflows for large I  I to 3600 is safe
        end do
      end do

! end COORDS test

! counts, disps for gatherv and scatterv

      isum=1
      allocate(isxa(0:num_procs-1) )
      allocate(jsxa(0:num_procs-1) )
      allocate(iexa(0:num_procs-1) )
      allocate(jexa(0:num_procs-1) )
      do i = 0, num_procs - 1
        call para_range2(im,jm,numx,num_procs/numx,i,isx,iex,jsx,jex) 
        icnt(i) = ((jex-jsx)+1)*((iex-isx)+1)
        isxa(i)=isx
        iexa(i)=iex
        jsxa(i)=jsx
        jexa(i)=jex
         
        idsp(i)=isumm
        isumm=isumm+icnt(i)                       
        if(jsx .eq. 1 .or. jex .eq. jm)  then
          icnt2(i) = (iex-isx+1)
        else
          icnt2(i)=0
        endif  
        idsp2(i)=isumm2
        if(jsx .eq. 1 .or. jex .eq. jm)  isumm2=isumm2+(iex-isx+1)

! GWV  Create send buffer for scatter.  This is now needed because we are no
! longer sending contiguous slices of the im,jm full state arrays to the
! processors with scatter.   Instead we are sending a slice of I and a slice of J
! and so have to reshape the send buffer below to make it contiguous groups of
! isx:iex,jsx:jex arrays

        do jj=jsx,jex
          do ii=isx,iex
            ibuff(isum)=icoords(ii,jj)
            rbuff(isum)=rcoords(ii,jj)
            isum=isum+1
          end do
        end do
            
      end do ! enddo of num_procs 
!
!     extraction limits -- set to two rows    
!
      jsta_2l = max(jsta - 2,  1 )
      jend_2u = min(jend + 2, jm )
      if(modelname=='GFS') then
        ista_2l=max(ista-2,0)
        iend_2u=min(iend+2,im+1)
      else
        ista_2l=max(ista-2,1)
        iend_2u=min(iend+2,im)
      endif

! special for c-grid v
      jvend_2u = min(jend + 2, jm+1 )
!
!  NEW neighbors

      ileft = me - 1
      iright = me + 1
      iup=MPI_PROC_NULL
      idn=MPI_PROC_NULL

      if(mod(me,numx) .eq. 0) print *,' LEFT POINT',me
      if(mod(me+1,numx) .eq. 0) print *,' RIGHT  POINT',me
      if(mod(me,numx) .eq. 0) ileft=MPI_PROC_NULL
      if(mod(me,numx) .eq. 0) ileftb=me+numx-1
      if(mod(me+1,numx) .eq. 0 .or. me .eq. num_procs-1)  iright=MPI_PROC_NULL
      if(mod(me+1,numx) .eq. 0 .or. me .eq. num_procs-1)  irightb=me-numx+1
      if(me .ge. numx) idn=me-numx
      if(me+1  .le. num_procs-numx) iup=me+numx

      print 102,me,ileft,iright,iup,idn,num_procs,'GWVX BOUNDS'

!     allocate arrays

      ibsize = ( (iend-ista) +1) *  ( (jend-jsta)+1)
      allocate(ibcoords(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(rbcoords(ista_2l:iend_2u,jsta_2l:jend_2u))
      allocate(ibufs(ibsize))
      allocate(rbufs(ibsize))
      call mpi_scatterv(ibuff,icnt,idsp,mpi_integer &
                    ,ibufs,icnt(me),mpi_integer ,0,MPI_COMM_WORLD,j)
      call mpi_scatterv(rbuff,icnt,idsp,mpi_real &
                    ,rbufs,icnt(me),mpi_real ,0,MPI_COMM_WORLD,j)

!
!GWV   reshape the receive subdomain

      isum=1
      do j=jsta,jend
        do i=ista,iend
          ibcoords(i,j)=ibufs(isum)
          rbcoords(i,j)=rbufs(isum)
          isum=isum+1
        end do
      end do

!GWV  end reshape
      do j=jsta,jend
        do i=ista,iend
          ii=ibcoords(i,j)/10000
          jj=ibcoords( i,j)-(ii*10000)
          if(ii .ne. i .or. jj .ne. j) then
            print *,i,j,ii,jj,ibcoords(i,j),' GWVX FAIL '
          else
            continue
          endif
         end do
      end do

      allocate(ipoles(im,2),ipole(ista:iend))
      allocate(rpoles(im,2),rpole(ista:iend))
      ipole=9900000
      ipoles=-999999999
              
      do i=ista,iend
        if(me .lt. num_procs/2. .and. jsta_2l .le. 1 .and. icnt2(me) .gt. 0) ipole(i)=ibcoords(i,1)
        if(me .lt. num_procs/2. .and. jsta_2l .le. 1 .and. icnt2(me) .gt. 0) rpole(i)=rbcoords(i,1)
        if(me .gt. num_procs/2. .and. jend_2u .ge. jm .and. icnt2(me) .gt. 0) ipole(i)=ibcoords(i,jm)
        if(me .gt. num_procs/2. .and. jend_2u .ge. jm .and. icnt2(me) .gt. 0) rpole(i)=rbcoords(i,jm)

! check code to be removed upon debugging
        if(me .lt. num_procs/2. .and. jsx .eq. 1) then
          continue
        endif
        if(me .gt. num_procs/2. .and. jend_2u  .ge. jm) then                       
          continue
        endif
      end do ! end check code

!  test pole gather
      print 105,' GWVX GATHER DISP ',icnt2(me),idsp2(me),me
 105  format(a30,3i12)

      call mpi_gatherv(ipole(ista),icnt2(me),MPI_INTEGER,  ipoles,icnt2,idsp2,MPI_INTEGER,0,MPI_COMM_WORLD, ierr ) 
      call mpi_gatherv(rpole(ista),icnt2(me),MPI_REAL   ,  rpoles,icnt2,idsp2,MPI_REAL   ,0,MPI_COMM_WORLD, ierr ) 

      if(me .eq. 0) then
        do j=1,2
          do i=1,im
            ii=rpoles(i,j)/4000
            jj=rpoles(i,j) -ii*4000
            if(ii .ne. i .or.  jj .ne. 1 .and. jj .ne. jm ) then
               write(0,169)' IPOLES BAD POINT',rpoles(i,j),ii,i,jj,' jm= ',jm
            else
              continue
!             write(0,169)'  IPOLES GOOD POINT',rpoles(i,j),ii,i,jj,' jm= ',jm
            endif
          end do
        end do
      endif

 107    format(a20,10i10)
 169    format(a25,f20.1,3i10,a10,4i10)
!
      print *, ' me, jsta_2l, jend_2u = ',me,jsta_2l, jend_2u,  &
               'jvend_2u=',jvend_2u,'im=',im,'jm=',jm,'lm=',lm, &
               'lp1=',lp1
      write(0,'(A,5I10)') 'MPI_FIRST me,jsta,jend,ista,iend,=',me,jsta,jend,ista,iend

      end

!      subroutine sub(a)
!         return
!             end



      subroutine fullpole(a,rpoles)

      use ctlblk_mod, only: num_procs, jend, iup, jsta, idn, mpi_comm_comp, im,MODELNAME,numx,&
          icoords,ibcoords,rbcoords,bufs,ibufs,me, &  
              jsta_2l,jend_2u,ileft,iright,ista_2l,iend_2u,ista,iend,jm,icnt2,idsp2
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      include 'mpif.h'
!
      real,intent(inout) :: a ( ista_2l:iend_2u,jsta_2l:jend_2u ),rpoles(im,2)
      real, allocatable ::  rpole(:)

      integer status(MPI_STATUS_SIZE)
      integer ierr
      integer size,ubound,lbound
      integer i,ii,jj, ibl,ibu,jbl,jbu,icc,jcc 
      integer ifirst
      data ifirst/0/
      integer iwest,ieast
      data iwest,ieast/0,0/
      allocate(rpole(ista:iend))

      do i=ista,iend
        if(me .lt. num_procs/2. .and. jsta_2l .le. 1  .and. icnt2(me) .gt. 0) rpole(i)=a(i,1)
        if(me .ge. num_procs/2. .and. jend_2u .ge. jm .and. icnt2(me) .gt. 0) rpole(i)=a(i,jm)
      end do

      call mpi_allgatherv(rpole(ista),icnt2(me),MPI_REAL,rpoles,icnt2,idsp2,MPI_REAL, MPI_COMM_COMP,ierr)

      call mpi_barrier(mpi_comm_comp,ierr)
      ifirst=1

      end

