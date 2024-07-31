!> @file
!> @brief Subroutine that computes the velocity potential and
!> streamfunction from isobaric winds.
!>
!><pre>
!> This routine is based on the CFS genpsiandchi program that
!> computes velocity potential and streamfunction from the
!> isobaric wind components. The program was authored and provided
!> by Saha and H. Chuang.
!> Given the U-V wind components at P-points, this routine
!> collects the winds in the full IM,JM,LSM domain,
!> transforms them back to spectrum space and computes divergence,
!> vorticity, streamfunction and potential. The routine returns: 
!> 	  PSI: the streamfunction in global domain
!> 	  CHI: the velocity potential in global domain
!></pre>
!>   
!> @param[in] UISO real U-wind component (m/s) at all P-points.
!> @param[in] VISO real V-wind component (m/s) at all P-points.
!> @param[out] CHI real velocity potential (m^2/s) in full grid domain at all P-points.
!> @param[out] PSI real streamfunction (m^2/s) in full grid domain at all P-points
!>
!> ### Program history log:
!> Date       | Programmer   | Comments
!> -----------|--------------|---------
!> 2024-07-17 | Karina Asmar | Initial   
!> 2024-07-25 | Jesse Meng   | Fixes for mpi processes 
!>
!> @author Karina Asmar EMC/VPPPG @date 2024-07-17
!-----------------------------------------------------------------------
!> @brief Subroutine that computes velocity potential and streamfunction
!> from isobaric winds.  
!>   
!> @param[in] UISO real U-wind component (m/s) at all P-points.
!> @param[in] VISO real V-wind component (m/s) at all P-points.
!> @param[out] CHI real velocity potential (m^2/s) in full grid domain at P-points.
!> @param[out] PSI real streamfunction (m^2/s) in full grid domain at P-points
!-----------------------------------------------------------------------
      SUBROUTINE CALCHIPSI(UISO,VISO,CHI,PSI)
!
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
!     
      use gridspec_mod, only: IDRT
      use ctlblk_mod, only: ISTA, IEND, JSTA, JEND, IM, JM, LSM, ME, SPVAL, MPI_COMM_COMP,&
                            num_procs, icnt, idsp, isxa, iexa, jsxa, jexa
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none

      include 'mpif.h'
!
!     DECLARE VARIABLES.
!     
      integer :: JCAP, I, J, L, IERR
      REAL, dimension(ISTA:IEND,JSTA:JEND,LSM), intent(in) :: UISO, VISO
      REAL, dimension(IM,JM,LSM), intent(out) ::  CHI, PSI

      integer k, m
      real, allocatable :: CHI1(:),CHISUB(:),PSI1(:),PSISUB(:),COL_UWIND(:,:),COL_VWIND(:,:),  		&
      			   IN_UWIND(:,:,:),IN_VWIND(:,:,:),OUT_UWIND(:,:,:),OUT_VWIND(:,:,:), 		&
	    		   DIV(:,:,:),ZO(:,:,:),CHI_OUT(:,:,:),PSI_OUT(:,:,:)
      
!     
!***************************************************************************
!     START CALCHIPSI HERE.
!    
!     SAVE ALL P LEVELS OF U/V WINDS AT GLOBAL GRID  

      ALLOCATE(COL_UWIND(IM,JM))
      ALLOCATE(COL_VWIND(IM,JM))
      ALLOCATE(IN_UWIND(IM,JM,LSM))
      ALLOCATE(IN_VWIND(IM,JM,LSM))

      DO L=1,LSM
        CALL COLLECT_ALL(UISO(ISTA:IEND,JSTA:JEND,L),COL_UWIND)
        CALL COLLECT_ALL(VISO(ISTA:IEND,JSTA:JEND,L),COL_VWIND)
        DO J=1,JM
          DO I=1,IM
	    IN_UWIND(I,J,L)=COL_UWIND(I,J)
      	    IN_VWIND(I,J,L)=COL_VWIND(I,J)
          ENDDO
        ENDDO
      ENDDO

      DEALLOCATE(COL_UWIND)
      DEALLOCATE(COL_VWIND)

! FILL CHI/PSI VALUES
      DO L=1,LSM
        DO J=1,JM
	  DO I=1,IM
     	    CHI(I,J,L)= SPVAL
	    PSI(I,J,L)= SPVAL
     	   ENDDO
	ENDDO
      ENDDO

      IF (ME==0) THEN 

	      ALLOCATE(OUT_UWIND(IM,JM,LSM))
	      ALLOCATE(OUT_VWIND(IM,JM,LSM))
	      ALLOCATE(DIV(IM,JM,LSM))
	      ALLOCATE(ZO(IM,JM,LSM))
	      ALLOCATE(CHI_OUT(IM,JM,LSM))
	      ALLOCATE(PSI_OUT(IM,JM,LSM))
      
        ! SET MAX WAVELENGTH FOR SPECTRAL TRUNCATION
	      IF(IDRT == 0)THEN
	        JCAP = (JM-3)/2
	      ELSE
	        JCAP = JM-1
	      ENDIF

	      ! COMPUTE CHI/PSI FROM WIND VECTORS IN SPECTRAL SPACE
 	      CALL SPTRUNV(0,JCAP,IDRT,IM,						             &
  		           JM,IDRT,IM,JM,LSM,						             &
	                   0,0,0,0,							             &
	                   0,0,0,0,							             &
	                   IN_UWIND(1,1,1),IN_VWIND(1,1,1),				             &
	      	           .FALSE.,OUT_UWIND(1,1,1),OUT_VWIND(1,1,1),			  	     &
	                   .FALSE.,DIV,ZO,						             &
	                   .TRUE.,CHI_OUT(1,1,1),PSI_OUT(1,1,1))

	      DEALLOCATE(IN_UWIND)
       	      DEALLOCATE(IN_VWIND)
      	      DEALLOCATE(OUT_UWIND)
	      DEALLOCATE(OUT_VWIND)
              DEALLOCATE(DIV)
	      DEALLOCATE(ZO)

      ENDIF                             ! END OF ME=0 BLOCK

      CALL MPI_BARRIER(MPI_COMM_COMP, IERR)

      ALLOCATE(CHI1(im*jm))
      ALLOCATE(CHISUB(icnt(me)))
      ALLOCATE(PSI1(im*jm))
      ALLOCATE(PSISUB(icnt(me)))

      DO L=1,LSM

         IF (ME==0) THEN
           k=0
           DO m=0,num_procs-1
           DO J=jsxa(m),jexa(m)
           DO I=isxa(m),iexa(m)
              k=k+1
              CHI1(k)=CHI_OUT(I,J,L)
              PSI1(k)=PSI_OUT(I,J,L)
           ENDDO
           ENDDO
           ENDDO
         ENDIF
      
         CALL MPI_SCATTERV(CHI1,icnt,idsp,MPI_REAL, &
                           CHISUB,icnt(me),MPI_REAL,0,MPI_COMM_WORLD,IERR)
         CALL MPI_SCATTERV(PSI1,icnt,idsp,MPI_REAL, &
                           PSISUB,icnt(me),MPI_REAL,0,MPI_COMM_WORLD,IERR)

 
         k=0
         DO J=JSTA,JEND
         DO I=ISTA,IEND
            k=k+1
            CHI(I,J,L)=CHISUB(k)
            PSI(I,J,L)=PSISUB(k)
         ENDDO
         ENDDO

      ENDDO
  
      DEALLOCATE(CHI1)
      DEALLOCATE(CHISUB)
      DEALLOCATE(PSI1)
      DEALLOCATE(PSISUB)
      IF(ALLOCATED(CHI_OUT)) DEALLOCATE(CHI_OUT)
      IF(ALLOCATED(PSI_OUT)) DEALLOCATE(PSI_OUT)
!
!     
!     END OF ROUTINE.
      RETURN
      END
