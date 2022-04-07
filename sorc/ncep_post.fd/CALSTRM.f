!> @file
!> @brief Subroutine that computes geo streamfunction.
!>
!> This routine computes the geostrophic streamfunction,
!> PSI, from the passed geopotential height field, Z.  
!> The formule used it PSI = G*Z/F0, where G is the
!> gravitational acceleration constant and F0 is a 
!> constant Coriolis parameter. F0 is set to be the
!> valus of the Coriolis parameter near the center
!> of the model grid.
!>
!> @param[in] Z1D Geopotential height (m).
!> @param[out] STRM Geostrophic streamfunction.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1998-06-08 | T Black      | Conversion from 1-D TO 2-D
!> 2000-01-05 | Jim Tuccillo | MPI Version
!> 2002-06-13 | Mike Baldwin | WRF Version
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE CALSTRM(Z1D,STRM)

!     
!
!     
!     
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE OTHER PARAMETERS.
!     
!      use vrbls2d, only:
      use params_mod, only: g
      use ctlblk_mod, only: jsta, jend, im
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      real,PARAMETER :: OMEGA=7.292E-5,TWOMG=2*OMEGA
!
!     DECLARE VARIABLES.
!     
!      LOGICAL FIRST,OLDRD,RESTRT,RUN,SIGMA,STRD
      REAL, dimension(im,jsta:jend), intent(in)    ::  Z1D
      REAL, dimension(im,jsta:jend), intent(inout) ::  STRM
!
      LOGICAL OLDRD,STRD
      integer IMID,I,J
      real f0,gof0
!     
!***************************************************************************
!     START CALSTRM HERE.
!     
!     COMPUTE CORIOLIS PARAMETER AT 40N
!
      IMID=IM/2
      F0   = 1.454441e-4*sin(40.0*0.01745329)
      GOF0 = G/F0
!     
!     COMPUTE GEOSTROPHIC STREAMFUNCTION.
!$omp  parallel do
      DO J=JSTA,JEND
        DO I=1,IM
          STRM(I,J) = GOF0*Z1D(I,J)
        ENDDO
      ENDDO
!     
!     END OF ROUTINE.
      RETURN
      END
