!> @file
!> @brief Subroutine that computes Theta-E.
!>
!> This routine computes the equivalent potential temperature
!> given pressure, temperature, and specific humidity. The 
!> equations of Bolton (MWR,1980) are used.
!>
!> @param[in] P1D pressure (Pa).
!> @param[in] T1D temperature (K).
!> @param[in] Q1D specific humidity(kg/kg).
!> @param[out] THTE Theta-E (K).
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1993-06-18 | Russ Treadon | Initial
!> 1998-06-16 | T Black      | Convesion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI Version  
!> 2021-07-28 | W Meng       | Restrict computation from undefined grids
!> 2021-09-02 | Bo Cui       | Decompose UPP in X direction          
!>     
!> @author Russ Treadon W/NP2 @date 1993-06-18
!--------------------------------------------------------------------------------------
!> @brief Subroutine that computes Theta-E.
!> 
!> @param[in] P1D pressure (Pa).
!> @param[in] T1D temperature (K).
!> @param[in] Q1D specific humidity(kg/kg).
!> @param[out] THTE Theta-E (K).
!--------------------------------------------------------------------------------------
      SUBROUTINE CALTHTE(P1D,T1D,Q1D,THTE)

!
!     
      use params_mod, only: d00, eps, oneps, d01, h1m12, p1000, h1
      use ctlblk_mod, only: jsta, jend, im, spval, ista, iend
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      real,PARAMETER :: KG2G=1.E3
      real,PARAMETER :: D35=3.5,D4805=4.805,H2840=2840.,H55=55.
      real,PARAMETER :: D2845=0.2845,D00028=0.00028,D3376=3.376
      real,PARAMETER :: D00254=0.00254,D00081=0.00081,D81=0.81
      real,PARAMETER :: D28=0.28,H2675=2675.
!
!     DECLARE VARIABLES.
!     
      REAL,dimension(ista:iend,jsta:jend),intent(in)    :: P1D,T1D,Q1D
      REAL,dimension(ista:iend,jsta:jend),intent(inout) :: THTE

      integer I,J
      real P,T,Q,EVP,RMX,CKAPA,RKAPA,ARG,DENOM,TLCL,PLCL,FAC,   &
           ETERM,THETAE
!     
!***************************************************************
!     START CALTHTE.
!     
!     ZERO THETA-E ARRAY
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          THTE(I,J) = D00
        ENDDO
      ENDDO
!     
!     COMPUTE THETA-E.
!
!      DO J=JSTA_M,JEND_M
!      DO I=ISTA_M,IEND_M
!$omp parallel do private(i,j,p,t,q,evp,rmx,ckapa,rkapa,arg,denom,tlcl,plcl,fac,eterm,thetae)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
          IF(P1D(I,J)<spval.and.T1D(I,J)<spval.and.Q1D(I,J)<spval)THEN
          P        = P1D(I,J)
          T        = T1D(I,J)
          Q        = Q1D(I,J)
          EVP      = P*Q/(EPS+ONEPS*Q)
          RMX      = EPS*EVP/(P-EVP)
          CKAPA    = D2845*(1.-D28*RMX)
          RKAPA    = 1./CKAPA
          ARG      = max(H1M12, EVP*D01)
          DENOM    = D35*LOG(T) - LOG(EVP*D01) - D4805
          TLCL     = H2840/DENOM + H55
          PLCL     = P*(TLCL/T)**RKAPA
          FAC      = (P1000/P)**CKAPA
          ETERM    = (D3376/TLCL-D00254)*(RMX*KG2G*(H1+D81*RMX))
          THETAE   = T*FAC*EXP(ETERM)
          THTE(I,J)= THETAE
          ENDIF
        ENDDO
      ENDDO
!     
!     END OF ROUTINE.
!
      RETURN
      END
