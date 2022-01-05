!> @file
!> @brief Subroutine that computes drag cofficient.
!     
!> This rountine computes a surface layer drag coefficient using
!> equation (7.4.1A) in "An introduction to boundary layer
!> meteorology" by Stull (1988, Kluwer Academic Publishers).
!>     
!> @param[out] DRAGCO surface layer drag coefficient
!>
!> Program history
!> - 93-09-01  Russ Treadon
!> - 98-06-15  T Black - Conversion from 1-D to 2-D
!> - 00-01-04  Jim Tuccillo - MPI version           
!> - 02-01-15  Mike Baldwin - WRF version
!> - 05-02-22  H Chuang - Add WRF NMM components 
!>
!> @author Russ Treadon W/NP2 @date 1993-09-01
      SUBROUTINE CALDRG(DRAGCO)

!     
!
      use vrbls3d, only: uh, vh
      use vrbls2d, only: uz0, vz0, ustar, u10, v10
      use masks, only: lmh
      use params_mod, only: d00, d50, d25
      use ctlblk_mod, only: jsta, jend, jsta_m, jend_m, modelname, spval, im, jm,  &
                            jsta_2l, jend_2u, ista, iend, ista_m, iend_m, ista_2l, iend_2u
      use gridspec_mod, only: gridtype
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!     INCLUDE/SET PARAMETERS.
!     
!     DECLARE VARIABLES.
      REAL,intent(inout) ::  DRAGCO(ista_2l:iend_2u,jsta_2l:jend_2u)
      INTEGER IHE(JM),IHW(JM)
      integer I,J,LHMK,IE,IW,LMHK
      real UBAR,VBAR,WSPDSQ,USTRSQ,SUMU,SUMV,ULMH,VLMH,UZ0H,VZ0H
!     
!********************************************************************
!     START CALDRG HERE.
!     
!     INITIALIZE DRAG COEFFICIENT ARRAY TO ZERO.
!     
!$omp parallel do private(i,j)
      DO J=JSTA,JEND
        DO I=ISTA,IEND
!          DRAGCO(I,J) = D00
          DRAGCO(I,J) = 0.0 

        ENDDO
      ENDDO
!

      IF(gridtype=='A')THEN 
       DO J=JSTA,JEND
       DO I=ISTA,IEND
!     

       IF (USTAR(I,J) /= SPVAL) THEN

!        LMHK=NINT(LMH(I,J))
!     
!        COMPUTE A MEAN MASS POINT WIND SPEED BETWEEN THE
!        FIRST ATMOSPHERIC ETA LAYER AND Z0.  ACCORDING TO
!        NETCDF OUTPUT, UZ0 AND VZ0 ARE AT MASS POINTS. (MEB 6/11/02)
!
!        UBAR=D50*(UH(I,J,LMHK)+UZ0(I,J))
!        VBAR=D50*(VH(I,J,LMHK)+VZ0(I,J))
!        WSPDSQ=UBAR*UBAR+VBAR*VBAR

! dong use 10m wind
         WSPDSQ=U10(I,J)*U10(I,J)+V10(I,J)*V10(I,J)
!     
!        COMPUTE A DRAG COEFFICIENT.
!
        USTRSQ=USTAR(I,J)*USTAR(I,J)
        IF(WSPDSQ > 1.0) DRAGCO(I,J)=USTRSQ/WSPDSQ

       END IF
       ENDDO
       ENDDO
      ELSE IF(gridtype=='E')THEN
      
       DO J=JSTA_M,JEND_M
        IHE(J)=MOD(J+1,2)
        IHW(J)=IHE(J)-1
       ENDDO
       
       DO J=JSTA_M,JEND_M
       DO I=ISTA_M,IEND_M
!
!        COMPUTE A MEAN MASS POINT WIND IN THE
!        FIRST ATMOSPHERIC ETA LAYER.
!
        LMHK=NINT(LMH(I,J))
        IE=I+IHE(J)
        IW=I+IHW(J)
        SUMU=UH(IE,J,LMHK)+UH(IW,J,LMHK)+UH(I,J-1,LMHK)          &
          +UH(I,J+1,LMHK)
        SUMV=VH(IE,J,LMHK)+VH(IW,J,LMHK)+VH(I,J-1,LMHK)          &
          +VH(I,J+1,LMHK)
        ULMH=D25*SUMU
        VLMH=D25*SUMV
!
!        COMPUTE A MEAN MASS POINT WIND AT HEIGHT Z0.
!
        UZ0H=D25*(UZ0(IE,J)+UZ0(IW,J)+UZ0(I,J-1)+UZ0(I,J+1))
        VZ0H=D25*(VZ0(IE,J)+VZ0(IW,J)+VZ0(I,J-1)+VZ0(I,J+1))
!
!        COMPUTE A MEAN MASS POINT WIND SPEED BETWEEN THE
!        FIRST ATMOSPHERIC ETA LAYER AND Z0.
!
        UBAR=D50*(ULMH+UZ0H)
        VBAR=D50*(VLMH+VZ0H)
        WSPDSQ=UBAR*UBAR+VBAR*VBAR
!jjt  WSPDSQ=MIN(WSPDSQ,0.1)
!
!        COMPUTE A DRAG COEFFICIENT.
!
        USTRSQ=USTAR(I,J)*USTAR(I,J)
        IF(WSPDSQ > 1.0E-6)DRAGCO(I,J)=USTRSQ/WSPDSQ
!
       END DO
       END DO
      ELSE IF(gridtype=='B')THEN 
       DO J=JSTA_M,JEND_M
       DO I=ISTA_M,IEND_M
!
!        COMPUTE A MEAN MASS POINT WIND IN THE
!        FIRST ATMOSPHERIC ETA LAYER.
!
        LMHK=NINT(LMH(I,J))
        IE=I
        IW=I-1
        SUMU=UH(IE,J,LMHK)+UH(IW,J,LMHK)+UH(I,J-1,LMHK)          &
          +UH(IW,J-1,LMHK)
        SUMV=VH(IE,J,LMHK)+VH(IW,J,LMHK)+VH(I,J-1,LMHK)          &
          +VH(IW,J-1,LMHK)
        ULMH=D25*SUMU
        VLMH=D25*SUMV
!
!        COMPUTE A MEAN MASS POINT WIND AT HEIGHT Z0.
!
! NEMS-NMMB is now putting uz0 and vz0 on mass points to save double interpolation time in
! the model
        IF(MODELNAME == 'NMM')THEN
          UZ0H=UZ0(I,J)
          VZ0H=VZ0(I,J)
        ELSE   
          UZ0H=D25*(UZ0(IE,J)+UZ0(IW,J)+UZ0(I,J-1)+UZ0(IW,J-1))
          VZ0H=D25*(VZ0(IE,J)+VZ0(IW,J)+VZ0(I,J-1)+VZ0(IW,J-1))
        END IF  
!
!        COMPUTE A MEAN MASS POINT WIND SPEED BETWEEN THE
!        FIRST ATMOSPHERIC ETA LAYER AND Z0.
!
        UBAR=D50*(ULMH+UZ0H)
        VBAR=D50*(VLMH+VZ0H)
        WSPDSQ=UBAR*UBAR+VBAR*VBAR
!jjt  WSPDSQ=MIN(WSPDSQ,0.1)
!
!        COMPUTE A DRAG COEFFICIENT.
!
        USTRSQ=USTAR(I,J)*USTAR(I,J)
        IF(WSPDSQ > 1.0E-6)DRAGCO(I,J)=USTRSQ/WSPDSQ
!
       END DO
       END DO
      ELSE 
      
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            DRAGCO(I,J) = SPVAL
          ENDDO
        ENDDO

      END IF 
!     
!     END OF ROUTINE.
!     
      RETURN
      END
