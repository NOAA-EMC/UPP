!> @file
!> @brief Subroutine that computes T, Q, U, V on the flight levels (FD).
!>
!> This routine computes temperature, spec. hum, u wind component,
!> and v wind component on the NFD=6 FD levels. The 
!> height of these levels (in meters) is given in the 
!> data statement below. The alogrithm proceeds as
!> follows. (AGL-Above ground level in parentheses)
!> 
!> At each mass point move up vertically from the LM-TH (lowest
!> atmospheric) ETA layer. Find the ETA layers whose  
!> height (above ground) bounds the target FD level height.
!> Vertically interpolate to get temperature at this FD
!> level. Average the four surrounding winds 
!> to get a mass point wind. Vertically interpolate these 
!> mass point winds to the target FD level. Continue this 
!> process until all NFD=6 FD levels have been processed.
!> Move on to the next mass point.  
!>     
!> Averaging the four above ground winds to the mass point
!> was found to smooth the field and reduce the occurrence
!> of point peak winds far in excess of the winds at 
!> adjacent points. Mass point values are returned.
!>
!> @param[in] ITYPE Flag that determines whether MSL (1) or AGL (2) Levels are used.
!> @param[out] TFD Temperature (K) on FD levels.
!> @param[out] QFD Spec hum on FD levels.
!> @param[out] UFD U wind (m/s) on FD levels.
!> @param[out] VFD V wind (m/s) on FD levels.
!>     
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1993-11-23 | Russ Treadon | Corrected routine to compute FD levels with respect to mean sea level
!> 1994-01-04 | Mike Baldwin | Include options for computing either AGL or MSL
!> 1998-06-15 | T Black      | Conversion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI version           
!> 2002-01-15 | Mike Baldwin | WRF version           
!> 2011-12-14 | Sarah Lu     | Add GOCART aerosol AERFD
!> 2021-10-15 | JESSE MENG   | 2D DECOMPOSITION
!> 2022-09-22 | Li(Kate) Zhang   | Remove Dust=> AERFD
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE FDLVL(ITYPE,TFD,QFD,UFD,VFD,PFD,ICINGFD)

!     
!
      use vrbls3d,    only: ZMID, T, Q, PMID, ICING_GFIP, UH, VH
      use vrbls2d,    only: FIS
      use masks,      only: LMH
      use params_mod, only: GI, G
      use ctlblk_mod, only: JSTA, JEND, SPVAL, JSTA_2L, JEND_2U, LM, JSTA_M, &
                            JEND_M, HTFD, NFD, IM, JM, NBIN_DU, gocart_on,   &
                            MODELNAME, ISTA, IEND, ISTA_2L, IEND_2U, ISTA_M, IEND_M
      use gridspec_mod, only: GRIDTYPE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!     SET NUMBER OF FD LEVELS.
!jw      integer,intent(in) :: NFD ! coming from calling subroutine
!     
!     DECLARE VARIABLES
!     
      integer,intent(in) ::  ITYPE(NFD)
!jw      real,intent(in) :: HTFD(NFD)
      real,dimension(ISTA:IEND,JSTA:JEND,NFD),intent(out) :: TFD,QFD,UFD,VFD,PFD,ICINGFD
!
      INTEGER LVL(NFD),LHL(NFD)
      INTEGER IVE(JM),IVW(JM)
      REAL DZABV(NFD), DZABH(NFD)
      LOGICAL DONEH, DONEV
!jw
      integer I,J,JVS,JVN,IE,IW,JN,JS,JNT,L,LLMH,IFD,N
      integer ISTART,ISTOP,JSTART,JSTOP
      real htt,htsfc,httuv,dz,rdz,delt,delq,delu,delv,z1,z2,htabv,htabh,htsfcv
!
!     SET FD LEVEL HEIGHTS IN METERS.
!      DATA HTFD  / 30.E0,50.E0,80.E0,100.E0,305.E0,457.E0,610.E0,914.E0,1524.E0,  &
!          1829.E0,2134.E0,2743.E0,3658.E0,4572.E0,6000.E0/
!     
!****************************************************************
!     START FDLVL HERE
!     
!     INITIALIZE ARRAYS.
!     
!$omp  parallel do
      DO IFD = 1,NFD
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            TFD(I,J,IFD)     = SPVAL
            QFD(I,J,IFD)     = SPVAL
            UFD(I,J,IFD)     = SPVAL
            VFD(I,J,IFD)     = SPVAL
            PFD(I,J,IFD)     = SPVAL
            ICINGFD(I,J,IFD) = SPVAL
          ENDDO
        ENDDO
      ENDDO

      IF(gridtype == 'E') THEN
        JVN =  1
        JVS = -1
        do J=JSTA,JEND
          IVE(J) = MOD(J,2)
          IVW(J) = IVE(J)-1
        enddo
      END IF

      IF(gridtype /= 'A')THEN
        CALL EXCH(FIS(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U))
        DO L=1,LM
          CALL EXCH(ZMID(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,L))
        END DO
        ISTART = ISTA_M
        ISTOP  = IEND_M
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE
        ISTART = ISTA
        ISTOP  = IEND
        JSTART = JSTA
        JSTOP  = JEND
      END IF
      DO IFD = 1, NFD
!
!     MSL FD LEVELS
!
        IF (ITYPE(IFD)==1) THEN
!	write(6,*) 'computing above MSL'
!     
!       LOOP OVER HORIZONTAL GRID.
!	  
          DO J=JSTART,JSTOP
            DO I=ISTART,ISTOP
              HTSFC = FIS(I,J)*GI
              LLMH  = NINT(LMH(I,J))
!             IFD = 1
!     
!        LOCATE VERTICAL INDICES OF T,Q,U,V, LEVEL JUST 
!        ABOVE EACH FD LEVEL.
!
!        DO 22 IFD = 1, NFD
              DONEH=.FALSE.
              DONEV=.FALSE.
              DO L = LM,1,-1
                HTT = ZMID(I,J,L)
                IF(gridtype == 'E') THEN
                  IE = I+IVE(J)
                  IW = I+IVW(J)
                  JN = J+JVN
                  JS = J+JVS
                  HTTUV = 0.25*(ZMID(IW,J,L)                          &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
                ELSE IF(gridtype=='B')THEN
                  IE = I+1
                  IW = I
                  JN = J+1
                  JS = J
                  HTTUV = 0.25*(ZMID(IW,J,L)                          &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L)) 
                ELSE
                  HTTUV = HTT
                END IF
    
                IF (.NOT. DONEH .AND. HTT>HTFD(IFD)) THEN
                  LHL(IFD)   = L
                  DZABH(IFD) = HTT-HTFD(IFD)
                  DONEH = .TRUE.
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
                  IF(HTSFC > HTFD(IFD)) THEN
!mp
                     LHL(IFD) = LM+1  ! CHUANG: changed to lm+1
!mp
                  ENDIF
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
!               IFD        = IFD + 1
!               IF (IFD>NFD) GOTO 30
                END IF   
     
                IF (.NOT. DONEV .AND. HTTUV>HTFD(IFD)) THEN
                  LVL(IFD)   = L
                  DZABV(IFD) = HTTUV-HTFD(IFD)
                  DONEV=.TRUE.
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
                  IF(HTSFC>HTFD(IFD)) THEN
!mp
                    LVL(IFD)=LM+1  ! CHUANG: changed to lm+1
!mp
                  ENDIF
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
!               IFD        = IFD + 1
!               IF (IFD>NFD) GOTO 30
                ENDIF
    
                IF(DONEH .AND. DONEV) exit
              enddo           ! end of l loop
! 22     CONTINUE   	
!     
!        COMPUTE T, Q, U, AND V AT FD LEVELS.
!
!         DO 40 IFD = 1,NFD
 
              L = LHL(IFD)
              IF (L < LM) THEN
                DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                RDZ  = 1./DZ
                DELT = T(I,J,L)-T(I,J,L+1)
                DELQ = Q(I,J,L)-Q(I,J,L+1)
                TFD(I,J,IFD) = T(I,J,L) - DELT*RDZ*DZABH(IFD)
                QFD(I,J,IFD) = Q(I,J,L) - DELQ*RDZ*DZABH(IFD)
                PFD(I,J,IFD) = PMID(I,J,L) - (PMID(I,J,L)-PMID(I,J,L+1))*RDZ*DZABH(IFD)
                ICINGFD(I,J,IFD) = ICING_GFIP(I,J,L) - &
                 (ICING_GFIP(I,J,L)-ICING_GFIP(I,J,L+1))*RDZ*DZABH(IFD)
              ELSEIF (L == LM) THEN
                TFD(I,J,IFD) = T(I,J,L)
                QFD(I,J,IFD) = Q(I,J,L)
                PFD(I,J,IFD) = PMID(I,J,L)
                ICINGFD(I,J,IFD) = ICING_GFIP(I,J,L)
              ENDIF
    
              L = LVL(IFD)
              IF (L < LM) THEN
                IF(gridtype == 'E')THEN
                  IE = I+IVE(J)
                  IW = I+IVW(J)
                  JN = J+JVN
                  JS = J+JVS
                  Z1 = 0.25*(ZMID(IW,J,L)                              &
                     + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
                  Z2 = 0.25*(ZMID(IW,J,L+1)                            &
                     + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(I,JS,L+1))
                  DZ = Z1-Z2
       
                ELSE IF(gridtype=='B')THEN
                  IE =I+1
                  IW = I
                  JN = J+1
                  JS = J
                  Z1 = 0.25*(ZMID(IW,J,L)                              &
                     + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L))
                  Z2 = 0.25*(ZMID(IW,J,L+1)                            &
                     + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(IE,JN,L+1))
                  DZ = Z1-Z2
                ELSE 
                  DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                END IF 
                RDZ  = 1./DZ
                DELU = UH(I,J,L) - UH(I,J,L+1)
                DELV = VH(I,J,L) - VH(I,J,L+1)
                UFD(I,J,IFD) = UH(I,J,L) - DELU*RDZ*DZABV(IFD)
                VFD(I,J,IFD) = VH(I,J,L) - DELV*RDZ*DZABV(IFD)
              ELSEIF (L==LM) THEN
                UFD(I,J,IFD)=UH(I,J,L)
                VFD(I,J,IFD)=VH(I,J,L)
              ENDIF
! 40      CONTINUE
!     
!     COMPUTE FD LEVEL T, Q, U, AND V AT NEXT K.
!
            enddo        ! end of i loop
          enddo          ! end of j loop
!     END OF MSL FD LEVELS
        ELSE
!          write(6,*) 'computing above AGL'
!
!     AGL FD LEVELS 
!
!     
!     LOOP OVER HORIZONTAL GRID.
!     
          DO J=JSTART,JSTOP
            DO I=ISTART,ISTOP
              HTSFC = FIS(I,J)*GI
              IF(gridtype == 'E') THEN
                IE = I+IVE(J)
                IW = I+IVW(J)
                JN = J+JVN
                JS = J+JVS
                HTSFCV = (FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(I,JS))*(0.25/G)
              ELSE IF(gridtype == 'B')THEN
                IE = I+1
                IW = I
                JN = J+1
                JS = J
                HTSFCV = (FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(IE,JN))*(0.25/G)
              END IF
              LLMH  = NINT(LMH(I,J))
!             IFD   = 1
!     
!        LOCATE VERTICAL INDICES OF T,U,V, LEVEL JUST 
!        ABOVE EACH FD LEVEL.
!
!             DO 222 IFD = 1, NFD
              DONEH=.FALSE.
              DONEV=.FALSE.
              DO L = LLMH,1,-1
                HTABH = ZMID(I,J,L)-HTSFC
!                if(i==245.and.j==813)print*,'Debug FDL HTABH= ',htabh,zmid(i,j,l),htsfc
                IF(gridtype=='E')THEN
                  HTABV = 0.25*(ZMID(IW,J,L)                        &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))-HTSFCV
                ELSE IF(gridtype=='B')THEN
                  HTABV = 0.25*(ZMID(IW,J,L)                        &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L))-HTSFCV
                ELSE
                  HTABV = HTABH
                END IF
    
                IF (.NOT. DONEH .AND. HTABH>HTFD(IFD)) THEN
                  LHL(IFD)   = L
                  DZABH(IFD) = HTABH-HTFD(IFD)
                  DONEH=.TRUE.
!                 IFD        = IFD + 1
!                 IF (IFD>NFD) GOTO 230
                ENDIF

                IF (.NOT. DONEV .AND. HTABV>HTFD(IFD)) THEN
                  LVL(IFD)   = L
                  DZABV(IFD) = HTABV-HTFD(IFD)
                  DONEV = .TRUE.
!                 IFD        = IFD + 1
!                 IF (IFD>NFD) GOTO 230
                ENDIF
                IF(DONEH .AND. DONEV) exit
              enddo        ! end of l loop
!     
!        COMPUTE T, Q, U, AND V AT FD LEVELS.
!
! 222     CONTINUE
!
!             DO 240 IFD = 1,NFD
               L = LHL(IFD)
               IF (L<LM) THEN
                 DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                 RDZ  = 1./DZ
                 DELT = T(I,J,L)-T(I,J,L+1)
                 DELQ = Q(I,J,L)-Q(I,J,L+1)
                 TFD(I,J,IFD) = T(I,J,L) - DELT*RDZ*DZABH(IFD)
                 QFD(I,J,IFD) = Q(I,J,L) - DELQ*RDZ*DZABH(IFD)
                 PFD(I,J,IFD) = PMID(I,J,L) - (PMID(I,J,L)-PMID(I,J,L+1))*RDZ*DZABH(IFD)
                 ICINGFD(I,J,IFD) = ICING_GFIP(I,J,L) - &
                   (ICING_GFIP(I,J,L)-ICING_GFIP(I,J,L+1))*RDZ*DZABH(IFD)
               ELSE
                 TFD(I,J,IFD) = T(I,J,L)
                 QFD(I,J,IFD) = Q(I,J,L)
                 PFD(I,J,IFD) = PMID(I,J,L)
                 ICINGFD(I,J,IFD) = ICING_GFIP(I,J,L)
               ENDIF

               L = LVL(IFD)
               IF (L < LM) THEN
                 IF(gridtype == 'E')THEN
                   IE = I+IVE(J)
                   IW = I+IVW(J)
                   JN = J+JVN
                   JS = J+JVS
                   Z1 = 0.25*(ZMID(IW,J,L)                          &
                      + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
                   Z2 = 0.25*(ZMID(IW,J,L+1)                        &
                      + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(I,JS,L+1))
                   DZ = Z1-Z2
                 ELSE IF(gridtype=='B')THEN
                   IE = I+1
                   IW = I
                   JN = J+1
                   JS = J
                   Z1 = 0.25*(ZMID(IW,J,L)                          &
                      + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L))
                   Z2 = 0.25*(ZMID(IW,J,L+1)                        &
                      + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(IE,JN,L+1))
                   DZ = Z1-Z2
                 ELSE
                   DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                 END IF
                 RDZ  = 1./DZ
                 DELU = UH(I,J,L)-UH(I,J,L+1)
                 DELV = VH(I,J,L)-VH(I,J,L+1)
                 UFD(I,J,IFD) = UH(I,J,L) - DELU*RDZ*DZABV(IFD)
                 VFD(I,J,IFD) = VH(I,J,L) - DELV*RDZ*DZABV(IFD)
               ELSE
                 UFD(I,J,IFD) = UH(I,J,L)
                 VFD(I,J,IFD) = VH(I,J,L)
              ENDIF
! 240     CONTINUE
!     
!     COMPUTE FD LEVEL T, U, AND V AT NEXT K.
!
            enddo      ! end of i loop
          enddo        ! end of j loop
!     END OF AGL FD LEVELS
        ENDIF
      enddo          ! end of IFD loop

!  safety check to avoid tiny QFD values
     !KRF: Need NCAR and NMM WRF cores in this check as well?
     IF(MODELNAME=='RAPR' .OR. MODELNAME=='NCAR' .OR. MODELNAME=='NMM') THEN   !
       DO 420 IFD = 1,NFD
         DO J=JSTA,JEND
         DO I=ISTA,IEND
            if(QFD(I,J,IFD) < 1.0e-8) QFD(I,J,IFD)=0.0
         ENDDO
         ENDDO
420    CONTINUE
     endif
!
!     END OF ROUTINE.
!
      RETURN
      END

!> Computes FD level for u,v.
!>
!> This routine computes u/v wind component on NFD FD levels.
!> The height of these levels (in meters) is passed as an 
!> input parameter. The alogrithm proceeds as
!> follows. (AGL-Above ground level in parentheses)
!> 
!> At each mass point move up vertically from the LM-TH (lowest
!> atmospheric) ETA layer. Find the ETA layers whose  
!> height (above ground) bounds the target FD level height.
!> Vertically interpolate to get temperature at this FD
!> level. Average the four surrounding winds 
!> to get a mass point wind. Vertically interpolate these 
!> mass point winds to the target FD level. Continue this 
!> process until all NFD FD levels have been processed.
!> Move on to the next mass point.  
!>     
!> Averaging the four above ground winds to the mass point
!> was found to smooth the field and reduce the occurrence
!> of point peak winds far in excess of the winds at 
!> adjacent points. Mass point values are returned.
!>
!> @param[in] ITYPE Flag that determines whether MSL (1) or AGL (2) Levels are used.
!> @param[in] NFD Number of FD levels.
!> @param[in] HTFD FD levels.
!> @param[out] UFD U wind (m/s) on FD levels.
!> @param[out] VFD V wind (m/s) on FD levels.
!>     
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1993-11-23 | Russ Treadon | Corrected routine to compute FD levels with respect to mean sea level
!> 1994-01-04 | Mike Baldwin | Include options for computing either AGL or MSL
!> 1998-06-15 | T Black      | Conversion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI version           
!> 2002-01-15 | Mike Baldwin | WRF version           
!> 2011-12-14 | Sarah Lu     | Add GOCART aerosol AERFD
!> 2019-09-25 | Y Mao        | Seperate U/V from mass
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE FDLVL_UV(ITYPE,NFD,HTFD,UFD,VFD)
!     
!
      use vrbls3d,    only: ZMID, PMID, UH, VH
      use vrbls2d,    only: FIS
      use masks,      only: LMH
      use params_mod, only: GI, G
      use ctlblk_mod, only: JSTA, JEND, SPVAL, JSTA_2L, JEND_2U, LM, JSTA_M, &
                            JEND_M, IM, JM, MODELNAME, &
                            ISTA, IEND, ISTA_2L, IEND_2U, ISTA_M, IEND_M
      use gridspec_mod, only: GRIDTYPE
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!     DECLARE VARIABLES
!     
      integer,intent(in) ::  ITYPE(NFD)
      integer,intent(in) :: NFD ! coming from calling subroutine
      real,intent(in) :: HTFD(NFD)
      real,dimension(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,NFD),intent(out) :: UFD,VFD
!
      INTEGER LVL(NFD)
      INTEGER IVE(JM),IVW(JM)
      REAL  DZABV(NFD)
!jw
      integer I,J,JVS,JVN,IE,IW,JN,JS,L,LLMH,IFD,N
      integer ISTART,ISTOP,JSTART,JSTOP
      real htt,htsfc,httuv,dz,rdz,delu,delv,z1,z2,htabv,htabh,htsfcv
!
!****************************************************************
!     START FDLVL_UV HERE
!     
!     INITIALIZE ARRAYS.
!     
!$omp  parallel do
      DO IFD = 1,NFD
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            UFD(I,J,IFD)     = SPVAL
            VFD(I,J,IFD)     = SPVAL
          ENDDO
        ENDDO
      ENDDO

      IF(gridtype == 'E') THEN
        JVN =  1
        JVS = -1
        do J=JSTA,JEND
          IVE(J) = MOD(J,2)
          IVW(J) = IVE(J)-1
        enddo
      END IF

      IF(gridtype /= 'A')THEN
        CALL EXCH(FIS(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U))
        DO L=1,LM
          CALL EXCH(ZMID(ISTA_2L:IEND_2U,JSTA_2L:JEND_2U,L))
        END DO
        ISTART = ISTA_M
        ISTOP  = IEND_M
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE
        ISTART = ISTA
        ISTOP  = IEND
        JSTART = JSTA
        JSTOP  = JEND
      END IF
      DO IFD = 1, NFD
!
!     MSL FD LEVELS
!
        IF (ITYPE(IFD) == 1) THEN
!	write(6,*) 'computing above MSL'
!     
!       LOOP OVER HORIZONTAL GRID.
!	  
          DO J=JSTART,JSTOP
            DO I=ISTART,ISTOP
              HTSFC = FIS(I,J)*GI
              LLMH  = NINT(LMH(I,J))
!     
!        LOCATE VERTICAL INDICES OF U,V, LEVEL JUST 
!        ABOVE EACH FD LEVEL.
!
              DO L = LM,1,-1
                HTT = ZMID(I,J,L)
                IF(gridtype == 'E') THEN
                  IE = I+IVE(J)
                  IW = I+IVW(J)
                  JN = J+JVN
                  JS = J+JVS
                  HTTUV = 0.25*(ZMID(IW,J,L)                          &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
                ELSE IF(gridtype=='B')THEN
                  IE = I+1
                  IW = I
                  JN = J+1
                  JS = J
                  HTTUV = 0.25*(ZMID(IW,J,L)                          &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L)) 
                ELSE
                  HTTUV = HTT
                END IF

                IF (HTTUV > HTFD(IFD)) THEN
                  LVL(IFD)   = L
                  DZABV(IFD) = HTTUV-HTFD(IFD)
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
                  IF(HTSFC > HTFD(IFD)) THEN
!mp
                    LVL(IFD)=LM+1  ! CHUANG: changed to lm+1
!mp
                  ENDIF
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
    
                  exit
                ENDIF
              enddo           ! end of l loop
!     
!        COMPUTE U V AT FD LEVELS.
!
              L = LVL(IFD)
              IF (L < LM) THEN
                IF(gridtype == 'E')THEN
                  IE = I+IVE(J)
                  IW = I+IVW(J)
                  JN = J+JVN
                  JS = J+JVS
                  Z1 = 0.25*(ZMID(IW,J,L)                              &
                     + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
                  Z2 = 0.25*(ZMID(IW,J,L+1)                            &
                     + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(I,JS,L+1))
                  DZ = Z1-Z2
       
                ELSE IF(gridtype=='B')THEN
                  IE =I+1
                  IW = I
                  JN = J+1
                  JS = J
                  Z1 = 0.25*(ZMID(IW,J,L)                              &
                     + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L))
                  Z2 = 0.25*(ZMID(IW,J,L+1)                            &
                     + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(IE,JN,L+1))
                  DZ = Z1-Z2
                ELSE 
                  DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                END IF 
                RDZ  = 1./DZ
                DELU = UH(I,J,L) - UH(I,J,L+1)
                DELV = VH(I,J,L) - VH(I,J,L+1)
                UFD(I,J,IFD) = UH(I,J,L) - DELU*RDZ*DZABV(IFD)
                VFD(I,J,IFD) = VH(I,J,L) - DELV*RDZ*DZABV(IFD)
              ELSEIF (L == LM) THEN
                UFD(I,J,IFD)=UH(I,J,L)
                VFD(I,J,IFD)=VH(I,J,L)
              ELSE ! Underground
                UFD(I,J,IFD)=UH(I,J,LM)
                VFD(I,J,IFD)=VH(I,J,LM)
              ENDIF
!
            enddo        ! end of i loop
          enddo          ! end of j loop
!     END OF MSL FD LEVELS
        ELSE
!          write(6,*) 'computing above AGL'
!
!     AGL FD LEVELS 
!
!     
!     LOOP OVER HORIZONTAL GRID.
!     
          DO J=JSTART,JSTOP
            DO I=ISTART,ISTOP
              HTSFC = FIS(I,J)*GI
              IF(gridtype == 'E') THEN
                IE = I+IVE(J)
                IW = I+IVW(J)
                JN = J+JVN
                JS = J+JVS
                HTSFCV = (FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(I,JS))*(0.25/G)
              ELSE IF(gridtype == 'B')THEN
                IE = I+1
                IW = I
                JN = J+1
                JS = J
                HTSFCV = (FIS(IW,J)+FIS(IE,J)+FIS(I,JN)+FIS(IE,JN))*(0.25/G)
              END IF
              LLMH  = NINT(LMH(I,J))
!     
!        LOCATE VERTICAL INDICES OF U,V, LEVEL JUST 
!        ABOVE EACH FD LEVEL.
!
              DO L = LLMH,1,-1
                HTABH = ZMID(I,J,L)-HTSFC
                IF(gridtype=='E')THEN
                  HTABV = 0.25*(ZMID(IW,J,L)                        &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))-HTSFCV
                ELSE IF(gridtype=='B')THEN
                  HTABV = 0.25*(ZMID(IW,J,L)                        &
                        + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L))-HTSFCV
                ELSE
                  HTABV = HTABH
                END IF

                IF (HTABV > HTFD(IFD)) THEN
                  LVL(IFD)   = L
                  DZABV(IFD) = HTABV-HTFD(IFD)
!                 IFD        = IFD + 1
                  exit
                ENDIF
              enddo        ! end of l loop
!     
!        COMPUTE U V AT FD LEVELS.
!
               L = LVL(IFD)
               IF (L < LM) THEN
                 IF(gridtype == 'E')THEN
                   IE = I+IVE(J)
                   IW = I+IVW(J)
                   JN = J+JVN
                   JS = J+JVS
                   Z1 = 0.25*(ZMID(IW,J,L)                          &
                      + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(I,JS,L))
                   Z2 = 0.25*(ZMID(IW,J,L+1)                        &
                      + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(I,JS,L+1))
                   DZ = Z1-Z2
                 ELSE IF(gridtype=='B')THEN
                   IE = I+1
                   IW = I
                   JN = J+1
                   JS = J
                   Z1 = 0.25*(ZMID(IW,J,L)                          &
                      + ZMID(IE,J,L)+ZMID(I,JN,L)+ZMID(IE,JN,L))
                   Z2 = 0.25*(ZMID(IW,J,L+1)                        &
                      + ZMID(IE,J,L+1)+ZMID(I,JN,L+1)+ZMID(IE,JN,L+1))
                   DZ = Z1-Z2
                 ELSE
                   DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                 END IF
                 RDZ  = 1./DZ
                 DELU = UH(I,J,L)-UH(I,J,L+1)
                 DELV = VH(I,J,L)-VH(I,J,L+1)
                 UFD(I,J,IFD) = UH(I,J,L) - DELU*RDZ*DZABV(IFD)
                 VFD(I,J,IFD) = VH(I,J,L) - DELV*RDZ*DZABV(IFD)
               ELSE
                 UFD(I,J,IFD) = UH(I,J,L)
                 VFD(I,J,IFD) = VH(I,J,L)
              ENDIF
!     
!     COMPUTE FD LEVEL T, U, AND V AT NEXT K.
!
            enddo      ! end of i loop
          enddo        ! end of j loop
!     END OF AGL FD LEVELS
        ENDIF
      enddo          ! end of IFD loop

      RETURN
      END

!> Computes FD level for mass variables.
!>
!> This routine computes mass variables (temperature, spec. hum...)
!> on NFD FD levels. The height of these levels (in meters) is 
!> passed as an input parameter. The alogrithm proceeds as
!> follows. (AGL-Above ground level in parentheses)
!> 
!> At each mass point move up vertically from the LM-TH (lowest
!> atmospheric) ETA layer. Find the ETA layers whose  
!> height (above ground) bounds the target FD level height.
!> Vertically interpolate to get temperature at this FD
!> level. Average the four surrounding winds 
!> to get a mass point wind. Vertically interpolate these 
!> mass point winds to the target FD level. Continue this 
!> process until all NFD FD levels have been processed.
!> Move on to the next mass point.
!>     
!> Averaging the four above ground winds to the mass point
!> was found to smooth the field and reduce the occurrence
!> of point peak winds far in excess of the winds at 
!> adjacent points. Mass point values are returned.
!>
!> NOTES for Q fields by Y Mao:
!> The following safety check should be executed by the caller of FDLVL subroutines.
!> Safety check to avoid tiny QFD values.
!> KRF: Need NCAR and NMM WRF cores in this check as well?
!> @code
!>    IF(MODELNAME=='RAPR' .OR. MODELNAME=='NCAR' .OR. MODELNAME=='NMM') THEN   !
!>      DO IFD = 1,NFD
!>        DO J=JSTA,JEND
!>        DO I=1,IM
!>           if(QFD(I,J,IFD) < 1.0e-8) QFD(I,J,IFD)=0.0
!>        ENDDO
!>        ENDDO
!>      ENDDO
!>    endif
!> @endcode 
!>
!> @param[in] ITYPE Flag that determines whether MSL (1) or AGL (2) Levels are used.
!> @param[in] NFD Number of FD levels.
!> @param[in] PTFD FD pressure levels.
!> @param[in] HTFD FD height levels.
!> @param[in] NIN Number of input fields.
!> @param[in] QIN Array of mass point value on model levels.
!> @param[in] QTYPE Charater array of variable type to differentiate underground interpolation.
!><pre>
!>                   C-5 Cloud Species
!>                   K-TURBULENT KINETIC ENERGY
!>                   Q-Specific Humidity
!>                   T-Temperature, 
!>                   W-Vertical Velocity or Omega
!></pre>
!> @param[out] QFD Array of mass point value on FD levels.
!>     
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-12-22 | Russ Treadon | Initial
!> 1993-11-23 | Russ Treadon | Corrected routine to compute FD levels with respect to mean sea level
!> 1994-01-04 | Mike Baldwin | Include options for computing either AGL or MSL
!> 1998-06-15 | T Black      | Conversion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI version           
!> 2002-01-15 | Mike Baldwin | WRF version           
!> 2011-12-14 | Sarah Lu     | Add GOCART aerosol AERFD
!> 2017-06-01 | Y Mao        | Add FD levels for GTG(EDPARM CATEDR MWTURB) and allow levels input from control file
!> 2019-09-25 | Y Mao        | Seperate mass from UV allow array of mass input to interpolate multiple fields with the same levels at one time. Dust=> AERFD can be processed when NIN=NBIN_DU
!> 2020-11-10 | Jesse Meng   | Use UPP_PHYSICS module
!>
!> @author Russ Treadon W/NP2 @date 1992-12-22
      SUBROUTINE FDLVL_MASS(ITYPE,NFD,PTFD,HTFD,NIN,QIN,QTYPE,QFD)
      use vrbls3d,    only: T,Q,ZMID,PMID,PINT,ZINT
      use vrbls2d,    only: FIS
      use masks,      only: LMH
      use params_mod, only: GI, G, GAMMA,PQ0, A2, A3, A4, RHMIN,RGAMOG
      use ctlblk_mod, only: JSTA, JEND, SPVAL, JSTA_2L, JEND_2U, LM, JSTA_M, &
                            JEND_M, IM, JM,global,MODELNAME, &
                            ISTA, IEND, ISTA_2L, IEND_2U, ISTA_M, IEND_M
      use gridspec_mod, only: GRIDTYPE
      use physcons_post,only: CON_FVIRT, CON_ROG, CON_EPS, CON_EPSM1
      use upp_physics,  only: FPVSNEW
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!     SET NUMBER OF FD LEVELS.
!     
!     DECLARE VARIABLES
!     
      real,parameter:: zshul=75.,tvshul=290.66

      integer,intent(in) ::  ITYPE(NFD)
      integer,intent(in) :: NFD ! coming from calling subroutine
      real, intent(in) :: PTFD(NFD)
      real,intent(in) :: HTFD(NFD)
      integer,intent(in) :: NIN
      real,intent(in) :: QIN(ISTA:IEND,JSTA:JEND,LM,NIN)
      character, intent(in) :: QTYPE(NIN)
      real,intent(out) :: QFD(ISTA:IEND,JSTA:JEND,NFD,NIN)

!
      INTEGER LHL(NFD)
      REAL DZABH(NFD)
!jw
      integer I,J,L,LLMH,IFD,N
      integer ISTART,ISTOP,JSTART,JSTOP
      real htt,htsfc,dz,rdz,delq,htabh

      real :: tvu,tvd,gammas,part,ES,QSAT,RHL,PL,ZL,TL,QL
      real :: TVRL,TVRBLO,TBLO,QBLO
!
!****************************************************************
!     START FDLVL_MASS HERE
!     
!     INITIALIZE ARRAYS.
!     
!$omp  parallel do
      DO N=1,NIN
      DO IFD = 1,NFD
        DO J=JSTA,JEND
          DO I=ISTA,IEND
            QFD(I,J,IFD,N)     = SPVAL
          ENDDO
        ENDDO
      ENDDO
      ENDDO

      IF(gridtype /= 'A')THEN
        ISTART = ISTA_M
        ISTOP  = IEND_M
        JSTART = JSTA_M
        JSTOP  = JEND_M
      ELSE
        ISTART = ISTA
        ISTOP  = IEND
        JSTART = JSTA
        JSTOP  = JEND
      END IF

      DO IFD = 1, NFD

!
!     MSL FD LEVELS
!
        IF (ITYPE(IFD) == 1) THEN
!	write(6,*) 'computing above MSL'
!     
!       LOOP OVER HORIZONTAL GRID.
!	  
          DO J=JSTART,JSTOP
            DO I=ISTART,ISTOP
              HTSFC = FIS(I,J)*GI
              LLMH  = NINT(LMH(I,J))
!     
!        LOCATE VERTICAL INDICES OF Q, LEVEL JUST 
!        ABOVE EACH FD LEVEL.
!
              DO L = LM,1,-1
                HTT = ZMID(I,J,L)
    
                IF (HTT > HTFD(IFD)) THEN
                  LHL(IFD)   = L
                  DZABH(IFD) = HTT-HTFD(IFD)
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL
                  IF(HTSFC > HTFD(IFD)) THEN
!mp
                     LHL(IFD) = LM+1  ! CHUANG: changed to lm+1
!mp
                  ENDIF
! THIS SHOULD SET BELOW GROUND VALUES TO SPVAL

                  exit
                END IF   

              ENDDO           ! end of L loop
!     
!        COMPUTE Q AT FD LEVELS.
!
              L = LHL(IFD)
              IF (L < LM) THEN
                DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                RDZ  = 1./DZ
                DO N = 1, NIN
                   if(QIN(I,J,L,N)<SPVAL) then
                      QFD(I,J,IFD,N)=QIN(I,J,L+1,N)
                   elseif(QIN(I,J,L+1,N)<SPVAL) then
                      QFD(I,J,IFD,N)=QIN(I,J,L,N)
                   else
                      QFD(I,J,IFD,N) = QIN(I,J,L,N) - &
                         (QIN(I,J,L,N)-QIN(I,J,L+1,N))*RDZ*DZABH(IFD)
                   endif
                ENDDO
              ELSEIF (L == LM) THEN
                DO N = 1, NIN
                   QFD(I,J,IFD,N) = QIN(I,J,L,N)
                ENDDO
              ELSE ! Underground
                 DO N = 1, NIN
                   ! Deduce T and Q differently by different models
                   IF(MODELNAME == 'GFS')THEN ! GFS deduce T using Shuell
                      if(QTYPE(N) == "T" .or. QTYPE(N) == "Q") then
                         tvu = T(I,J,LM) * (1.+con_fvirt*Q(I,J,LM))
                         if(ZMID(I,J,LM) > zshul) then
                            tvd = tvu + gamma*ZMID(I,J,LM)
                            if(tvd > tvshul) then
                               if(tvu > tvshul) then
                                  tvd = tvshul - 5.e-3*(tvu-tvshul)*(tvu-tvshul)
                               else
                                  tvd = tvshul
                               endif
                            endif
                            gammas = (tvu-tvd)/ZMID(I,J,LM)
                         else
                            gammas = 0.
                         endif
                         part = con_rog*(LOG(PTFD(IFD))-LOG(PMID(I,J,LM)))
                         part = ZMID(I,J,LM) - tvu*part/(1.+0.5*gammas*part)
                         part = T(I,J,LM) - gamma*(part-ZMID(I,J,LM))

                         if(QTYPE(N) == "T") QFD(I,J,IFD,N) = part

                         if(QTYPE(N) == "Q") then

! Compute RH at lowest model layer because Iredell and Chuang decided to compute
! underground GFS Q to maintain RH
                            ES   = min(FPVSNEW(T(I,J,LM)), PMID(I,J,LM))
                            QSAT = CON_EPS*ES/(PMID(I,J,LM)+CON_EPSM1*ES)
                            RHL  = Q(I,J,LM)/QSAT
! compute saturation water vapor at isobaric level
                            ES   = min(FPVSNEW(part), PTFD(IFD))
                            QSAT = CON_EPS*ES/(PTFD(IFD)+CON_EPSM1*ES)
!     Q at isobaric level is computed by maintaining constant RH  
                            QFD(I,J,IFD,N) = RHL*QSAT
                         endif
                      endif

                   ELSE
                      if(QTYPE(N) == "T" .or. QTYPE(N) == "Q") then
                         PL = PINT(I,J,LM-1)
                         ZL = ZINT(I,J,LM-1)
                         TL = 0.5*(T(I,J,LM-2)+T(I,J,LM-1))
                         QL = 0.5*(Q(I,J,LM-2)+Q(I,J,LM-1))

                         QSAT = PQ0/PL*EXP(A2*(TL-A3)/(TL-A4))
                         RHL  = QL/QSAT
!
                         IF(RHL > 1.)THEN
                            RHL = 1.
                            QL  = RHL*QSAT
                         ENDIF
!
                         IF(RHL < RHmin)THEN
                            RHL = RHmin
                            QL  = RHL*QSAT
                         ENDIF
!
                         TVRL   = TL*(1.+0.608*QL)
                         TVRBLO = TVRL*(PTFD(IFD)/PL)**RGAMOG
                         TBLO   = TVRBLO/(1.+0.608*QL)

                         QSAT     = PQ0/PTFD(IFD)*EXP(A2*(TBLO-A3)/(TBLO-A4))
                         if(QTYPE(N) == "T") QFD(I,J,IFD,N) = TBLO
                         QBLO     = RHL*QSAT
                         if(QTYPE(N) == "Q") QFD(I,J,IFD,N) = MAX(1.E-12,QBLO)
                      endif
                   END IF       ! endif loop for deducing T and Q differently for GFS  

                   if(QTYPE(N) == "W") QFD(I,J,IFD,N)=QIN(I,J,LM,N) ! W OMGA
                   if(QTYPE(N) == "K") QFD(I,J,IFD,N)= max(0.0,0.5*(QIN(I,J,LM,N)+QIN(I,J,LM-1,N))) ! TKE
                   if(QTYPE(N) == "C") QFD(I,J,IFD,N)=0.0 ! Hydrometeor fields
                 END DO

              ENDIF ! Underground
    
!     
!     COMPUTE FD LEVEL Q AT NEXT K.
!
            enddo        ! end of i loop
          enddo          ! end of j loop
!     END OF MSL FD LEVELS
        ELSE
!          write(6,*) 'computing above AGL'
!
!     AGL FD LEVELS 
!
!     
!     LOOP OVER HORIZONTAL GRID.
!     
          DO J=JSTART,JSTOP
            DO I=ISTART,ISTOP
              HTSFC = FIS(I,J)*GI
              LLMH  = NINT(LMH(I,J))
!     
!        LOCATE VERTICAL INDICES OF Q, LEVEL JUST 
!        ABOVE EACH FD LEVEL.
!
              DO L = LLMH,1,-1
                HTABH = ZMID(I,J,L)-HTSFC
    
                IF ( HTABH > HTFD(IFD)) THEN
                  LHL(IFD)   = L
                  DZABH(IFD) = HTABH-HTFD(IFD)

                  exit
                ENDIF
              enddo        ! end of l loop
!     
!        COMPUTE Q AT FD LEVELS.
!
               L = LHL(IFD)
               IF (L < LM) THEN
                 DZ   = ZMID(I,J,L)-ZMID(I,J,L+1)
                 RDZ  = 1./DZ
                 DO N = 1, NIN
                    if(QIN(I,J,L,N)<SPVAL) then
                       QFD(I,J,IFD,N)=QIN(I,J,L+1,N)
                    elseif(QIN(I,J,L+1,N)<SPVAL) then
                       QFD(I,J,IFD,N)=QIN(I,J,L,N)
                    else
                       QFD(I,J,IFD,N) = QIN(I,J,L,N) - &
                       (QIN(I,J,L,N)-QIN(I,J,L+1,N))*RDZ*DZABH(IFD)
                    endif
                 ENDDO
               ELSE
                  DO N = 1, NIN
                     QFD(I,J,IFD,N) = QIN(I,J,L,N)
                  ENDDO
               ENDIF

!     
!     COMPUTE FD LEVEL Q AT NEXT K.
!
            enddo      ! end of i loop
          enddo        ! end of j loop
!     END OF AGL FD LEVELS
        ENDIF
      enddo          ! end of IFD loop

!
!     END OF ROUTINE.
!
      RETURN
      END
