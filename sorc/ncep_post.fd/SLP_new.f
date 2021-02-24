      SUBROUTINE MEMSLP(TPRES,QPRES,FIPRES)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .
!   SUBROUTINE:  MEMSLP      MEMBRANE SLP REDUCTION
!
! ABSTRACT:  THIS ROUTINE COMPUTES THE SEA LEVEL PRESSURE
!            REDUCTION USING THE MESINGER RELAXATION
!            METHOD FOR SIGMA COORDINATES.
!            A BY-PRODUCT IS THE
!            SET OF VALUES FOR THE UNDERGROUND TEMPERATURES
!            ON THE SPECIFIED PRESSURE LEVELS
!
! PROGRAM HISTORY LOG:
!   99-09-23  T BLACK - REWRITTEN FROM ROUTINE SLP (ETA
!                       COORDINATES)
!   02-07-26  H CHUANG - PARALLIZE AND MODIFIED FOR WRF A/C GRIDS
!                        ALSO REDUCE S.O.R. COEFF FROM 1.75 to 1.25 
!                        BECAUSE THERE WAS NUMERICAL INSTABILITY  
!   02-08-21  H CHUANG - MODIFIED TO ALWAYS USE OLD TTV FOR RELAXATION
!                        SO THAT THERE WAS BIT REPRODUCIBILITY BETWEEN
!                        USING ONE AND MULTIPLE TASKS	      
!   11-04-29  H CHUANG - FIX GFS GIBSING BY USING LM-1 STATE VARIABLES
!                        TO DERIVE SLP HYDROSTATICALLY
!   13-12-06  H CHUANG - REMOVE EXTRA SMOOTHING OF SLP ITSELF  
!                        CHANGES TO AVOID RELAXATION FOR ABOVE G GIBSING
!                        ARE COMMENTED OUT FOR NOW
!   19-10-30  Bo CUI - REMOVE "GOTO" STATEMENT
!
! USAGE:  CALL SLPSIG FROM SUBROUITNE ETA2P
!
!   INPUT ARGUMENT LIST:
!     PD   - SFC PRESSURE MINUS PTOP
!     FIS  - SURFACE GEOPOTENTIAL
!     T    - TEMPERATURE 
!     Q    - SPECIFIC HUMIDITY
!     FI   - GEOPOTENTIAL
!     PT   - TOP PRESSURE OF DOMAIN
!
!   OUTPUT ARGUMENT LIST:
!     PSLP - THE FINAL REDUCED SEA LEVEL PRESSURE ARRAY
!
!   SUBPROGRAMS CALLED:
!     UNIQUE:
!             NONE
!
!-----------------------------------------------------------------------
      use vrbls3d,    only: pint, zint, t, q
      use vrbls2d,    only: pslp, fis
      use masks,      only: lmh
      use params_mod, only: overrc, ad05, cft0, g, rd, d608, h1, kslpd
      use ctlblk_mod, only: jend, jsta, spval, spl, num_procs, mpi_comm_comp, lsmp1, &
                            jsta_m, jend_m, lm, im, jsta_2l, jend_2u, lsm, jm,&
                            im_jm
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!      
      INCLUDE "mpif.h"
!-----------------------------------------------------------------------
      integer,PARAMETER   :: NFILL=0,NRLX1=500,NRLX2=100
      real,parameter:: def_of_mountain=2.0
!-----------------------------------------------------------------------
      real,dimension(IM,JSTA_2L:JEND_2U,LSM),intent(in) :: QPRES
      real,dimension(IM,JSTA_2L:JEND_2U,LSM),intent(inout) :: TPRES,FIPRES
      REAL  ::  TTV(IM,JSTA_2L:JEND_2U),TNEW(IM,JSTA_2L:JEND_2U)        &
        ,      P1(IM,JSTA_2L:JEND_2U),HTM2D(IM,JSTA_2L:JEND_2U)
      REAL  :: HTMO(IM,JSTA_2L:JEND_2U,LSM)    
      real  :: P2,TLYR,GZ1,GZ2,SPLL,PSFC,PCHK,SLOPE,TVRTC,DIS,TVRT,tem
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
      INTEGER :: KMNTM(LSM),IMNT(IM_JM,LSM),JMNT(IM_JM,LSM)             &
         ,       LMHO(IM,JSTA_2L:JEND_2U)
      INTEGER :: IHE(JM),IHW(JM),IVE(JM),IVW(JM),IHS(JM),IHN(JM)
      integer    ii,jj,I,J,L,N,LLMH,KM,KS,IHH2,KOUNT,KMN,NRLX,LHMNT,    &
                 LMHIJ,LMAP1,KMM,LP,LXXX,IERR
! dong
      real a1,a2,a3,a4,a5,a6,a7,a8
!-----------------------------------------------------------------------
      LOGICAL :: DONE(IM,JSTA_2L:JEND_2U)
!-----------------------------------------------------------------------
!***
!***  CALCULATE THE I-INDEX EAST-WEST INCREMENTS
!***
!
      ii = IM/2
      jj = (JEND-JSTA)/2
      DO J=1,JM
        IHE(J) =  1
        IHW(J) = -1
        IHS(J) = -1
        IHN(J) =  1
        IVE(J) = MOD(J,2)
        IVW(J) = IVE(J)-1
      ENDDO
!      print*,'relaxation coeff= ',OVERRC
!-----------------------------------------------------------------------
!***
!***  INITIALIZE ARRAYS.  LOAD SLP ARRAY WITH SURFACE PRESSURE.
!***
!$omp parallel do  private(i,j,llmh)
      DO J=JSTA,JEND
        DO I=1,IM
          LLMH      = NINT(LMH(I,J))
          PSLP(I,J) = PINT(I,J,LLMH+1)
! dong
!          TTV(I,J)  = 0.
          TTV(I,J)  = spval 
          TNEW(I,J)  = spval


          LMHO(I,J) = 0
          DONE(I,J) = .FALSE.
        ENDDO
      ENDDO
!
!--------------------------------------------------------------------
!***
!***  CREATE A 3-D "HEIGHT MASK" FOR THE SPECIFIED PRESSURE LEVELS
!***  (1 => ABOVE GROUND) AND A 2-D INDICATOR ARRAY THAT SAYS 
!***  WHICH PRESSURE LEVEL IS THE LOWEST ONE ABOVE THE GROUND
!***
      DO L=1,LSM
        SPLL = SPL(L)
!       
!$omp parallel do private(j,i,psfc,pchk)
        DO J=JSTA,JEND
          DO I=1,IM

           if(PSLP(I,J)<spval) then

            PSFC = PSLP(I,J)
            PCHK = PSFC
            IF(NFILL > 0) THEN
              PCHK = PINT(I,J,NINT(LMH(I,J))+1-NFILL)
            ENDIF
            IF(FIS(I,J) < 1.) PCHK = PSFC
!
            IF(SPLL < PCHK) THEN
              HTMO(I,J,L) = 1.
            ELSE
              HTMO(I,J,L) = 0.
              IF(L > 1 .AND. HTMO(I,J,L-1) > 0.5) LMHO(I,J) = L-1
            ENDIF
            IF(L == LSM .AND. HTMO(I,J,L) > 0.5) LMHO(I,J) = LSM
!
! test new idea of filtering above-ground pressure levels for Gibsing
!           IF(L==LSM.AND.HTMO(I,J,L)>0.5)THEN
!	      IF(FIS(I,J)>0.)THEN 
!	        LMHO(I,J)=LSM
!	      ELSE
!	        LMHO(I,J)=LSM-2
!	        HTMO(I,J,LSM)=0.
!	        HTMO(I,J,LSM-1)=0. 
!	      END IF
!	    END IF  
!           if(i==ii.and.j==jj)print*,'Debug: HTMO= ',HTMO(I,J,L)

           endif  !if pslp

          ENDDO
        ENDDO
!
      ENDDO
!      if(jj>=jsta.and.jj<=jend) print*,'Debug: LMHO=',LMHO(ii,jj)
!--------------------------------------------------------------------
!***
!***  WE REACH THIS LINE IF WE WANT THE MESINGER ETA SLP REDUCTION
!***  BASED ON RELAXATION TEMPERATURES.  THE FIRST STEP IS TO
!***  FIND THE HIGHEST LAYER CONTAINING MOUNTAINS.
!***
      LHMNT = LSM 
      LOOP210: DO L=LSM,1,-1

        DO J=JSTA,JEND
          DO I=1,IM
           if(PSLP(I,J)<spval) then
            IF(HTMO(I,J,L) < 0.5) CYCLE LOOP210
           endif
          ENDDO
        ENDDO
        LHMNT = L+1
        EXIT LOOP210
        ENDDO LOOP210
 !210  continue
 !220  continue

!      print*,'Debug in SLP: LHMNT=',LHMNT

      if ( num_procs > 1 ) then
        CALL MPI_ALLREDUCE                                      &
          (LHMNT,LXXX,1,MPI_INTEGER,MPI_MIN,MPI_COMM_COMP,IERR)
        LHMNT = LXXX
      end if
   
      IF(LHMNT == LSMP1) GO TO 325

!      print*,'Debug in SLP: LHMNT A ALLREDUCE=',LHMNT
!***
!***  NOW GATHER THE ADDRESSES OF ALL THE UNDERGROUND POINTS.
!***
!!$omp parallel do private(kmn,kount)
      DO L=LHMNT,LSM
        KMN      = 0
        KMNTM(L) = 0
        KOUNT    = 0
!       DO 240 J=JSTA_M2,JEND_M2
        DO J=JSTA_M,JEND_M
          DO I=2,IM-1
           if(PSLP(I,J)<spval) then
            KOUNT = KOUNT + 1
            IMNT(KOUNT,L) = 0
            JMNT(KOUNT,L) = 0
            IF(HTMO(I,J,L) > 0.5) cycle
            KMN         = KMN + 1
            IMNT(KMN,L) = I
            JMNT(KMN,L) = J
           endif
          enddo
        enddo
        KMNTM(L) = KMN
      enddo
!
!
!***  CREATE A TEMPORARY TV ARRAY, AND FOLLOW BY SEQUENTIAL
!***  OVERRELAXATION, DOING NRLX PASSES.
!
!     IF(NTSD==1)THEN
        NRLX = NRLX2
!     ELSE
!       NRLX=NRLX2
!     ENDIF
!
!!$omp parallel do private(i,j,ttv,tem,kmma (Can this loop be threaded?))
      DO L=LHMNT,LSM
!
!$omp parallel do private(i,j)
        DO J=JSTA,JEND
          DO I=1,IM
! dong
!            if (QPRES(I,J,LSM) < spval) then 
           if(PSLP(I,J)<spval) then
            TTV(I,J)   = TPRES(I,J,L)
            HTM2D(I,J) = HTMO(I,J,L)
           end if ! spval if
!            end if ! spval if
!     IF(TTV(I,J)<150. .and. TTV(I,J)>325.0)print*                &  
!       ,'abnormal IC for T relaxation',i,j,TTV(I,J)
          enddo
        enddo
!
!***  FOR GRID BOXES NEXT TO MOUNTAINS, COMPUTE TV TO USE AS
!***  BOUNDARY CONDITIONS FOR THE RELAXATION UNDERGROUND
!
        CALL EXCH(HTM2D(1,JSTA_2L))  !ONLY NEED TO EXCHANGE ONE ROW FOR A/C GRID
!       DO J=JSTA_M2,JEND_M2
!$omp parallel do private(i,j,tem)
        DO J=JSTA_M,JEND_M
          DO I=2,IM-1

           if(PSLP(I,J)<spval) then

!HC        IF(HTM2D(I,J,L)>0.5.AND.
!HC     1     HTM2D(I+IHW(J),J-1,L)*HTM2D(I+IHE(J),J-1,L)
!HC     2    *HTM2D(I+IHW(J),J+1,L)*HTM2D(I+IHE(J),J+1,L)
!HC     3    *HTM2D(I-1     ,J  ,L)*HTM2D(I+1     ,J  ,L)
!HC     4    *HTM2D(I       ,J-2,L)*HTM2D(I       ,J+2,L)<0.5)THEN
!HC MODIFICATION FOR C AND A GRIDS

            tem = HTM2D(I-1,J)*HTM2D(I+1,J)*HTM2D(I,J-1)*HTM2D(I,J+1)         &
                * HTM2D(I-1,J-1)*HTM2D(I+1,J-1)*HTM2D(I-1,J+1)*HTM2D(I+1,J+1)
            IF(HTM2D(I,J) > 0.5 .AND. tem < 0.5) then
              TTV(I,J) = TPRES(I,J,L)*(1.+0.608*QPRES(I,J,L))
            ENDIF
!           if(i==ii.and.j==jj)print*,'Debug:L,TTV B SMOO= ',l,TTV(I,J) 
           end if ! spval
          ENDDO
        ENDDO
!
        KMM = KMNTM(L)

!      print*,'Debug:L,KMM=',L,KMM
!
        DO N=1,NRLX
          CALL EXCH(TTV(1,JSTA_2L))
!$omp parallel do private(i,j,km)
          DO KM=1,KMM
            I = IMNT(KM,L)
            J = JMNT(KM,L)

            if(PSLP(I,J)<spval) then

!HC      TTV(I,J)=AD05*(4.*(TTV(I+IHW(J),J-1)+TTV(I+IHE(J),J-1)
!HC     1                  +TTV(I+IHW(J),J+1)+TTV(I+IHE(J),J+1))
!HC     2                  +TTV(I-1,J)       +TTV(I+1,J)
!HC     3                  +TTV(I,J-2)       +TTV(I,J+2))
!HC     4                  -CFT0*TTV(I,J)
!HC MODIFICATION FOR C AND A GRIDS
! eight point relaxation using updated TTV to the lower and left
!      TTV(I,J)=AD05*(4.*(TTV(I-1,J)+TTV(I+1,J)
!     1                  +TTV(I,J-1)+TTV(I,J+1))
!     2                  +TTV(I-1,J-1)+TTV(I+1,J-1)
!     3                  +TTV(I-1,J+1)+TTV(I+1,J+1))
!     4                  -CFT0*TTV(I,J)
! eight point relaxation using old TTV
            a1=TTV(I-1,J)
            a2=TTV(I+1,J)
            a3=TTV(I,J-1)
            a4=TTV(I,J+1)
            a5=TTV(I-1,J-1)
            a6=TTV(I+1,J-1)
            a7=TTV(I-1,J+1)
            a8=TTV(I+1,J+1)
!            if ((a1-spval) <= 1e-10) a1=TTV(I,J)
!            if ((a2-spval) <= 1e-10) a2=TTV(I,J)
!            if ((a3-spval) <= 1e-10) a3=TTV(I,J)
!            if ((a4-spval) <= 1e-10) a4=TTV(I,J)
!            if ((a5-spval) <= 1e-10) a5=TTV(I,J)
!            if ((a6-spval) <= 1e-10) a6=TTV(I,J)
!            if ((a7-spval) <= 1e-10) a7=TTV(I,J)
!            if ((a8-spval) <= 1e-10) a8=TTV(I,J)

             if ((a1 < spval) .and.   &
                (a2 < spval) .and.   & 
                (a3 < spval) .and.  &
                (a4 < spval) .and.  & 
                (a5 < spval) .and.  &
                (a6 < spval) .and.  &
                (a7 < spval) .and.  &
                (a8 < spval) .and. (TTV(I,J) < spval))  then                    

!            TNEW(I,J) = AD05*(4.*(a1  +a2   +a3    &          
!                                 +a4) +a5 +a6  &        
!                                 +a7+a8)-TTV(I,J)*CFT0          

            TNEW(I,J) = AD05*(4.*(TTV(I-1,J)  +TTV(I+1,J)   +TTV(I,J-1)    &
                                 +TTV(I,J+1)) +TTV(I-1,J-1) +TTV(I+1,J-1)  &
                                 +TTV(I-1,J+1)+TTV(I+1,J+1))-TTV(I,J)*CFT0
            else 
            TNEW(I,J) = TTV(I,J)
            end if ! spval

! four point relaxation using old TTV
!      TNEW(I,J)=TTV(I,J)+1.0*((TTV(I-1,J)+TTV(I+1,J)
!     1          +TTV(I,J-1)+TTV(I,J+1)-4.0*TTV(I,J))/4.0)
! four point relaxation using updated TTV to the lower and left 
!      TTV(I,J)=TTV(I,J)+1.0*((TTV(I-1,J)+TTV(I+1,J)
!     1          +TTV(I,J-1)+TTV(I,J+1)-4.0*TTV(I,J))/4.0)
!
!     if(i==ii.and.j==jj)print*,'Debug: L,TTV A S'
!    1,l,TTV(I,J),N
!     1,l,TNEW(I,J),N

           end if ! spval

          enddo
!
!$omp parallel do private(i,j,km)
          DO KM=1,KMM
            I = IMNT(KM,L)
            J = JMNT(KM,L)
            if(PSLP(I,J)<spval .and. TNEW(I,J)< SPVAL/100.) then
              TTV(I,J) = TNEW(I,J)
            end if ! spval
          END DO
        END DO              ! NRLX loop
!
!$omp parallel do private(i,j,km)
        DO KM=1,KMM
          I = IMNT(KM,L)
          J = JMNT(KM,L)

          if(PSLP(I,J)<spval) then

! dong try to fix missing value for hgtprs at 1000 mb
          TPRES(I,J,L) = TTV(I,J)
          end if ! spval

!          if (QPRES(I,J,L) < 1000) TPRES(I,J,L) = TTV(I,J)
!          if (QPRES(I,J,L) < 1000) TPRES(I,J,L) = 1 

        END DO
      enddo          ! end of l loop
!----------------------------------------------------------------
!***
!***  CALCULATE THE SEA LEVEL PRESSURE AS PER THE NEW SCHEME.
!***  INTEGRATE THE HYDROSTATIC EQUATION DOWNWARD FROM THE
!***  GROUND THROUGH EACH OUTPUT PRESSURE LEVEL (WHERE TV
!***  IS NOW KNOWN) TO FIND GZ AT THE NEXT MIDPOINT BETWEEN
!***  PRESSURE LEVELS.  WHEN GZ=0 IS REACHED, SOLVE FOR THE
!***  PRESSURE.
!***
!
!***  BEFORE APPLYING RELAXATION FOR UNDERGROUND POINTS,
!***  FIRST FIND GRID POINTS AT/NEAR/BELOW SEA LEVEL AND DERIVE
!***  SEA LEVEL PRESSURE TO AVOID MEMBRANE RELAXATION
!***  AT THESE GRID POINTS.  E.G. HURRICANE CENTER NEAR COAST 
!
      KOUNT = 0
      DO J=JSTA,JEND
        DO I=1,IM

         if(PSLP(I,J)<spval) then
!         P1(I,J)=SPL(NINT(LMH(I,J)))
!         DONE(I,J)=.FALSE.

          IF(ABS(FIS(I,J)) < 1.) THEN
            PSLP(I,J) = PINT(I,J,NINT(LMH(I,J))+1)
            DONE(I,J) = .TRUE.
            KOUNT     = KOUNT + 1
!           if(i==ii.and.j==jj)print*,'Debug:DONE,PSLP A S1='        &  
!            ,done(i,j),PSLP(I,J)
          ELSE IF(FIS(I,J) < -1.0) THEN
            DO L=LM,1,-1
              IF(ZINT(I,J,L) > 0.)THEN
!               PSLP(I,J)=PINT(I,J,L)/EXP(-ZINT(I,J,L)*G                &
!               /(RD*T(I,J,L)*(Q(I,J,L)*D608+1.0)))

                tem = 0.5*(T(I,J,L)+T(I,J,L-1))*(1.0+0.5*D608*(Q(I,J,L)+Q(I,J,L-1)))
                PSLP(I,J) = PINT(I,J,L-1)/EXP(-ZINT(I,J,L-1)*G/(rd*tem))
                DONE(I,J) = .TRUE.
!               if(i==ii.and.j==jj)print*                           &
!               ,'Debug:DONE,PINT,PSLP A S1='                           &
!                ,done(i,j),PINT(I,J,L),PSLP(I,J)
                exit
              END IF
            END DO
          ENDIF

         end if ! spval

        ENDDO
      ENDDO
!
      KMM = KMNTM(LSM)
!!$omp parallel do private(gz1,gz2,i,j,lmap1,p1,p2),shared(pslp)
LOOP320:DO KM=1,KMM
        I = IMNT(KM,LSM)
        J = JMNT(KM,LSM)

        if(PSLP(I,J)<spval) then

        IF(DONE(I,J)) cycle
        LMHIJ   = LMHO(I,J)
        GZ1     = FIPRES(I,J,LMHIJ)
        P1(I,J) = SPL(LMHIJ)
!
        LMAP1 = LMHIJ+1
        DO L=LMAP1,LSM
          P2            = SPL(L)
          TLYR          = 0.5*(TPRES(I,J,L)+TPRES(I,J,L-1))
          GZ2           = GZ1 + RD*TLYR*LOG(P1(I,J)/P2)
          FIPRES(I,J,L) = GZ2
!         if(i==ii.and.j==jj)print*,'Debug:L,FI A S2=',L,GZ2
          IF(GZ2 <= 0.)THEN
            PSLP(I,J) = P1(I,J)/EXP(-GZ1/(RD*TPRES(I,J,L-1)))
!           if(i==ii.and.j==jj)print*,'Debug:PSLP A S2=',PSLP(I,J)
            DONE(I,J) = .TRUE.
            KOUNT     = KOUNT + 1
            CYCLE LOOP320            
          ENDIF
          P1(I,J) = P2
          GZ1     = GZ2
        ENDDO
!HC EXPERIMENT
        LP = LSM
        SLOPE     = -6.6E-4 
        TLYR      = TPRES(I,J,LP)-0.5*FIPRES(I,J,LP)*SLOPE
        PSLP(I,J) = spl(lp)/EXP(-FIPRES(I,J,LP)/(RD*TLYR))
        DONE(I,J) = .TRUE.
!     if(i==ii.and.j==jj)print*,'Debug:spl,FI,TLYR,PSLPA3='   &
!         ,spl(lp),FIPRES(I,J,LP),TLYR,PSLP(I,J)       
!HC EXPERIMENT
       end if ! spval

ENDDO LOOP320
 320  CONTINUE
!
!***  WHEN SEA LEVEL IS BELOW THE LOWEST OUTPUT PRESSURE LEVEL,
!***  SOLVE THE HYDROSTATIC EQUATION BY CHOOSING A TEMPERATURE
!***  AT THE MIDPOINT OF THE LAYER BETWEEN THAT LOWEST PRESSURE
!***  LEVEL AND THE GROUND BY EXTRAPOLATING DOWNWARD FROM T ON
!***  THE LOWEST PRESSURE LEVEL USING THE DT/DFI BETWEEN THE
!***  LOWEST PRESSURE LEVEL AND THE ONE ABOVE IT.
!
!      TOTAL=(IM-2)*(JM-4)
!
!HC      DO 340 LP=LSM,1,-1
!      IF(KOUNT==TOTAL)GO TO 350
!HC MODIFICATION FOR SMALL HILL HIGH PRESSURE SITUATION
!HC IF SURFACE PRESSURE IS CLOSER TO SEA LEVEL THAN LWOEST
!HC OUTPUT PRESSURE LEVEL, USE SURFACE PRESSURE TO DO EXTRAPOLATION

 325  CONTINUE 
      LP = LSM
      DO J=JSTA,JEND
        DO I=1,IM

         if(PSLP(I,J)<spval) then

!         if(i==ii.and.j==jj)print*,'Debug: with 330 loop'
          IF(DONE(I,J)) cycle

!         if(i==ii.and.j==jj)print*,'Debug: still within 330 loop'
!HC Comment out the following line for situation with terrain 
!HC at boundary (ie FIPRES<0)
!HC because they were not counted as undergound point for 8 pt
!HC relaxation
!HC      IF(FIPRES(I,J,LP)<0.)GO TO 330
!      IF(FIPRES(I,J,LP)<0.)THEN  
!       DO LP=LSM,1,-1
!        IF (FIPRES(I,J) <= 0)

!      IF(FIPRES(I,J,LP)<0..OR.DONE(I,J))GO TO 330
!     SLOPE=(TPRES(I,J,LP)-TPRES(I,J,LP-1))
!     & /(FIPRES(I,J,LP)-FIPRES(I,J,LP-1))     

          SLOPE = -6.6E-4
          IF(PINT(I,J,NINT(LMH(I,J))+1) > SPL(LP))THEN
            LLMH      = NINT(LMH(I,J))
            TVRT      = T(I,J,LLMH)*(H1+D608*Q(I,J,LLMH))
            DIS       = ZINT(I,J,LLMH+1)-ZINT(I,J,LLMH)+0.5*ZINT(I,J,LLMH+1)
            TLYR      = TVRT-DIS*G*SLOPE
            PSLP(I,J) = PINT(I,J,LLMH+1)*EXP(ZINT(I,J,LLMH+1)*G              &  
                         /(RD*TLYR))
!           if(i==ii.and.j==jj)print*,'Debug:PSFC,zsfc,TLYR,PSLPA3='
!           1,PINT(I,J,LLMH+1),ZINT(I,J,LLMH+1),TLYR,PSLP(I,J)
          ELSE
            TLYR=TPRES(I,J,LP)-0.5*FIPRES(I,J,LP)*SLOPE
            PSLP(I,J)=spl(lp)/EXP(-FIPRES(I,J,LP)/(RD*TLYR))
!           if(i==ii.and.j==jj)print*,'Debug:spl,FI,TLYR,PSLPA3='      &
!          ,spl(lp),FIPRES(I,J,LP),TLYR,PSLP(I,J)
          END IF
          DONE(I,J) = .TRUE.
          KOUNT     = KOUNT + 1
         end if ! spval

        enddo
      enddo
!HC  340 CONTINUE
!
! 350 CONTINUE
!--------------------------------------------------------------------
      RETURN
      END
