!> @file
!> @brief Subroutine that computes dewpoints from vapor pressure.
!>
!> This routine is to  computes the dewpoints for the N values
!> of vapor pressure in array VP.
!> The forumla:
!>
!> VP = 0.611 * (X**A) * EXP( (A+B)*(1-X) )
!>
!> is used to get dewpoint temperature T, where
!>
!> X = T3/T,                   T3=Triple PT temperature,
!> VP=Vapor pressure in CBS,   0.611=VP at T3,
!> A=(Spec. HT. of WATER-CSUBP of vapor)/gas const of vapor
!>                           and
!> B=Latent heat at T3/(gas const of vapor times T3).
!>   
!> on the first call, a table TDP is constructed giving
!> dewpoint as a function of vapor pressure.
!>  
!> Values of VP less than the first table entry
!> (RVP1 in the code) will be given dewpoints for
!> that beginning valus. Similarly, VP vaules that
!> exceed the maximum table value (RVP2 in the code)
!> will be assigned dewpoints for that maximum value.
!>   
!> The values 0.02 and 8.0 for RVP1 and RVP2 yield
!> dewpoints of 233.6K and 314.7K,respectively.
!>
!> @param[in] VP Array of N vapor pressures(centibars).
!> @param[out] TD Dewpoint in degrees absolute.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1990-05-19 | Jim Tuccillo | Initial
!> 1993-05-12 | R Treadon    | Expanded table size and reset range of pressures covered by table.
!> 1998-06-12 | T Black      | Conversion from 1-D to 2-D
!> 2000-01-04 | Jim Tuccillo | MPI Version
!> 2021-07-26 | W Meng       | Restrict computation from undefined grids
!>
!> @author Jim Tuccillo W/NP2 @date 1990-05-19
      SUBROUTINE DEWPOINT( VP, TD)

       use ctlblk_mod, only: jsta, jend, im, spval
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
!          NT IS THE TABLE SIZE
      integer,PARAMETER :: NT=2000
!...TRANSLATED BY FPP 3.00Z36 11/09/90  14:48:53  
!...SWITCHES: OPTON=I47,OPTOFF=VAE0
      real,intent(out) :: TD(IM,jsta:jend)
      real,intent(in) ::  VP(IM,jsta:jend)
      real TDP(NT)
!jw
      integer NN,I,J,JNT
      real rvp1,rvp2,rt3,rvp3,rlog3,ra,rb,rapb,rtest,rnt,rdvp
      real rgs,rvp,rlvp,rn,rd,rch,rt,w1,w2
      real A,B,DNTM1

      logical :: jcontinue=.true.

!          PREPARE THE TABLE (TDP=DEWPT AS FCN OF VAPOR PRESS).
!          RANGE IN CENTIBARS IS FROM RVP1 THRU RVP2
      rvp1  = 0.0001E0
      rvp2  = 10.E0
!          THE TRIPLE POINT
      RT3   = 273.16E0
!          VAPOR PRESS AT THE TRIPLE POINT
      RVP3  = 0.611E0
      RLOG3 = LOG(RVP3)
!          (SPEC HT OF WATER -CSUBP OF VAPOR)/GAS CONST OF VAPOR.
      RA    = 5.0065E0
!          LATENT HEAT AT T3/(GAS CONST OF VAPOR * TRIPLE PT TEMP).
      RB    = 19.83923E0
      RAPB  = RA + RB
!          CRITERION FOR CONVERGENCE OF NEWTON ITERATION
      RTEST = 1.E-6
!MEB  RTEST=1.E-8  !  PROBABLY WON'T CONVERGE WITH 32-BIT AT THIS CRITERION
!
      RNT   = FLOAT(NT)
!          TABLE INCREMENT IN VAPOR PRESS
      RDVP  = (RVP2-RVP1)/(RNT-1.E0)
!          RGS WILL BE THE GUESSED VALUE OF (T3  /  DEWPOINT)
      RGS   = 1.E0
      RVP   = RVP1-RDVP
!
      DO 20 NN=1,NT
        RVP=RVP+RDVP
        RLVP=LOG(RVP)-RLOG3-RAPB
!     ***** ENTER NEWTON ITERATION LOOP
        jcontinue=.true.
        do while (jcontinue)
   10   RN=RA*LOG(RGS)-RAPB*RGS-RLVP
!          THAT WAS VALUE OF FUNCTION
!          NOW GET ITS DERIVATIVE
        RD=(RA/RGS)-RAPB
!          THE DESIRED CHANGE IN THE GUESS
        RCH=RN/RD
        IF( ABS(RCH) < RTEST ) jcontinue=.false.
!            NEED MORE ITERATIONS
          DO WHILE (ABS(RCH) >= RTEST)
            RGS=RGS-RCH
            EXIT
          ENDDO
        ENDDO
!          *****
!          HAVE ACCURATE ENUF VALUE OF RGS=T3/DEWPOINT.
   15   RT=RT3/RGS
        TDP(NN)=RT
!
   20 CONTINUE
!      PRINT 25,RVP1,RVP2,TDP(1),TDP(NT)
!  25  FORMAT(/'0', 'IN SUBROUTINE DEWPOINT, THE DEWPT TABLE ',
!    1             'HAS RVP1=', 1PE13.6, ', RVP2=', 1PE13.6,
!    2             ', TDP(1)=', 1PE13.6, ', AND TDP(NT)=',
!    3             1PE13.6, '.'/)
!           CONSTANTS FOR USING THE TABLE
      A     = 1./RDVP
      B     = 1. - A*RVP1
      DNTM1 = FLOAT(NT) -.01
!
!X      END IF
!
!          *********** ENTER TO USE THE TABLE.  ************
!
!$omp parallel do private(i,j,w1,w2,jnt)
      DO J=JSTA,JEND
        DO I=1,IM
          IF(VP(I,J)<spval)THEN
          W1  = MIN(MAX((A*VP(I,J)+B),1.0),DNTM1)
          W2  = AINT(W1)
          JNT = INT(W2)
          TD(I,J) = TDP(JNT) + (W1-W2)*(TDP(JNT+1)-TDP(JNT))
          ELSE
          TD(I,J) = spval
          ENDIF
        ENDDO
      ENDDO
!
!
      RETURN
      END
