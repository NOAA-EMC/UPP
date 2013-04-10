#if LINUX
!----------------------------------------------------------------------
! From Numerical Recipes
! added two Numerical Recipes routines for matrix inversion
! LUBKSB - solves a system of linear equations, follows call to LUDCMP
! LUDCMP - replaces an NxN matrix a with the LU decomposition
! 
!  2012-11-05 E.Mirvis separated this generic LU from the splat.F
!----------------------------------------------------------------------
      SUBROUTINE LUBKSB(A,N,NP,INDX,B)
      REAL A(NP,NP),B(N)
      INTEGER INDX(N)
      II=0  
      DO 12 I=1,N
        LL=INDX(I)
        SUM=B(LL)
        B(LL)=B(I)
        IF (II.NE.0)THEN
          DO 11 J=II,I-1
            SUM=SUM-A(I,J)*B(J)
   11     CONTINUE
        ELSE IF (SUM.NE.0.) THEN
          II=I
        ENDIF
        B(I)=SUM
   12 CONTINUE
      DO 14 I=N,1,-1  
        SUM=B(I)
        IF(I.LT.N)THEN
          DO 13 J=I+1,N
            SUM=SUM-A(I,J)*B(J)
   13 CONTINUE
        ENDIF
        B(I)=SUM/A(I,I)
   14 CONTINUE
      RETURN
      END

      SUBROUTINE LUDCMP(A,N,NP,INDX)
C      PARAMETER (NMAX=400,TINY=1.0E-20)
      PARAMETER (TINY=1.0E-20)
C==EM==^^^
C
      REAL A(NP,NP),VV(N),D
C     REAL A(NP,NP),VV(NMAX),D
C==EM==^^^
      INTEGER INDX(N)
      D=1.  
      DO 12 I=1,N
        AAMAX=0.
        DO 11 J=1,N
          IF (ABS(A(I,J)).GT.AAMAX) AAMAX=ABS(A(I,J))
   11   CONTINUE
        IF (AAMAX.EQ.0.) print *, 'SINGULAR MATRIX.'
        VV(I)=1./AAMAX
   12 CONTINUE
      DO 19 J=1,N
        IF (J.GT.1) THEN
          DO 14 I=1,J-1
            SUM=A(I,J)
            IF (I.GT.1)THEN
              DO 13 K=1,I-1
                SUM=SUM-A(I,K)*A(K,J)
   13         CONTINUE
              A(I,J)=SUM
            ENDIF
   14     CONTINUE
        ENDIF
        AAMAX=0.
        DO 16 I=J,N
          SUM=A(I,J)
          IF (J.GT.1)THEN
            DO 15 K=1,J-1
              SUM=SUM-A(I,K)*A(K,J)
   15       CONTINUE  
            A(I,J)=SUM
          ENDIF
          DUM=VV(I)*ABS(SUM)
          IF (DUM.GE.AAMAX) THEN
            IMAX=I
            AAMAX=DUM
          ENDIF
   16   CONTINUE
        IF (J.NE.IMAX)THEN
          DO 17 K=1,N
            DUM=A(IMAX,K)
            A(IMAX,K)=A(J,K)
            A(J,K)=DUM
   17     CONTINUE
          D=-D
          VV(IMAX)=VV(J)
        ENDIF
        INDX(J)=IMAX
        IF(J.NE.N)THEN
          IF(A(J,J).EQ.0.)A(J,J)=TINY
          DUM=1./A(J,J)
          DO 18 I=J+1,N
            A(I,J)=A(I,J)*DUM
   18     CONTINUE
        ENDIF
   19 CONTINUE
      IF(A(N,N).EQ.0.)A(N,N)=TINY
      RETURN
      END
#endif
