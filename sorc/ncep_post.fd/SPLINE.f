!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SPLINE(JTB,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)
!     ******************************************************************
!     * *
!     * THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE        *
!     * PROGRAMED FOR A SMALL SCALAR MACHINE.                         *
!     * *
!     * PROGRAMER Z. JANJIC, YUGOSLAV FED. HYDROMET. INST., BEOGRAD   *
!     * *
!     * *
!     * NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION.  MUST BE GE 3. *
!     * XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
!     * FUNCTION ARE GIVEN.  MUST BE IN ASCENDING ORDER.       *
!     * YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD.   *
!     * Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD.  IF NATURAL *
!     * SPLINE IS FITTED Y2(1)=0. AND Y2(NOLD)=0. MUST BE      *
!     * SPECIFIED.                                             *
!     * NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.     *
!     * XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
!     * FUNCTION ARE CALCULATED.  XNEW(K) MUST BE GE XOLD(1)   *
!     * AND LE XOLD(NOLD).                                     *
!     * YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED.           *
!     * P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2.                *
!     * *
!     ******************************************************************
!
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      implicit none
!
      integer,intent(in) :: JTB,NOLD,NNEW
      real,dimension(JTB),intent(in) ::  XOLD,YOLD,XNEW 
      real,dimension(JTB),intent(inout) :: P,Q,Y2
      real,dimension(JTB),intent(out) ::  YNEW
!
      integer NOLDM1,K,K1,K2,KOLD
      real DXL,DXR,DYDXL,DYDXR,RTDXC,DXC,DEN,XK,Y2K,Y2KP1,DX,RDX,       &
            Ak,BK,CK,X,XSQ
!-----------------------------------------------------------------------
      NOLDM1=NOLD-1
!
      DXL=XOLD(2)-XOLD(1)
      DXR=XOLD(3)-XOLD(2)
      DYDXL=(YOLD(2)-YOLD(1))/DXL
      DYDXR=(YOLD(3)-YOLD(2))/DXR
      RTDXC=.5/(DXL+DXR)
!
      P(1)= RTDXC*(6.*(DYDXR-DYDXL)-DXL*Y2(1))
      Q(1)=-RTDXC*DXR
!
!-----------------------------------------------------------------------
      DO K=3, NOLD-1
      DXL=DXR
      DYDXL=DYDXR
      DXR=XOLD(K+1)-XOLD(K)
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR
      DXC=DXL+DXR
      DEN=1./(DXL*Q(K-2)+DXC+DXC)
!
      P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))
      Q(K-1)=-DEN*DXR
      END DO
!-----------------------------------------------------------------------
      DO K=NOLDM1, 2, -1
      Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)
      END DO
!-----------------------------------------------------------------------
      K=0
      DO K1=1,NNEW
      XK=XNEW(K1)
!
      KOLD=0
      DO K2=2,NOLD
      IF(XOLD(K2)>XK) THEN
      KOLD=K2-1
      EXIT
      ENDIF
      END DO
!
      IF(KOLD==0) THEN
      YNEW(K1)=YOLD(NOLD)
      CYCLE
      ENDIF
!
      IF(K/=KOLD) THEN
      K=KOLD
!
      Y2K=Y2(K)
      Y2KP1=Y2(K+1)
      DX=XOLD(K+1)-XOLD(K)
      RDX=1./DX
!
      AK=.1666667*RDX*(Y2KP1-Y2K)
      BK=.5*Y2K
      CK=RDX*(YOLD(K+1)-YOLD(K))-.1666667*DX*(Y2KP1+Y2K+Y2K)
      ENDIF
!
      X=XK-XOLD(K)
      XSQ=X*X
!
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)
      END DO
!-----------------------------------------------------------------------
                              RETURN
                              END
