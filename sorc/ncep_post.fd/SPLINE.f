!
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
      SUBROUTINE SPLINE(JTB,NOLD,XOLD,YOLD,Y2,NNEW,XNEW,YNEW,P,Q)
!     ******************************************************************
!     *                                                                *
!     *  THIS IS A ONE-DIMENSIONAL CUBIC SPLINE FITTING ROUTINE        *
!     *  PROGRAMED FOR A SMALL SCALAR MACHINE.                         *
!     *                                                                *
!     *  PROGRAMER Z. JANJIC, YUGOSLAV FED. HYDROMET. INST., BEOGRAD  *
!     *                                                                *
!     *                                                                *
!     *                                                                *
!     *  NOLD - NUMBER OF GIVEN VALUES OF THE FUNCTION.  MUST BE GE 3. *
!     *  XOLD - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
!     *         FUNCTION ARE GIVEN.  MUST BE IN ASCENDING ORDER.       *
!     *  YOLD - THE GIVEN VALUES OF THE FUNCTION AT THE POINTS XOLD.   *
!     *  Y2   - THE SECOND DERIVATIVES AT THE POINTS XOLD.  IF NATURAL *
!     *         SPLINE IS FITTED Y2(1)=0. AND Y2(NOLD)=0. MUST BE      *
!     *         SPECIFIED.                                             *
!     *  NNEW - NUMBER OF VALUES OF THE FUNCTION TO BE CALCULATED.     *
!     *  XNEW - LOCATIONS OF THE POINTS AT WHICH THE VALUES OF THE     *
!     *         FUNCTION ARE CALCULATED.  XNEW(K) MUST BE GE XOLD(1)   *
!     *         AND LE XOLD(NOLD).                                     *
!     *  YNEW - THE VALUES OF THE FUNCTION TO BE CALCULATED.           *
!     *  P, Q - AUXILIARY VECTORS OF THE LENGTH NOLD-2.                *
!     *                                                                *
!     ******************************************************************
!
! PROGRAM HISTORY LOG:
!   19-10-30  Bo CUI - REMOVE "GOTO" STATEMENT
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
      real DXL,DXR,DYDXL,DYDXR,RTDXC,DXC,DEN,XK,Y2K,Y2KP1,DX,RDX,      &
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
  
      loop700: do
!     IF(NOLD.EQ.3) GO TO 700
      IF(NOLD.EQ.3) exit loop700
!-----------------------------------------------------------------------
      K=3
!
 loop100: do
 100  DXL=DXR
      DYDXL=DYDXR
      DXR=XOLD(K+1)-XOLD(K)
      DYDXR=(YOLD(K+1)-YOLD(K))/DXR
      DXC=DXL+DXR
      DEN=1./(DXL*Q(K-2)+DXC+DXC)
!
      P(K-1)= DEN*(6.*(DYDXR-DYDXL)-DXL*P(K-2))
      Q(K-1)=-DEN*DXR
!
      K=K+1
!     IF(K.LT.NOLD) GO TO 100
      IF(K.LT.NOLD) cycle loop100
      exit loop100
      enddo loop100
!-----------------------------------------------------------------------
      exit loop700
      enddo loop700
 700  K=NOLDM1
!
 loop200: do
 200  Y2(K)=P(K-1)+Q(K-1)*Y2(K+1)
!
      K=K-1
!     IF(K.GT.1) GO TO 200
      IF(K.GT.1) cycle loop200
      exit loop200
      enddo loop200
!-----------------------------------------------------------------------
      K1=1
!
 loop300: do
 300  XK=XNEW(K1)
!
      loop600: do
      loop550: do
      loop500: do
      loop450: do
      DO 400 K2=2,NOLD
!     IF(XOLD(K2).LE.XK) GO TO 400
      IF(XOLD(K2).LE.XK) cycle         
      KOLD=K2-1
!     GO TO 450
      exit loop450
 400  CONTINUE
      YNEW(K1)=YOLD(NOLD)
!     GO TO 600
      exit loop600
!
      exit loop450
      enddo loop450
!450  IF(K1.EQ.1)   GO TO 500
!     IF(K.EQ.KOLD) GO TO 550
 450  IF(K1.EQ.1)   exit loop500
      IF(K.EQ.KOLD) exit loop550
!
      exit loop500
      enddo loop500
 500  K=KOLD
!
      Y2K=Y2(K)
      Y2KP1=Y2(K+1)
      DX=XOLD(K+1)-XOLD(K)
      RDX=1./DX
!
!VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
!     WRITE(6,5000) K,Y2K,Y2KP1,DX,RDX,YOLD(K),YOLD(K+1)
!5000 FORMAT(' K=',I4,' Y2K=',E12.4,' Y2KP1=',E12.4,' DX=',E12.4,' RDX='
!    2,E12.4,' YOK=',E12.4,' YOP1=',E12.4)
!AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
      AK=.1666667*RDX*(Y2KP1-Y2K)
      BK=.5*Y2K
      CK=RDX*(YOLD(K+1)-YOLD(K))-.1666667*DX*(Y2KP1+Y2K+Y2K)
!
      exit loop550
      enddo loop550
 550  X=XK-XOLD(K)
      XSQ=X*X
!
      YNEW(K1)=AK*XSQ*X+BK*XSQ+CK*X+YOLD(K)
!
      exit loop600
      enddo loop600
 600  K1=K1+1
!     IF(K1.LE.NNEW) GO TO 300
      IF(K1.LE.NNEW) cycle loop300
      exit loop300
      enddo loop300
!-----------------------------------------------------------------------
                              RETURN
                              END
