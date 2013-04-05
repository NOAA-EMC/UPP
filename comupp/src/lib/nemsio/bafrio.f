!-----------------------------------------------------------------------
      SUBROUTINE BAFRINDEX(LU,IB,LX,IX,do_byteswap)
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LU,IB
      INTEGER,INTENT(INOUT):: LX
      INTEGER,INTENT(OUT):: IX
      logical,intent(in) :: do_byteswap
      integer(kind=8) :: LONG_JB,LONG_JX,LONG_IX
!
      call BAFRINDEXL(LU,LONG_JB,LONG_JX,LONG_IX,do_byteswap)
      LX=LONG_JX
      IX=LONG_IX
      return
      end SUBROUTINE BAFRINDEX
!-----------------------------------------------------------------------
      SUBROUTINE BAFRINDEXL(LU,IB,LX,IX,do_byteswap)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAFRINDEX      BYTE-ADDRESSABLE FORTRAN RECORD INDEX
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1999-01-21
!
! ABSTRACT: THIS SUBPROGRAM EITHER READS AN UNFORMATTED FORTRAN RECORD
!   AND RETURN ITS LENGTH AND START BYTE OF THE NEXT FORTRAN RECORD;
!   OR GIVEN THE RECORD LENGTH, WITHOUT I/O IT DETERMINES THE START BYTE
!   OF THE NEXT FORTRAN RECORD.
!
! PROGRAM HISTORY LOG:
!   1999-01-21  IREDELL
!
! USAGE:    CALL BAFRINDEX(LU,IB,LX,IX)
!   INPUT ARGUMENTS:
!     LU           INTEGER LOGICAL UNIT TO READ
!                  IF LU<=0, THEN DETERMINE IX FROM LX
!     IB           INTEGER FORTRAN RECORD START BYTE
!                  (FOR THE FIRST FORTRAN RECORD, IB SHOULD BE 0)
!     LX           INTEGER RECORD LENGTH IN BYTES IF LU<=0
!
!   OUTPUT ARGUMENTS:
!     LX           INTEGER RECORD LENGTH IN BYTES IF LU>0,
!                  OR LX=-1 FOR I/O ERROR (PROBABLE END OF FILE),
!                  OR LX=-2 FOR I/O ERROR (INVALID FORTRAN RECORD)
!     IX           INTEGER START BYTE FOR THE NEXT FORTRAN RECORD
!                  (COMPUTED ONLY IF LX>=0)
!
! SUBPROGRAMS CALLED:
!   BAREAD         BYTE-ADDRESSABLE READ
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LU
      INTEGER(KIND=8),INTENT(IN):: IB
      INTEGER(KIND=8),INTENT(INOUT):: LX
      INTEGER(KIND=8),INTENT(OUT):: IX
      logical,intent(in) :: do_byteswap
      INTEGER(KIND=8),PARAMETER:: LBCW=4
      INTEGER(KIND=LBCW):: BCW1,BCW2
      INTEGER(KIND=8):: KR
      INTEGER(4) LBCW2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPARE FIRST BLOCK CONTROL WORD AND TRAILING BLOCK CONTROL WORD
      IF(LU.GT.0) THEN
        CALL BAREADL(LU,IB,LBCW,KR,BCW1)
!LLF+MP--
        if(do_byteswap) then
          LBCW2=LBCW
          CALL byteswap(BCW1,LBCW2,1)
        endif
!LLF+MP==
        IF(KR.NE.LBCW) THEN
          LX=-1
        ELSE
          CALL BAREADL(LU,IB+LBCW+BCW1,LBCW,KR,BCW2)
!LLF+MP--
          if(do_byteswap) CALL byteswap(BCW2,LBCW2,1)
!LLF+MP==
          IF(KR.NE.LBCW.OR.BCW1.NE.BCW2) THEN
            LX=-2
          ELSE
            LX=BCW1
          ENDIF
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE START BYTE FOR THE NEXT FORTRAN RECORD
      IF(LX.GE.0) IX=IB+LBCW+LX+LBCW
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE BAFRINDEXL
!-----------------------------------------------------------------------
      SUBROUTINE BAFRREAD(LU,IB,NB,KA,A,do_byteswap)
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LU,IB,NB
      INTEGER,INTENT(OUT):: KA
      CHARACTER,INTENT(OUT):: A(NB)
      logical,intent(in) :: do_byteswap
      INTEGER(KIND=8) :: LONG_IB,LONG_NB,LONG_KA
!
        if(IB<0 .or. NB<0 ) THEN
          print *,'WRONG: in BAFRREAD starting postion IB or read '//    &
     & 'data size NB < 0, STOP! Consider use bafreadl and long integer'
          KA=0
          return
        ENDIF
        LONG_IB=IB
        LONG_NB=NB
        CALL BAFRREADL(LU,LONG_IB,LONG_NB,LONG_KA,A,do_byteswap)
        KA=LONG_KA
      END SUBROUTINE BAFRREAD
!-----------------------------------------------------------------------
      SUBROUTINE BAFRREADL(LU,IB,NB,KA,A,do_byteswap)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAFRREAD       BYTE-ADDRESSABLE FORTRAN RECORD READ
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1999-01-21
!
! ABSTRACT: THIS SUBPROGRAM READS AN UNFORMATTED FORTRAN RECORD
!
! PROGRAM HISTORY LOG:
!   1999-01-21  IREDELL
!
! USAGE:    CALL BAFRREAD(LU,IB,NB,KA,A)
!   INPUT ARGUMENTS:
!     LU           INTEGER LOGICAL UNIT TO READ
!     IB           INTEGER FORTRAN RECORD START BYTE
!                  (FOR THE FIRST FORTRAN RECORD, IB SHOULD BE 0)
!     NB           INTEGER NUMBER OF BYTES TO READ
!
!   OUTPUT ARGUMENTS:
!     KA           INTEGER NUMBER OF BYTES IN FORTRAN RECORD
!                  (IN WHICH CASE THE NEXT FORTRAN RECORD
!                  SHOULD HAVE A START BYTE OF IB+KA),
!                  OR KA=-1 FOR I/O ERROR (PROBABLE END OF FILE),
!                  OR KA=-2 FOR I/O ERROR (INVALID FORTRAN RECORD),
!                  OR KA=-3 FOR I/O ERROR (REQUEST LONGER THAN RECORD)
!     A            CHARACTER*1 (NB) DATA READ
!
! SUBPROGRAMS CALLED:
!   BAFRINDEX      BYTE-ADDRESSABLE FORTRAN RECORD INDEX
!   BAREAD         BYTE-ADDRESSABLE READ
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LU
      INTEGER(kind=8),INTENT(IN):: IB,NB
      INTEGER(kind=8),INTENT(OUT):: KA
      CHARACTER,INTENT(OUT):: A(NB)
      logical,intent(in) :: do_byteswap
      INTEGER(kind=8),PARAMETER:: LBCW=4
      INTEGER(kind=8):: LX,IX
      INTEGER(kind=8):: KR
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  VALIDATE FORTRAN RECORD
      CALL BAFRINDEXL(LU,IB,LX,IX,do_byteswap)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  READ IF VALID
      IF(LX.LT.0) THEN
        KA=LX
      ELSEIF(LX.LT.NB) THEN
        KA=-3
      ELSE
        CALL BAREADL(LU,IB+LBCW,NB,KR,A)
        IF(KR.NE.NB) THEN
          KA=-1
        ELSE
          KA=LBCW+LX+LBCW
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE BAFRREADL
!-----------------------------------------------------------------------
      SUBROUTINE BAFRWRITE(LU,IB,NB,KA,A,do_byteswap)
!
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LU,IB,NB
      INTEGER,INTENT(OUT):: KA
      CHARACTER,INTENT(OUT):: A(NB)
      logical,intent(in) :: do_byteswap
      INTEGER(KIND=8) :: LONG_IB,LONG_NB,LONG_KA
!
        if(IB<0 .or. NB<0 ) THEN
          print *,'WRONG: in BAFRREAD starting postion IB or read '//    &
     &   'data size NB <0, STOP! ' //                                    &
     &   'Consider use bafrrwritel and long integer'
          KA=0
          return
        ENDIF
        LONG_IB=IB
        LONG_NB=NB
        CALL BAFRWRITEL(LU,LONG_IB,LONG_NB,LONG_KA,A,do_byteswap)
        KA=LONG_KA
!
      END SUBROUTINE BAFRWRITE
!-----------------------------------------------------------------------
      SUBROUTINE BAFRWRITEL(LU,IB,NB,KA,A,do_byteswap)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAFRWRITE      BYTE-ADDRESSABLE FORTRAN RECORD WRITE
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1999-01-21
!
! ABSTRACT: THIS SUBPROGRAM WRITES AN UNFORMATTED FORTRAN RECORD
!
! PROGRAM HISTORY LOG:
!   1999-01-21  IREDELL
!
! USAGE:    CALL BAFRWRITE(LU,IB,NB,KA,A)
!   INPUT ARGUMENTS:
!     LU           INTEGER LOGICAL UNIT TO WRITE
!     IB           INTEGER FORTRAN RECORD START BYTE
!                  (FOR THE FIRST FORTRAN RECORD, IB SHOULD BE 0)
!     NB           INTEGER NUMBER OF BYTES TO WRITE
!     A            CHARACTER*1 (NB) DATA TO WRITE
!
!   OUTPUT ARGUMENTS:
!     KA           INTEGER NUMBER OF BYTES IN FORTRAN RECORD
!                  (IN WHICH CASE THE NEXT FORTRAN RECORD
!                  SHOULD HAVE A START BYTE OF IB+KA),
!                  OR KA=-1 FOR I/O ERROR
!
! SUBPROGRAMS CALLED:
!   BAWRITE        BYTE-ADDRESSABLE WRITE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      IMPLICIT NONE
      INTEGER,INTENT(IN):: LU
      INTEGER(KIND=8),INTENT(IN):: IB,NB
      INTEGER(kind=8),INTENT(OUT):: KA
      CHARACTER,INTENT(IN):: A(NB)
      logical,intent(in) :: do_byteswap
      INTEGER(kind=8),PARAMETER:: LBCW=4
      INTEGER(kind=LBCW):: BCW
      INTEGER(kind=8):: KR
!LLF+MP--
      INTEGER(LBCW)::BCW2,LBCW2
!LLF+MP==
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  WRITE DATA BRACKETED BY BLOCK CONTROL WORDS
      BCW=NB
!      print *,'in writel,bcw=',bcw,'do_byteswap=',do_byteswap
!LLF+MP--
      if(do_byteswap) then
        LBCW2=LBCW
        CALL byteswap(BCW,LBCW2,1)
      endif
!LLF+MP==
      CALL BAWRITEL(LU,IB,LBCW,KR,BCW)
!LLF+MP--
      if(do_byteswap) CALL byteswap(BCW,LBCW2,1)
!LLF+MP==
      IF(KR.NE.LBCW) THEN
        KA=-1
      ELSE
        CALL BAWRITEL(LU,IB+LBCW,NB,KR,A)
        IF(KR.NE.NB) THEN
          KA=-1
        ELSE
!LLF+MP--
          BCW2=BCW
          if(do_byteswap) CALL byteswap(BCW,LBCW2,1)
!LLF+MP==
          CALL BAWRITEL(LU,IB+LBCW+BCW2,LBCW,KR,BCW)
!          CALL BAWRITEL(LU,IB+LBCW+BCW,LBCW,KR,BCW)
!LLF+MP--
          if(do_byteswap) CALL byteswap(BCW,LBCW2,1)
!LLF+MP==
          IF(KR.NE.LBCW) THEN
            KA=-1
          ELSE
            KA=LBCW+BCW+LBCW
          ENDIF
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE  BAFRWRITEL

