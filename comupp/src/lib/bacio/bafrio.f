      SUBROUTINE BAFRINDEX(LU,IB,LX,IX)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    bafrindex
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-09-01  lueken - added subprogram doc block
!
!   input argument list:
!    LU,IB
!    LX
!
!   output argument list:
!    LX
!    IX
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
      use kinds, only: i_kind,i_llong
      IMPLICIT NONE
      INTEGER(i_kind),INTENT(IN):: LU,IB
      INTEGER(i_kind),INTENT(INOUT):: LX
      INTEGER(i_kind),INTENT(OUT):: IX
      integer(i_llong) :: LONG_IB,LONG_LX ,LONG_IX
!
      LONG_IB=IB
      LONG_LX=LX
      call BAFRINDEXL(LU,LONG_IB,LONG_LX,LONG_IX)
      LX=LONG_LX
      IX=LONG_IX

      return
      end SUBROUTINE BAFRINDEX
!-----------------------------------------------------------------------

      SUBROUTINE BAFRINDEXL(LU,IB,LX,IX)
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
      use kinds, only: i_kind,i_llong
      use constants, only: izero,ione
      IMPLICIT NONE
      INTEGER(i_kind),INTENT(IN):: LU
      INTEGER(i_llong),INTENT(IN):: IB
      INTEGER(i_llong),INTENT(INOUT):: LX
      INTEGER(i_llong),INTENT(OUT):: IX
      INTEGER(i_llong),PARAMETER:: LBCW=4
      INTEGER(i_kind):: BCW1,BCW2
      INTEGER(i_llong):: KR
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPARE FIRST BLOCK CONTROL WORD AND TRAILING BLOCK CONTROL WORD
      IF(LU.GT.izero) THEN
        CALL BAREADL(LU,IB,LBCW,KR,BCW1)
        IF(KR.NE.LBCW) THEN
          LX=-ione
        ELSE
          CALL BAREADL(LU,IB+LBCW+BCW1,LBCW,KR,BCW2)
          IF(KR.NE.LBCW.OR.BCW1.NE.BCW2) THEN
            LX=-2
          ELSE
            LX=BCW1
          ENDIF
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  COMPUTE START BYTE FOR THE NEXT FORTRAN RECORD
      IF(LX.GE.izero) IX=IB+LBCW+LX+LBCW
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE BAFRINDEXL
!-----------------------------------------------------------------------

      SUBROUTINE BAFRREAD(LU,IB,NB,KA,A)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    bafrread
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-09-01  lueken - added subprogram doc block
!
!   input argument list:
!    LU,IB,NB
!
!   output argument list:
!    KA
!    A
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
      use kinds, only: i_kind,i_llong
      use constants, only: izero
      IMPLICIT NONE
      INTEGER(i_kind),INTENT(IN):: LU,IB,NB
      INTEGER(i_kind),INTENT(OUT):: KA
      CHARACTER,INTENT(OUT):: A(NB)
      INTEGER(i_llong) :: LONG_IB,LONG_NB,LONG_KA
!
        if(IB<izero .or. NB<izero ) THEN
          print *,'WRONG: in BAFRREAD starting postion IB or read '//    &
       'data size NB < 0, STOP! Consider use bafreadl and long integer'
          KA=izero
          return
        ENDIF
        LONG_IB=IB
        LONG_NB=NB
        CALL BAFRREADL(LU,LONG_IB,LONG_NB,LONG_KA,A)
        KA=LONG_KA
      END SUBROUTINE BAFRREAD
!-----------------------------------------------------------------------

      SUBROUTINE BAFRREADL(LU,IB,NB,KA,A)
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
      use kinds, only: i_kind,i_llong
      use constants, only: izero,ione
      IMPLICIT NONE
      INTEGER(i_kind),INTENT(IN):: LU
      INTEGER(i_llong),INTENT(IN):: IB,NB
      INTEGER(i_llong),INTENT(OUT):: KA
      CHARACTER,INTENT(OUT):: A(NB)
      INTEGER(i_llong),PARAMETER:: LBCW=4
      INTEGER(i_llong):: LX,IX
      INTEGER(i_llong):: KR
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  VALIDATE FORTRAN RECORD
      CALL BAFRINDEXL(LU,IB,LX,IX)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  READ IF VALID
      IF(LX.LT.izero) THEN
        KA=LX
      ELSEIF(LX.LT.NB) THEN
        KA=-3
      ELSE
        CALL BAREADL(LU,IB+LBCW,NB,KR,A)
        IF(KR.NE.NB) THEN
          KA=-ione
        ELSE
          KA=LBCW+LX+LBCW
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE BAFRREADL
!-----------------------------------------------------------------------

      SUBROUTINE BAFRWRITE(LU,IB,NB,KA,A)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    bafrwrite
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-09-01  lueken - added subprogram doc block
!
!   input argument list:
!    LU,IB,NB
!
!   output argument list:
!    KA
!    A
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
      use kinds, only: i_kind,i_llong
      use constants, only: izero
      IMPLICIT NONE
      INTEGER(i_kind),INTENT(IN):: LU,IB,NB
      INTEGER(i_kind),INTENT(OUT):: KA
      CHARACTER,INTENT(OUT):: A(NB)
      INTEGER(i_llong) :: LONG_IB,LONG_NB,LONG_KA
!
        if(IB<izero .or. NB<izero ) THEN
          print *,'WRONG: in BAFRREAD starting postion IB or read '//    &
         'data size NB <0, STOP! ' //                                    &
         'Consider use bafrrwritel and long integer'
          KA=izero
          return
        ENDIF
        LONG_IB=IB
        LONG_NB=NB
        CALL BAFRWRITEL(LU,LONG_IB,LONG_NB,LONG_KA,A)
        KA=LONG_KA
!
      END SUBROUTINE BAFRWRITE
!-----------------------------------------------------------------------

      SUBROUTINE BAFRWRITEL(LU,IB,NB,KA,A)
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
      use kinds, only: i_kind,i_llong
      use constants, only: ione
      IMPLICIT NONE
      INTEGER(i_kind),INTENT(IN):: LU
      INTEGER(i_llong),INTENT(IN):: IB,NB
      INTEGER(i_llong),INTENT(OUT):: KA
      CHARACTER,INTENT(IN):: A(NB)
      INTEGER(i_llong),PARAMETER:: LBCW=4
      INTEGER(i_kind):: BCW
      INTEGER(i_llong):: KR
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  WRITE DATA BRACKETED BY BLOCK CONTROL WORDS
      BCW=NB
      CALL BAWRITEL(LU,IB,LBCW,KR,BCW)
      IF(KR.NE.LBCW) THEN
        KA=-ione
      ELSE
        CALL BAWRITEL(LU,IB+LBCW,NB,KR,A)
        IF(KR.NE.NB) THEN
          KA=-ione
        ELSE
          CALL BAWRITEL(LU,IB+LBCW+BCW,LBCW,KR,BCW)
          IF(KR.NE.LBCW) THEN
            KA=-ione
          ELSE
            KA=LBCW+BCW+LBCW
          ENDIF
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE  BAFRWRITEL
