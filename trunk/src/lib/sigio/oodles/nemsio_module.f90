!<><><>    add baciof.f here (from /climate/save/wx20wa/gfsio/bacio/sorc)
!-----------------------------------------------------------------------
      MODULE BACIO_MODULE
!$$$  F90-MODULE DOCUMENTATION BLOCK
!
! F90-MODULE: BACIO_MODULE   BYTE-ADDRESSABLE I/O MODULE
!   PRGMMR: IREDELL          ORG: NP23        DATE: 98-06-04
!
! ABSTRACT: MODULE TO SHARE FILE DESCRIPTORS
!   IN THE BYTE-ADDRESSABLE I/O PACKAGE.
!
! PROGRAM HISTORY LOG:
!   98-06-04  IREDELL
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind
      use constants, only: ione
      implicit none

      INTEGER(i_kind),EXTERNAL:: BACIO,BACIOL
      INTEGER(i_kind),DIMENSION(999),SAVE:: FD=999*0
      INTEGER(i_kind),DIMENSION(20),SAVE:: BAOPTS=0
  !   INCLUDE 'baciof.h'
      !><><< insert bacio.f manually here (from /climate/save/wx20wa/gfsio/bacio/sorc)
!     Include file to define variables for Fortran to C interface(s)
!     Robert Grumbine 16 March 1998
      INTEGER(i_kind),PARAMETER:: BACIO_OPENR=ione   ! Open file for read only
      INTEGER(i_kind),PARAMETER:: BACIO_OPENW=2      ! Open file for write only
      INTEGER(i_kind),PARAMETER:: BACIO_OPENRW=4     ! Open file for read or write
      INTEGER(i_kind),PARAMETER:: BACIO_CLOSE=8      ! Close file
      INTEGER(i_kind),PARAMETER:: BACIO_READ=16      ! Read from the file
      INTEGER(i_kind),PARAMETER:: BACIO_WRITE=32     ! Write to the file
      INTEGER(i_kind),PARAMETER:: BACIO_NOSEEK=64    ! Start I/O from previous spot
      INTEGER(i_kind),PARAMETER:: BACIO_OPENWT=128   ! Open for write only with truncation
      INTEGER(i_kind),PARAMETER:: BACIO_OPENWA=256   ! Open for write only with append
      END
!-----------------------------------------------------------------------

      SUBROUTINE BASETO(NOPT,VOPT)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BASETO         BYTE-ADDRESSABLE SET OPTIONS
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: SET OPTIONS FOR BYTE-ADDRESSABLE I/O.
!   ALL OPTIONS DEFAULT TO 0.
!   OPTION 1: BLOCKED READING OPTION
!             IF THE OPTION VALUE IS 1, THEN THE READING IS BLOCKED
!             INTO FOUR 4096-BYTE BUFFERS.  THIS MAY BE EFFICIENT IF
!             THE READS WILL BE REQUESTED IN MUCH SMALLER CHUNKS.
!             OTHERWISE, EACH CALL TO BAREAD INITIATES A PHYSICAL READ.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BASETO(NOPT,VOPT)
!   INPUT ARGUMENTS:
!     NOPT         INTEGER OPTION NUMBER
!     VOPT         INTEGER OPTION VALUE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind
      use constants, only: ione
      USE BACIO_MODULE
      implicit none

      INTEGER(i_kind), intent(in) :: NOPT,VOPT
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(NOPT.GE.ione.AND.NOPT.LE.20) BAOPTS(NOPT)=VOPT
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
!-----------------------------------------------------------------------
      SUBROUTINE BAOPEN(LU,CFN,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAOPEN         BYTE-ADDRESSABLE OPEN
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: OPEN A BYTE-ADDRESSABLE FILE.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BAOPEN(LU,CFN,IRET)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO OPEN
!     CFN          CHARACTER FILENAME TO OPEN
!                  (CONSISTING OF NONBLANK PRINTABLE CHARACTERS)
!   OUTPUT ARGUMENTS:
!     IRET         INTEGER RETURN CODE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$

      use kinds, only: i_kind,i_llong
      use constants, only: ione
      USE BACIO_MODULE
      implicit none
      integer(i_kind), intent(in) :: LU
      CHARACTER, intent(in) :: CFN*(*)
      integer(i_kind), intent(out) :: IRET
      CHARACTER(80) CMSG
      integer(i_llong) IB,JB,NB,KA
      CHARACTER A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LU.LT.001.OR.LU.GT.999) THEN
        IRET=6
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IRET=BACIOL(BACIO_OPENRW,IB,JB,ione,NB,KA,FD(LU),CFN,A)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
!-----------------------------------------------------------------------

      SUBROUTINE BAOPENR(LU,CFN,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAOPENR        BYTE-ADDRESSABLE OPEN
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: OPEN A BYTE-ADDRESSABLE FILE FOR READ ONLY.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BAOPENR(LU,CFN,IRET)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO OPEN
!     CFN          CHARACTER FILENAME TO OPEN
!                  (CONSISTING OF NONBLANK PRINTABLE CHARACTERS)
!   OUTPUT ARGUMENTS:
!     IRET         INTEGER RETURN CODE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: ione
      USE BACIO_MODULE
      implicit none
      CHARACTER, intent(in) :: CFN*(*)
      INTEGER(i_kind), intent(in) :: LU
      INTEGER(i_kind), intent(out) :: iret
      integer(i_llong) IB,JB,NB,KA
      CHARACTER A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LU.LT.001.OR.LU.GT.999) THEN
      print *,'in baopenr, lu=',lu
        IRET=6
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      print *,'in baopenr, before call baciol'
      IRET=BACIOL(BACIO_OPENR,IB,JB,ione,NB,KA,FD(LU),CFN,A)
      print *,'in baopenr, iret=',iret
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
!-----------------------------------------------------------------------

      SUBROUTINE BAOPENW(LU,CFN,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAOPENW        BYTE-ADDRESSABLE OPEN
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: OPEN A BYTE-ADDRESSABLE FILE FOR WRITE ONLY.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BAOPENW(LU,CFN,IRET)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO OPEN
!     CFN          CHARACTER FILENAME TO OPEN
!                  (CONSISTING OF NONBLANK PRINTABLE CHARACTERS)
!   OUTPUT ARGUMENTS:
!     IRET         INTEGER RETURN CODE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: ione
      USE BACIO_MODULE
      implicit none
      CHARACTER, intent(in) :: CFN*(*)
      INTEGER(i_kind), intent(in) :: LU
      INTEGER(i_kind), intent(out) :: iret
      integer(i_llong) IB,JB,NB,KA
      CHARACTER A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LU.LT.001.OR.LU.GT.999) THEN
        IRET=6
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IRET=BACIOL(BACIO_OPENW,IB,JB,ione,NB,KA,FD(LU),CFN,A)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
!-----------------------------------------------------------------------

      SUBROUTINE BAOPENWT(LU,CFN,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAOPENWT       BYTE-ADDRESSABLE OPEN
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: OPEN A BYTE-ADDRESSABLE FILE FOR WRITE ONLY WITH TRUNCATION.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BAOPENWT(LU,CFN,IRET)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO OPEN
!     CFN          CHARACTER FILENAME TO OPEN
!                  (CONSISTING OF NONBLANK PRINTABLE CHARACTERS)
!   OUTPUT ARGUMENTS:
!     IRET         INTEGER RETURN CODE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: ione
      USE BACIO_MODULE
      implicit none
      CHARACTER, intent(in) :: CFN*(*)
      INTEGER(i_kind), intent(in) :: LU
      INTEGER(i_kind), intent(out) :: iret
      integer(i_llong) IB,JB,NB,KA
      CHARACTER A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LU.LT.001.OR.LU.GT.999) THEN
        IRET=6
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IRET=BACIOL(BACIO_OPENWT,IB,JB,ione,NB,KA,FD(LU),CFN,A)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
!-----------------------------------------------------------------------

      SUBROUTINE BAOPENWA(LU,CFN,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAOPENWA       BYTE-ADDRESSABLE OPEN
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: OPEN A BYTE-ADDRESSABLE FILE FOR WRITE ONLY WITH APPEND.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BAOPENWA(LU,CFN,IRET)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO OPEN
!     CFN          CHARACTER FILENAME TO OPEN
!                  (CONSISTING OF NONBLANK PRINTABLE CHARACTERS)
!   OUTPUT ARGUMENTS:
!     IRET         INTEGER RETURN CODE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: ione
      USE BACIO_MODULE
      implicit none
      CHARACTER, intent(in) :: CFN*(*)
      INTEGER(i_kind), intent(in) :: LU
      INTEGER(i_kind), intent(out) :: iret
      integer(i_llong) IB,JB,NB,KA
      CHARACTER A
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LU.LT.001.OR.LU.GT.999) THEN
        IRET=6
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IRET=BACIOL(BACIO_OPENWA,IB,JB,ione,NB,KA,FD(LU),CFN,A)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
!-----------------------------------------------------------------------

      SUBROUTINE BACLOSE(LU,IRET)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BACLOSE        BYTE-ADDRESSABLE CLOSE
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: CLOSE A BYTE-ADDRESSABLE FILE.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BACLOSE(LU,IRET)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO CLOSE
!   OUTPUT ARGUMENTS:
!     IRET         INTEGER RETURN CODE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! REMARKS:  A BAOPEN MUST HAVE ALREADY BEEN CALLED.
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: izero,ione
      USE BACIO_MODULE
      implicit none
      INTEGER(i_kind), intent(in) :: LU
      INTEGER(i_kind), intent(out) :: iret
      integer(i_llong) IB,JB,NB,KA
      CHARACTER A
      CHARACTER CFN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(LU.LT.001.OR.LU.GT.999) THEN
        IRET=6
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IRET=BACIOL(BACIO_CLOSE,IB,JB,ione,NB,KA,FD(LU),CFN,A)
      IF(IRET.EQ.izero) FD(LU)=izero
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END
!-----------------------------------------------------------------------

      SUBROUTINE BAREAD(LU,IB,NB,KA,A)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    baread
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
        INTEGER(i_kind),INTENT(IN) :: LU,IB,NB
        INTEGER(i_kind),INTENT(OUT) :: KA
        CHARACTER,INTENT(OUT) :: A(NB)
        INTEGER(i_llong) :: LONG_IB,LONG_NB,LONG_KA
!
        if(IB<izero .or. NB<izero ) THEN
          print *,'WRONG: in BAFRREAD starting postion IB or read '//    &
         'data size NB < 0, STOP! '//                                    &
         'Consider using bafreadl and long integer'
          KA=izero
          return
        ENDIF
        LONG_IB=IB
        LONG_NB=NB
        CALL BAREADL(LU,LONG_IB,LONG_NB,LONG_KA,A)
        KA=LONG_KA

      END SUBROUTINE BAREAD
!-----------------------------------------------------------------------

      SUBROUTINE BAREADL(LU,IB,NB,KA,A)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAREAD         BYTE-ADDRESSABLE READ
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: READ A GIVEN NUMBER OF BYTES FROM AN UNBLOCKED FILE,
!   SKIPPING A GIVEN NUMBER OF BYTES.
!   THE PHYSICAL I/O IS BLOCKED INTO FOUR 4096-BYTE BUFFERS
!   IF THE BYTE-ADDRESSABLE OPTION 1 HAS BEEN SET TO 1 BY BASETO.
!   THIS BUFFERED READING IS INCOMPATIBLE WITH NO-SEEK READING.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BAREAD(LU,IB,NB,KA,A)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO READ
!     IB           INTEGER NUMBER OF BYTES TO SKIP
!                  (IF IB<0, THEN THE FILE IS ACCESSED WITH NO SEEKING)
!     NB           INTEGER NUMBER OF BYTES TO READ
!   OUTPUT ARGUMENTS:
!     KA           INTEGER NUMBER OF BYTES ACTUALLY READ
!     A            CHARACTER*1 (NB) DATA READ
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! REMARKS:  A BAOPEN MUST HAVE ALREADY BEEN CALLED.
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: izero,ione
      USE BACIO_MODULE
!    
      IMPLICIT NONE
      INTEGER(i_kind),intent(in)  :: LU
      INTEGER(i_llong),intent(in)  :: IB,NB
      INTEGER(i_llong),intent(out) :: KA
      CHARACTER,intent(out)       :: A(NB)
      CHARACTER CFN
      integer(i_llong),PARAMETER :: NY=4096,MY=4
      INTEGER(i_llong) NS(MY),NN(MY)
      INTEGER(i_llong) JB,LONG_0,KY,I,K,IY,JY,LUX
      INTEGER(i_kind) IRET
!      INTEGER LU,IB,NB,KA
      CHARACTER Y(NY,MY)
      DATA LUX/0/
      SAVE JY,NS,NN,Y,LUX
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(FD(LU).LE.izero) THEN
        KA=izero
        RETURN
      ENDIF
      IF(IB.LT.izero.AND.BAOPTS(1).EQ.ione) THEN
        KA=izero
        RETURN
      ENDIF
      IF(NB.LE.izero) THEN
        KA=izero
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LONG_0=izero                                                         !jw
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  UNBUFFERED I/O
      IF(BAOPTS(1).NE.ione) THEN
        IF(IB.GE.izero) THEN
          IRET=BACIOL(BACIO_READ,IB,JB,ione,NB,KA,FD(LU),CFN,A)
        ELSE
!          IRET=BACIOL(BACIO_READ+BACIO_NOSEEK,izero,JB,ione,NB,KA,FD(LU),CFN,A)
          IRET=BACIOL(BACIO_READ+BACIO_NOSEEK,LONG_0,JB,ione,NB,KA,        &
                      FD(LU),CFN,A)
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  BUFFERED I/O
!  GET DATA FROM PREVIOUS CALL IF POSSIBLE
      ELSE
        KA=izero
        IF(LUX.NE.LU) THEN
          JY=izero
          NS=izero
          NN=izero
        ELSE
          DO I=1,MY
            IY=MOD(JY+I-ione,MY)+ione
            KY=IB+KA-NS(IY)
            IF(KA.LT.NB.AND.KY.GE.LONG_0.AND.KY.LT.NN(IY)) THEN
              K=MIN(NB-KA,NN(IY)-KY)
              A(KA+ione:KA+K)=Y(KY+ione:KY+K,IY)
              KA=KA+K
            ENDIF
          ENDDO
        ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  SET POSITION AND READ BUFFER AND GET DATA
        IF(KA.LT.NB) THEN
          LUX=ABS(LU)
          JY=MOD(JY,MY)+ione
          NS(JY)=IB+KA
          IRET=BACIOL(BACIO_READ,NS(JY),JB,ione,NY,NN(JY),  &
                     FD(LUX),CFN,Y(1,JY))
          IF(NN(JY).GT.izero) THEN
            K=MIN(NB-KA,NN(JY))
            A(KA+ione:KA+K)=Y(1:K,JY)
            KA=KA+K
          ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  CONTINUE TO READ BUFFER AND GET DATA
          DO WHILE(NN(JY).EQ.NY.AND.KA.LT.NB)
            JY=MOD(JY,MY)+ione
            NS(JY)=NS(JY)+NN(JY)
            IRET=BACIOL(BACIO_READ+BACIO_NOSEEK,NS(JY),JB,ione,NY,NN(JY),  &
                       FD(LUX),CFN,Y(1,JY))
            IF(NN(JY).GT.izero) THEN
              K=MIN(NB-KA,NN(JY))
              A(KA+ione:KA+K)=Y(1:K,JY)
              KA=KA+K
            ENDIF
          ENDDO
        ENDIF
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE BAREADL
!-----------------------------------------------------------------------

      SUBROUTINE BAWRITE(LU,IB,NB,KA,A)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    bawrite
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-09-01  lueken - added subprogram doc block
!
!   input argument list:
!    LU,IB,NB
!    A
!
!   output argument list:
!    KA
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
        use kinds, only: i_kind,i_llong
        use constants, only: izero
        IMPLICIT NONE
        INTEGER(i_kind),INTENT(IN) :: LU,IB,NB
        INTEGER(i_kind),INTENT(OUT) :: KA
        CHARACTER,INTENT(IN) :: A(NB)
        INTEGER(i_llong) :: LONG_IB,LONG_NB,LONG_KA
!
        if(IB<izero .or. NB<izero ) THEN
          print *,'WRONG: in BAFRWRITEstarting postion IB or read '//    &
         'data size NB <0, STOP! ' //                                    &
         'Consider using bafrrwritel and long integer'
          KA=izero
          return
        ENDIF
!
        LONG_IB=IB
        LONG_NB=NB
        CALL BAWRITEL(LU,LONG_IB,LONG_NB,LONG_KA,A)
        KA=LONG_KA                                                     !jun

      END SUBROUTINE BAWRITE
!-----------------------------------------------------------------------

      SUBROUTINE BAWRITEL(LU,IB,NB,KA,A)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: BAWRITE        BYTE-ADDRESSABLE WRITE
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: WRITE A GIVEN NUMBER OF BYTES TO AN UNBLOCKED FILE,
!   SKIPPING A GIVEN NUMBER OF BYTES.
!
! PROGRAM HISTORY LOG:
!   1998-06-04  IREDELL
!
! USAGE:    CALL BAWRITE(LU,IB,NB,KA,A)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO WRITE
!     IB           INTEGER NUMBER OF BYTES TO SKIP
!                  (IF IB<0, THEN THE FILE IS ACCESSED WITH NO SEEKING)
!     NB           INTEGER NUMBER OF BYTES TO WRITE
!     A            CHARACTER*1 (NB) DATA TO WRITE
!   OUTPUT ARGUMENTS:
!     KA           INTEGER NUMBER OF BYTES ACTUALLY WRITTEN
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! REMARKS:  A BAOPEN MUST HAVE ALREADY BEEN CALLED.
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: izero,ione
      USE BACIO_MODULE
!
      IMPLICIT NONE
!
      INTEGER(i_kind),intent(in) :: LU
      INTEGER(i_llong),intent(in) :: IB,NB
      INTEGER(i_llong),intent(out):: KA
      CHARACTER,intent(in) ::  A(NB)
!
      CHARACTER CFN
      INTEGER(i_llong) :: JB,LONG_0
      INTEGER(i_kind) :: IRET
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(FD(LU).LE.izero) THEN
        KA=izero
        RETURN
      ENDIF
      IF(NB.LE.izero) THEN
        KA=izero
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LONG_0=izero
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(IB.GE.izero) THEN
        IRET=BACIOL(BACIO_WRITE,IB,JB,ione,NB,KA,FD(LU),CFN,A)
      ELSE
!        IRET=BACIOL(BACIO_WRITE+BACIO_NOSEEK,izero,JB,ione,NB,KA,FD(LU),CFN,A)
        IRET=BACIOL(BACIO_WRITE+BACIO_NOSEEK,LONG_0,JB,ione,NB,KA,         &
                    FD(LU),CFN,A)
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      END SUBROUTINE  BAWRITEL
!-----------------------------------------------------------------------

      SUBROUTINE WRYTE(LU,NB,A)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    wryte
!   prgmmr:
!
! abstract:
!
! program history log:
!   2009-09-01  lueken - added subprogram doc block
!
!   input argument list:
!    LU
!    NB
!    A
!
!   output argument list:
!
! attributes:
!   language: f90
!   machine:
!
!$$$ end documentation block
      use kinds, only: i_kind,i_llong
      use constants, only: izero
      USE BACIO_MODULE
!
      IMPLICIT NONE
!
      INTEGER(i_kind),intent(in) :: LU
      INTEGER(i_kind),intent(in) :: NB
      CHARACTER,intent(in) ::  A(NB)
      INTEGER(i_llong) :: LONG_NB
!
      IF(NB<izero) THEN
       PRINT *,'WRONG: NB: the number of bytes to write  <0, STOP!'
       RETURN
      ENDIF
      LONG_NB=NB
      print *,'in wryte,nb=',nb
      CALL WRYTEL(LU,LONG_NB,A)
!
      END SUBROUTINE WRYTE
!-----------------------------------------------------------------------

      SUBROUTINE WRYTEL(LU,NB,A)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: WRYTE          WRITE DATA OUT BY BYTES
!   PRGMMR: IREDELL          ORG: W/NMC23     DATE: 1998-06-04
!
! ABSTRACT: WRITE A GIVEN NUMBER OF BYTES TO AN UNBLOCKED FILE.
!
! PROGRAM HISTORY LOG:
!   92-10-31  IREDELL
!   95-10-31  IREDELL     WORKSTATION VERSION
!   1998-06-04  IREDELL   BACIO VERSION
!
! USAGE:    CALL WRYTE(LU,NB,A)
!   INPUT ARGUMENTS:
!     LU           INTEGER UNIT TO WHICH TO WRITE
!     NB           INTEGER NUMBER OF BYTES TO WRITE
!     A            CHARACTER*1 (NB) DATA TO WRITE
!
! MODULES USED:
!   BACIO_MODULE   BYTE-ADDRESSABLE I/O FORTRAN INTERFACE
!
! SUBPROGRAMS CALLED:
!   BACIO          BYTE-ADDRESSABLE I/O C PACKAGE
!
! REMARKS:  A BAOPEN MUST HAVE ALREADY BEEN CALLED.
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN 90
!
!$$$
      use kinds, only: i_kind,i_llong
      use constants, only: izero,ione
      USE BACIO_MODULE
!
      IMPLICIT NONE
      INTEGER(i_kind),intent(in) :: LU
      INTEGER(i_llong),intent(in) :: NB
      CHARACTER,INTENT(in)       :: A(NB)
      INTEGER(i_llong) :: LONG_0,JB,KA
      INTEGER(i_kind) :: IRET
      CHARACTER CFN
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      IF(FD(LU).LE.izero) THEN
        RETURN
      ENDIF
      IF(NB.LE.izero) THEN
        RETURN
      ENDIF
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      LONG_0=izero
      IRET=BACIOL(BACIO_WRITE+BACIO_NOSEEK,LONG_0,JB,ione,NB,KA,           &
                  FD(LU),CFN,A)
      print *,'in wrytel,nb=',nb,'jb=',jb,'ka=',ka,'iret=',iret
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END

!<><><>    add bafrio.f here (from /climate/save/wx20wa/gfsio/bacio/sorc)
!-----------------------------------------------------------------------

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
