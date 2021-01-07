      subroutine get_postfilename(fname)
!
! ABSTRACT: THIS SUBROUTINE GENERATE POST FILE NAME FROM THE DATSET IN
!           POST CONTROL FILE
!
!  Program log:
!     11-02        Jun Wang   generate code from subroutine gribit
!
      use ctlblk_mod,  only : ifhr, me, modelname, ifmin
      use rqstfld_mod, only : ritehd, datset, iget
!
      implicit none
!
      character(*),intent(inout) :: fname
!
!local vars
      integer IHR,KDAT,KENV,KTHR,NDIG
      CHARACTER*4  RESTHR,BLANK
      CHARACTER*10  DESCR2,DESCR3
      character CFHOUR*40,CFORM*40
      CHARACTER*180 ENVAR
      CHARACTER*255 PGBOUT,IPVOUT,D3DOUT
!
      DATA BLANK /'    '/
!
      IF (RITEHD) THEN
!
!        PUT FORECAST HOUR INTO DIR PREFIX FOR GRIB FILE.
         IHR = IFHR
!
!        GET FULL PATH FOR OUTPUT FILE FROM ENVIRONMENT VARIABLE
!        COMSP WHICH IS SET IN THE SCRIPT RUNNING THE MODEL.
!
!        CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE
         ENVAR = ' '
         RESTHR = ' '
         PGBOUT = ' '
         IPVOUT = ' '
         D3DOUT = ' '
         CALL GETENV('COMSP',ENVAR)
         CALL GETENV('tmmark',RESTHR)
         CALL GETENV('PGBOUT',PGBOUT)
         CALL GETENV('IPVOUT',IPVOUT)
         CALL GETENV('D3DOUT',D3DOUT)
         KDAT = INDEX(DATSET,' ') -1
         IF (KDAT<=0) KDAT = LEN(DATSET)
         KENV = INDEX(ENVAR,' ') -1
         IF (KENV<=0) KENV = LEN(ENVAR)
         KTHR = INDEX(RESTHR,' ') -1
         IF (KTHR<=0) KTHR = LEN(RESTHR)
         if(me==0) print *,'PGBOUT=',trim(PGBOUT)
!
      if(me==0)print *,'in get postfilename, ritehd=',ritehd,'ifhr=',ifhr,'modelname=',modelname, &
        'ENVAR(1:4)=',ENVAR(1:4),'RESTHR(1:4)=',RESTHR(1:4),'ifmin=',ifmin,'DATSET(1:KDAT)=', &
        DATSET(1:KDAT)
!
!        CONSTRUCT FULL PATH-FILENAME FOR OUTPUT FILE
         IF(MODELNAME=='GFS')THEN
          IF(D3DOUT(1:4)/=BLANK .AND. &
             ((IGET(354)>0).OR.(IGET(355)>0).OR.  &
             (IGET(356)>0).OR.(IGET(357)>0).OR.  &
             (IGET(358)>0).OR.(IGET(359)>0).OR.  &
             (IGET(360)>0).OR.(IGET(361)>0).OR.  &
             (IGET(362)>0).OR.(IGET(363)>0).OR.  &
             (IGET(364)>0).OR.(IGET(365)>0).OR.  &
             (IGET(366)>0).OR.(IGET(367)>0).OR.  &
             (IGET(368)>0).OR.(IGET(369)>0).OR.  &
             (IGET(370)>0).OR.(IGET(371)>0).OR.  &
             (IGET(372)>0).OR.(IGET(373)>0).OR.  &
             (IGET(374)>0).OR.(IGET(375)>0)))THEN
              FNAME = D3DOUT
              if(me==0)PRINT*,' FNAME FROM D3DOUT=',trim(FNAME)
          ELSE IF(IPVOUT(1:4)/=BLANK .AND.           &
              index(DATSET(1:KDAT),"IPV")>0 .AND.  &
             ((IGET(332)>0).OR.(IGET(333)>0).OR.  &
             (IGET(334)>0).OR.(IGET(335)>0).OR.  &
             (IGET(351)>0).OR.(IGET(352)>0).OR.  &
             (IGET(353)>0).OR.(IGET(378)>0)))THEN
              FNAME = IPVOUT
              if(me==0)PRINT*,' FNAME FROM IPVOUT=',trim(FNAME)
          ELSE IF(PGBOUT(1:4)/=BLANK)THEN
            FNAME = PGBOUT
            if(me==0)PRINT*,' FNAME FROM PGBOUT=',trim(FNAME)
          ELSE
              NDIG=MAX(LOG10(IHR+0.5)+1.,2.)
!          WRITE(CFORM,'("('.GrbF',I",I1,".",I1,")")') NDIG,NDIG
              WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
              WRITE(CFHOUR,CFORM) IHR
              FNAME = DATSET(1:KDAT) //'.GrbF'// CFHOUR
              if(me==0)print *,' FNAME=',trim(FNAME)
          END IF
!         IF(MODELNAME=='GFS'.AND.PGBOUT(1:4)/=BLANK)THEN
!          FNAME = PGBOUT
!          PRINT*,' FNAME FROM PGBOUT=',trim(FNAME)
!
         ELSEIF (ENVAR(1:4)==BLANK.AND.RESTHR(1:4)==BLANK) THEN
          IF(IFMIN >= 1)THEN
           WRITE(DESCR2,1011) IHR
           WRITE(DESCR3,1012) IFMIN
           FNAME = DATSET(1:KDAT) // TRIM(DESCR2)  //'.'// DESCR3(1:2)
          ELSE
           NDIG=MAX(LOG10(IHR+0.5)+1.,2.)
!          WRITE(CFORM,'("('.GrbF',I",I1,".",I1,")")') NDIG,NDIG
           WRITE(CFORM,'("(I",I1,".",I1,")")') NDIG,NDIG
           WRITE(CFHOUR,CFORM) IHR
           FNAME = DATSET(1:KDAT) //'.GrbF'// CFHOUR
           if(me==0)print *,' FNAME=',trim(FNAME)
!
!          IF(IHR<100)THEN
!           WRITE(DESCR2,1011) IHR
!          ELSE
!           WRITE(DESCR2,1013) IHR
!          END IF
 1011      FORMAT('.GrbF',I2.2)
!1013      FORMAT('.GrbF',I3.3)
!          FNAME = DATSET(1:KDAT) // DESCR2
          END IF
!
         ELSEIF(ENVAR(1:4)==BLANK.AND.RESTHR(1:4)/=BLANK) THEN
          IF(IFMIN >= 1)THEN
           WRITE(DESCR3,1012) IFMIN
           IF (IHR<100) THEN
              WRITE(DESCR2,1012) IHR
              FNAME = DATSET(1:KDAT) // DESCR2(1:2)  //'.'// DESCR3(1:2) &
                 //'.'// RESTHR
           ELSE
              WRITE(DESCR2,1014) IHR
              FNAME = DATSET(1:KDAT) // DESCR2(1:3)  //'.'// DESCR3(1:2) &
                 //'.'// RESTHR
           ENDIF
          ELSE
           IF (IHR<100) THEN
             WRITE(DESCR2,1012) IHR
             FNAME = DATSET(1:KDAT) // DESCR2(1:2)  //'.'// RESTHR
           ELSE
             WRITE(DESCR2,1014) IHR
             FNAME = DATSET(1:KDAT) // DESCR2(1:3)  //'.'// RESTHR
           ENDIF
          end if
         ELSE
          IF(IFMIN >= 1)THEN
           WRITE(DESCR3,1012) IFMIN
           IF (IHR<100) THEN
             WRITE(DESCR2,1012) IHR
             FNAME = ENVAR(1:KENV) // DATSET(1:KDAT) // DESCR2(1:2)  &
             //'.'// DESCR3(1:2) //'.'// RESTHR
           ELSE
             WRITE(DESCR2,1014) IHR
             FNAME = ENVAR(1:KENV) // DATSET(1:KDAT) // DESCR2(1:3)  &
             //'.'// DESCR3(1:2) //'.'// RESTHR
           ENDIF
          ELSE
           IF (IHR<100) THEN
             WRITE(DESCR2,1012) IHR
             FNAME = ENVAR(1:KENV) // DATSET(1:KDAT) // DESCR2(1:2) &
                    //'.'// RESTHR
 1012        FORMAT(I2.2)
 1014        FORMAT(I3.3)
           ELSE
             WRITE(DESCR2,1014) IHR
             FNAME = ENVAR(1:KENV) // DATSET(1:KDAT) // DESCR2(1:3) &
                    //'.'// RESTHR
           ENDIF
          end if
         ENDIF
!
      ENDIF
      if(me==0) then
        print*,'FNAME= ',trim(FNAME)
        print *,'end of get post filename'
      endif

      end subroutine get_postfilename 
