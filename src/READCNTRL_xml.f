      SUBROUTINE READCNTRL_xml(kth,kpv,pv)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    READCNTRL  READS CONTROL FILE
!   PRGRMMR: TREADON         ORG: W/NP2      DATE: 92-12-20       
!     
! ABSTRACT:
!     THIS ROUTINE READS THE CONTROL FILE SPECIFYING
!     DATA FORMAT(S) AND FIELD(S) TO POST.  THE
!     ORDER OF OPERATIONS IS 
!        (1) READ HEADER BLOCK OF CONTROL FILE,
!        (2) SET FLAGS, CLOSE OPEN UNITS
!        (3) READ BODY OF CONTROL FILE (FIELD SPECIFICATIONS)
!   .     
!     
! PROGRAM HISTORY LOG:
!   03-31-2010  Jun Wang - GRIB2 VERSION
!     
! USAGE:    CALL READCNTRL(IEOF)
!   INPUT ARGUMENT LIST:
!     NONE
!
!   OUTPUT ARGUMENT LIST: 
!     IEOF     - INTEGER FLAG FOR EOF IN CONTROL FILE.
!                IEOF=0 WHEN AN EOF IS READ IN THE
!                CONTROL FILE.  IEOF=1 OTHERWISE.
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!
!     LIBRARY:
!       COMMON   - RQSTFLD 
!                  CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : CRAY C-90
!$$$  
!
!     
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
       use xml_data_postcntrl_t
       use grib2_module, only: num_pset,pset,nrecout,grib_info_init
       use lookup_mod,only: ITB,JTB,ITBQ,JTBQ
       use ctlblk_mod
       use rqstfld_mod
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     DECLARE VARIABLES.
!     
      LOGICAL NORTH
      CHARACTER*2  CHAR2
      CHARACTER*4  CHAR4
      CHARACTER*80 LINE
      REAL EGRID1(IM,JM), EGRID2(IM,JM)
!
      integer, intent(in) :: KTH,KPV
      real,intent(in) :: pv(kpv)
!
      integer ISUM,L,IFLD,MFLD,IAVBL,IEOF,IREC,I,J,IND1,IND2
      CHARACTER*50 AVBLGRB_NAME
      logical :: FOUND_FLD,LSTATS,LSIGMA,first_grbtbl,LHYBPRES
!
!******************************************************************************
!     START READCNTRL HERE.
!     
      IF(ME.EQ.0)THEN
        WRITE(6,*)'READCNTRL:  POSTING FCST HR ',IFHR,' FROM ',         &
             IHRST,'UTC ',SDAT(1),'-',SDAT(2),'-',SDAT(3),' RUN'
      ENDIF
!     
!     INITIALIZE VARIABLES.
!        IEOF IS THE END OF FILE FLAG FOR THE CONTROL FILE.
!        ARRAY IGET IS THE "GET FIELD" FLAG ARRAY.
!
      IEOF=0
      DO 100 IFLD=1,MXFLD
        IGET(IFLD)=-1
 100  CONTINUE
!
!     SET FLAG TO OPEN NEW OUTPUT FILE
!
      RITEHD = .TRUE.
      DO J=1,size(LVLS,2)
      DO I=1,size(LVLS,1)
        LVLS(I,J)=0
        LVLSXML(I,J)=0
      ENDDO
      ENDDO
!
!     READ post cntrl file
      first_grbtbl=.false.
      if(npset==1) then
        call read_xml_file_postcntrl_t( 'postcntrl.xml')
        num_pset=size(paramset)
        first_grbtbl=.true.
      endif
      pset=paramset(npset)
      datset=pset%datset
      print *,'in readxml, num_pset=',num_pset,'datset=',trim(pset%datset)
! 
!     NOW READ WHICH FIELDS ON 
!     WHICH LEVELS TO INTERPOLATE TO THE OUTPUT GRID.  THE
!     CHARACTER STRING "DONE" MARKS THE END OF THE OUTPUT
!     FIELD SPECIFICATIONS.
!
      call grib_info_init(first_grbtbl)
      MFLD = size(pset%param)
      print *,'start reading control file,MFLD=',MFLD,'datset=',datset,MXFLD
      IFLD=1
      irec=0
      DO I=1, MFLD
         
!     SEE IF REQUESTED FIELD IS AVAILABLE.  IF NOT, 
!     WRITE MESSAGE TO 6 AND DECREMENT FIELD 
!     COUNTER BY ONE.  THEN READ NEXT REQUESTED FIELD.
!     
      doavbl:   DO 20 IAVBL = 1,MXFLD
             
            IND1=INDEX(trim(AVBLGRB2(IAVBL))," ON ")
            AVBLGRB_NAME=AVBLGRB2(IAVBL)(1:IND1-1)
            IND2=INDEX(trim(AVBLGRB_NAME)," ",back=.true.)
            IF(IND2>0) THEN
              AVBLGRB_NAME=AVBLGRB2(IAVBL)(IND2+1:IND1-1)
            ENDIF
!            if(me==0.and.AVBLGRB_NAME=="TMP".and.trim(paramset(npset)%param(i)%pname)=='TMP'  &
!             .and.trim(paramset(npset)%param(i)%fixed_sfc1_type)=='depth_bel_land_sfc') &
!              print *,'fldname=',trim(paramset(npset)%param(i)%pname),trim(AVBLGRB2(IAVBL)),   &
!             'level=',trim(paramset(npset)%param(i)%fixed_sfc1_type),   &
!             'i=',i,'iavbl=',iavbl,'ind1=',ind1,'ind2=',ind2,'ind3=',ind3
 
            IF (trim(AVBLGRB_NAME)==trim(paramset(npset)%param(i)%pname).and.    &
              INDEX(trim(AVBLGRB2(IAVBL)),trim(paramset(npset)%param(i)%fixed_sfc1_type))/=0) then

             LSTATS=trim(paramset(npset)%param(i)%stats_proc)/=''

             LSIGMA=trim(paramset(npset)%param(i)%fixed_sfc1_type)=='sigma_lvl' .and. &
                    size(paramset(npset)%param(i)%level)==1 .and. &
                    nint(paramset(npset)%param(i)%level(1))==9950
!
	     LHYBPRES=trim(paramset(npset)%param(i)%fixed_sfc2_type)=='hybrid_lvl'.and. &
                    (paramset(npset)%param(i)%level2(1)==1.and. &
                     paramset(npset)%param(i)%level2(2)==-9999 .and. &
                     INDEX(trim(AVBLGRB2(IAVBL)),'1L')>0 ).or. &
                     (paramset(npset)%param(i)%level2(1)==-9999.and. &
                     paramset(npset)%param(i)%level2(2)==-9999 .and. &
                     INDEX(trim(AVBLGRB2(IAVBL)),'LLM')>0)

             if(LHYBPRES ) then
               IGET(IAVBL) = IFLD
               IDENT(IFLD) = IAVBL
               IAVBLFLD(IFLD)=I
               FOUND_FLD=.true.
               irec=irec+1
               LVLS(1,IFLD)=1
               LVLSXML(1,IFLD)=1
               IFLD=IFLD+1
               exit doavbl
             endif 
               
!            if(me==0.and.trim(paramset(npset)%param(i)%pname)=='CDLYR'.and. &
!               INDEX(trim(AVBLGRB2(IAVBL)),trim(paramset(npset)%param(i)%pname))/=0 ) &
!               print *,'CDLYR,LSTATS=',LSTATS,trim(paramset(npset)%param(i)%stats_proc), &
!               INDEX(trim(AVBLGRB2(IAVBL)),trim(paramset(npset)%param(i)%stats_proc))

             if (LSTATS .and. INDEX(trim(AVBLGRB2(IAVBL)),trim(paramset(npset)%param(i)%stats_proc))/=0 &
              .or. .NOT.LSTATS) then
!
              if(IGET(IAVBL)==-1) THEN
                IGET(IAVBL) = IFLD
                IDENT(IFLD) = IAVBL
                IAVBLFLD(IFLD)=I
                FOUND_FLD=.true.
                if((IAVBL==571.or.IAVBL==115).and.                  &
                  trim(paramset(npset)%param(i)%fixed_sfc2_type)/='') then
                 IGET(IAVBL) = -1
                endif
!for rh, tmp,pres,spf_h
                if((IAVBL==066.or.IAVBL==091.or.IAVBL==092.or.IAVBL==093.or.IAVBL==095  &
                    .or.IAVBL==096.or.IAVBL==210).and.                  &
                  trim(paramset(npset)%param(i)%fixed_sfc1_type)=='sigma_lvl') then
                 IGET(IAVBL) = -1
                endif
              endif
              if(size(paramset(npset)%param(i)%level)==0) then
                irec=irec+1
                LVLS(1,IFLD)=1
                LVLSXML(1,IFLD)=1
              else if(size(paramset(npset)%param(i)%level)==1.and.  &
               .not.LSIGMA) then
                irec=irec+1
                LVLS(1,IFLD)=1
                LVLSXML(1,IFLD)=1
              else
                irec=irec+size(paramset(npset)%param(i)%level)
                call getlvls(paramset(npset)%param(i),i,ifld,FOUND_FLD,kpv,pv)
                if(me==0.and.trim(paramset(npset)%param(i)%pname)=='RH'.and.   &
                 paramset(npset)%param(i)%fixed_sfc1_type=='sigma_lvl')  &   
                 print *,'for rh sigma,ifld=',ifld,'i=',i,'iget(66)=',iget(66),&
                  'IAVBLFLD(iget(66)=',IAVBLFLD(iget(66))
              endif
!                if(me==0.and.trim(paramset(npset)%param(i)%pname)=='MNTSF'.and.   &
!                 paramset(npset)%param(i)%fixed_sfc1_type=='isentropic_lvl')  &   
!                 print *,'ifld=',ifld,'iget(353)=',iget(353),'lvls(1,ifld)=',lvls(1,ifld), &
!                 'lvlsxml(1,ifld)=',LVLSXML(1,IFLD),'size(lvl)=',size(paramset(npset)%param(i)%level), &
!                 'lsigma=',lsigma,'ind1=',ind1,'ind2=',ind2,'ind3=',ind3
              if(FOUND_FLD)then
                IFLD=IFLD+1
                exit doavbl
              endif
!
             endif
            ENDIF
 20      ENDDO doavbl
         IF(ME.EQ.0.and..not.FOUND_FLD)THEN
           WRITE(6,*)'FIELD ',trim(paramset(npset)%param(i)%pname)//trim(        &
     &        paramset(npset)%param(i)%fixed_sfc1_type),' NOT AVAILABLE'
         ENDIF
         if(me==0) &
          print *,'in readxml,i=',i,'ifld=',ifld-1,'irec=',irec,  &
          trim(paramset(npset)%param(i)%pname),trim(paramset(npset)%param(i)%fixed_sfc1_type), &
          'lvl=',size(paramset(npset)%param(i)%level),'lvlsxml(1,ifld)=',LVLSXML(1,IFLD-1)
!
      ENDDO
!     
!     ALL DONE READING REQUESTED FIELDS FOR CURRENT OUTPUT GRID.
!     SET NFLD TO TOTAL NUMBER OF REQUESTED OUTPUT FIELDS THAT 
!     ARE AVAILABLE.
!
      NFLD = IFLD-1
      NRECOUT = IREC
      allocate(fld_info(NRECOUT))
      do i=1,nrecout
         fld_info(i)%ifld=0
        fld_info(i)%lvl=0
        fld_info(i)%lvl1=0
        fld_info(i)%lvl2=0
        fld_info(i)%ntrange=0
        fld_info(i)%tinvstat=0
      enddo
      write(0,*)'in readxml. nfld=',nfld,'nrecout=',nrecout
!
! skip creating ipv files if kth=0 and no isobaric fields are requested in ctl file      
      if(kth==0 .and. iget(013)<=0)go to 999
!     
!     ECHO OUTPUT FIELDS/LEVELS TO 6.
!
      IF(ME.EQ.0)THEN
        WRITE(6,*)'BELOW ARE FIELD/LEVEL/SMOOTHING ',       &
             'SPECIFICATIONS.,NFLD=',NFLD,'MXLVL=',MXLVL,'nrecout=',nrecout
      ENDIF
      DO 50 IFLD = 1,NFLD
        IF(ME.EQ.0)THEN
         i=IAVBLFLD(IFLD)
         write(0,*)'readxml,ifld=',ifld,'iget(',IDENT(ifld),')=',iget(ident(ifld)),'iavbl=',IAVBLFLD(iget(ident(ifld))),'postvar=',trim(paramset(npset)%param(i)%pname),  &
             trim(paramset(npset)%param(i)%fixed_sfc1_type),'lvls=',LVLS(:,ifld)
         if(size(paramset(npset)%param(i)%level)>0) then
           WRITE(0,*) paramset(npset)%param(i)%level
         endif
        ENDIF
 50   CONTINUE
!     
!     WE HAVE AN OUTPUT GRID AND THE FIELDS TO GENERATE ON IT.
!     SKIP OVER THE FOLLOWING EOF MESSAGE TO EXIT THIS ROUTINE.
!     
!     
!     END OF ROUTINE.
!     
 60   CONTINUE
 999  CONTINUE

!       write(0,*)'end of read_postcntrl_xml'
      RETURN
      END
