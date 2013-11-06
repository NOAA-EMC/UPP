      SUBROUTINE SET_OUTFLDS(kth,kpv,pv)
!
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!                .      .    .     
! SUBPROGRAM:    READCNTRLgrb2_xml  READS POST xml CONTROL FILE
!   PRGRMMR: J. WANG         ORG: NCEP/EMC   DATE: 12-01-27       
!     
! ABSTRACT:
!     THIS ROUTINE READS THE CONTROL FILE IN XML FORMAT SPECIFYING
!     FIELD(S) TO POST, AND SAVE ALL THE FIELD INFORMATION IN
!     A DATATYPE array PSET 
!     
! PROGRAM HISTORY LOG:
!   01_27_2012  Jun Wang - INITIAL CODE
!     
! USAGE:    CALL READCNTRL_XML(kth,kpv,pv)
!   INPUT ARGUMENT LIST:
!     KTH 
!     KPV
!     PV
!
!   OUTPUT ARGUMENT LIST: 
!     NONE     - 
!     
!   OUTPUT FILES:
!     NONE
!     
!   SUBPROGRAMS CALLED:
!     UTILITIES:
!
!     LIBRARY:
!       COMMON   - RQSTFLDGRB2
!                  CTLBLK
!     
!   ATTRIBUTES:
!     LANGUAGE: FORTRAN
!     MACHINE : IBM      
!$$$  
!
!     
!     INCLUDE ETA GRID DIMENSIONS.  SET/DERIVE PARAMETERS.
!
       use xml_data_post_t,only: paramset,post_avblflds
       use grib2_module, only: num_pset,pset,nrecout,first_grbtbl,grib_info_init
       use lookup_mod,only: ITB,JTB,ITBQ,JTBQ
       use ctlblk_mod, only: npset, me, fld_info
!       use rqstfld_mod,only: num_post_afld,MXFLD,MXLVL,lvls,lvlsxml,iget
       use rqstfld_mod, only: mxfld, iget, ritehd, lvlsxml, datset, ident,& 
              iavblfld, nfld, lvls
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
       implicit none
!
!     DECLARE VARIABLES.
!     
      integer, intent(in) :: KTH,KPV
      real,intent(in) :: pv(kpv)
!
      integer L,IFLD,MFLD,IAVBL,IREC,I,J
      CHARACTER*50 AVBLGRB_NAME
      logical :: FOUND_FLD
!
!******************************************************************************
!     START READCNTRL_XML HERE.
!     
!      IF(ME.EQ.0)THEN
!        WRITE(6,*)'READCNTRL_XML:  POSTING FCST HR ',IFHR,' FROM ',         &
!             IHRST,'UTC ',SDAT(1),'-',SDAT(2),'-',SDAT(3),' RUN'
!      ENDIF
!     
!     INITIALIZE VARIABLES.
!        ARRAY IGET IS THE "GET FIELD" FLAG ARRAY.
!
      DO 100 IFLD=1,MXFLD
        IGET(IFLD)=-1
 100  CONTINUE
!
!     SET FLAG TO OPEN NEW OUTPUT FILE
!
      LVLS=0
      RITEHD = .TRUE.
      DO J=1,size(LVLSXML,2)
      DO I=1,size(LVLSXML,1)
        LVLSXML(I,J)=0
      ENDDO
      ENDDO
!
      pset=paramset(npset)
      datset=pset%datset
      print *,'in readxml, num_pset=',num_pset,'datset=',trim(pset%datset),'npset=',npset
! 
!     NOW READ WHICH FIELDS ON 
!     WHICH LEVELS TO INTERPOLATE TO THE OUTPUT GRID.  THE
!     CHARACTER STRING "DONE" MARKS THE END OF THE OUTPUT
!     FIELD SPECIFICATIONS.
!
      call grib_info_init()
      MFLD = size(pset%param)
      print *,'start reading control file,MFLD=',MFLD,'datset=',datset,MXFLD
      if(size(post_avblflds%param)<=0) then
         print *,'WRONG: post available fields not ready!!!'
         return
      endif
!
      IFLD=0
      irec=0
      DO I=1, MFLD
         
!     SEE IF REQUESTED FIELD IS AVAILABLE.  IF NOT, 
!     WRITE MESSAGE TO 6 AND DECREMENT FIELD 
!     COUNTER BY ONE.  THEN READ NEXT REQUESTED FIELD.
!     
!     GET POST AVAILBLE FIELD INDEX NUMBER FOR EACH OUTPUT FIELDS IN PSET
!
         FOUND_FLD=.false.
         print *,'cntfile,i=',i,'fld shortname=',trim(pset%param(i)%shortname)
         doavbl:   DO 20 J = 1,size(post_avblflds%param)
             
            if(trim(pset%param(i)%shortname)==trim(post_avblflds%param(j)%shortname)) then
!
              IFLD=IFLD+1
              IAVBL=post_avblflds%param(j)%post_avblfldidx
              IGET(IAVBL)=IFLD
              IDENT(IFLD)=IAVBL
              IAVBLFLD(IFLD)=I
              FOUND_FLD=.true.
              call fill_psetfld(pset%param(i),post_avblflds%param(j))
              call set_lvlsxml(pset%param(i),ifld,irec,kpv,pv)
              exit doavbl
!
            endif

 20      ENDDO doavbl
!
         IF(ME.EQ.0.and..not.FOUND_FLD)THEN
           WRITE(0,*)'FIELD ',trim(pset%param(i)%pname)//trim(        &
     &        pset%param(i)%fixed_sfc1_type),' NOT AVAILABLE'
         ENDIF
!         if(me==0) &
!          print *,'in readxml,i=',i,'ifld=',ifld,'irec=',irec,  &
!          trim(pset%param(i)%pname),trim(pset%param(i)%fixed_sfc1_type), &
!          'lvl=',size(pset%param(i)%level),'lvlsxml(1,ifld)=',LVLSXML(1,IFLD), &
!          'ident(ifld)=',ident(ifld),'iget(ident(ifld))=',iget(ident(ifld))
!
      ENDDO
!     
!     ALL DONE FOUNDING REQUESTED FIELDS FOR current OUTPUT GRID.
!     SET NFLD TO TOTAL NUMBER OF REQUESTED OUTPUT FIELDS THAT 
!     ARE AVAILABLE., SET NRECOUT to total number of OUTPUT records 
!     NOTE: here NFLD i s total number of fields found in post_avblfld_table,
!           while nrecoutis the total number of grib messages that go 
!           into the output file. One fieldmay contain many different levels,
!           which each different level will be counted as one record
!
      NFLD = IFLD
      NRECOUT = IREC
      allocate(fld_info(NRECOUT+100))
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
!      IF(ME.EQ.0)THEN
!        WRITE(6,*)'BELOW ARE FIELD/LEVEL/SMOOTHING ',       &
!             'SPECIFICATIONS.,NFLD=',NFLD,'MXLVL=',MXLVL,'nrecout=',nrecout
!      ENDIF
!      DO 50 IFLD = 1,NFLD
!        IF(ME.EQ.0)THEN
!         i=IAVBLFLD(IFLD)
!         write(0,*)'readxml,ifld=',ifld,'iget(',IDENT(ifld),')=',iget(ident(ifld)),'iavbl=',IAVBLFLD(iget(ident(ifld))),'postvar=',trim(pset%param(i)%pname),  &
!             trim(pset%param(i)%fixed_sfc1_type),'lvls=',LVLS(:,ifld)
!         if(size(pset%param(i)%level)>0) then
!           WRITE(0,*) pset%param(i)%level
!         endif
!        ENDIF
! 50   CONTINUE
!     
!     END OF ROUTINE.
!     
 999  CONTINUE

       print *,'end of read_postcntrl_xml'
      RETURN
      END
