!> @file
!> @brief get_bits() computes number of bits and round field.
!>
!> The number of bits requited to pack a given field
!> at a particular decimal scaling is computed using the field range.
!> The field is rounded off to the decimal scaling for packing.
!> The minimum and maximum rounded field values are also returned.
!> Grib bitmap masking for valid data is optionally used.
!>
!> @param[in] IBM Integer bitmap flag (=0 for no bitmap).
!> @param[in] SGDS Maximum significant digits to keep.
!><pre>
!> (E.G. SGDS=3.0 keeps 3 significant digits)
!>  or binary precision if <0
!> (E.G. SGDS=-2.0 keeps field to nearest 1/4
!>            -3.0 keeps field to nearest 1/8
!>            2**SGDS precision)
!></pre>
!> @param[in] LEN Integer length of the field and bitmap.
!> @param[in] MG Integet (LEN) bitmap if IBM=1 (0 to skip, 1 tp keep).
!> @param[in] G Real (LEN) field.
!> @param[out] ISCALE Integer decimal scaling.
!> @param[out] GROUND Real (LEN) field rounded to decimal scaling.
!> @param[out] GMIN Real minimum valid rounded field value.
!> @param[out] GMAX Real maximum valid rounded field value.
!> @param[out] NBIT Integer number of bits to pack.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-10-31 | Iredell | Initial
!> 1995-04-14 | Baldwin | Modify following Keith Brill's code to use sig digits to compute DEC scale
!>
!> @author Iredell W/NP23 @date 1992-10-31
!-----------------------------------------------------------------------
!> get_bits() computes number of bits and round field.
!> 
!> @param[in] IBM integer bitmap flag (=0 for no bitmap).
!> @param[in] SGDS real Maximum significant digits to keep.
!> @param[in] LEN integer Length of the field and bitmap.
!> @param[in] MG integer (LEN) Bitmap if IBM=1 (0 to skip, 1 tp keep).
!> @param[in] G real (LEN) Field.
!> @param[inout] ISCALE integer Decimal scaling.
!> @param[inout] GROUND real (LEN) Field rounded to decimal scaling.
!> @param[out] GMIN real Minimum valid rounded field value.
!> @param[out] GMAX real Maximum valid rounded field value.
!> @param[inout] NBIT integer Number of bits to pack.
!>
      SUBROUTINE GET_BITS(IBM,SGDS,LEN,MG,G,ISCALE,GROUND,           &
                          GMIN,GMAX,NBIT)

!
      implicit none
!
      REal,DIMENSION(LEN),intent(in):: G
      real,DIMENSION(LEN),intent(inout) ::  GROUND
      integer,DIMENSION(LEN),intent(in):: MG
      integer,intent(in) :: IBM,LEN
      integer,intent(inout) :: ISCALE,NBIT
      real,intent(out) :: GMAX,GMIN
      integer I1,I,IRETT
      real SGDS
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  DETERMINE EXTREMES WHERE BITMAP IS ON
!
      IF(IBM==0) THEN
        GMAX=G(1)
        GMIN=G(1)
        DO I=2,LEN
          GMAX=MAX(GMAX,G(I))
          GMIN=MIN(GMIN,G(I))
        ENDDO
      ELSE
        I1=0
        DO I=1,LEN
          IF(MG(I)/=0.AND.I1==0) I1=I
        ENDDO
        IF(I1>0.AND.I1<=LEN) THEN
          GMAX=G(I1)
          GMIN=G(I1)
          DO I=I1+1,LEN
            IF(MG(I)/=0) THEN
              GMAX=MAX(GMAX,G(I))
              GMIN=MIN(GMIN,G(I))
            ENDIF
          ENDDO
        ELSE
          GMAX=0.
          GMIN=0.
        ENDIF
      ENDIF
!
!
!
      CALL FNDBIT  ( GMIN, GMAX, SGDS, NBIT, ISCALE, IRETT)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      RETURN
      END
!> fndbit() computes the number of packing bits given the
!> maximum number of significant digits to preserve or the binary
!> precision to store the data.  The binary precision is given as a
!> negative integer, ISCALE will always be zero in this case.
!>
!> The binary precision translates as follows:
!> <pre>
!>     -1  =>  store data to nearest 1/2
!>     -2  =>  store data to nearest 1/4
!>     -3  =>  store data to nearest 1/8
!> </pre>
!>
!> Note that a fractional number of significant digits is allowed.
!>
!> @param[in] RMIN Real Minimum value.
!> @param[in] RMAX Real Maximum value.
!> @param[in] RDB Real Maximum # of significant digits OR binary precision if < 0.
!> @param[out] NMBTS Integer Number of bits for packing.
!> @param[out] ISCALE Integer Power of 10 scaling to use.
!> @param[out] IRET Integer Return code. 0 = normal return.
!>
!> ### Program History Log
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1992-06-?? | K Brill   | Initial
!> 1995-12-?? | K Brill   | Added binary precision
!> 1996-10-?? | M Baldwin | Added fix for negative nmbts
!>
!> @author K. Brill NMC @date 1992-06-??
    SUBROUTINE FNDBIT  ( rmin, rmax, rdb, nmbts, iscale, iret )
    implicit none
!
    integer,intent(inout) ::  iscale,nmbts
    real,intent(inout)    ::  rmin,rmax,rdb
    real                  ::  range,rr,rng2,po,rln2
    integer               ::  iret,icnt,ipo,le,ibin
!    
	DATA		rln2/0.69314718/
!-----------------------------------------------------------------------
	iret = 0
	icnt = 0
	iscale = 0
	range = rmax - rmin
	IF ( range <= 0.00 ) THEN
	    nmbts = 8
	    RETURN
	END IF
!*
	IF ( rdb == 0.0 ) THEN
	    nmbts = 8
	    RETURN
	ELSE IF ( rdb > 0.0 ) THEN
	    ipo = INT (ALOG10 ( range ))
	    IF ( range < 1.00 ) ipo = ipo - 1
	    po = float(ipo) - rdb + 1.
	    iscale = - INT ( po )
	    rr = range * 10. ** ( -po )
	    nmbts = INT ( ALOG ( rr ) / rln2 ) + 1
	ELSE
	    ibin = NINT ( -rdb )
	    rng2 = range * 2. ** ibin
	    nmbts = INT ( ALOG ( rng2 ) / rln2 ) + 1
	END IF
!*
        IF(NMBTS<=0) THEN
          NMBTS=0
          IF(ABS(RMIN)>=1.) THEN
            ISCALE=-INT(ALOG10(ABS(RMIN)))
          ELSE IF (ABS(RMIN)<1.0.AND.ABS(RMIN)>0.0) THEN
            ISCALE=-INT(ALOG10(ABS(RMIN)))+1
          ELSE
            ISCALE=0
          ENDIF
        ENDIF
	RETURN
	END
