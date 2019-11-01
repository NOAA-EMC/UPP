      SUBROUTINE CALVESSEL(ICEG)
! Algorithm for calculating ice growth rate
      use vrbls2d, only: sst, u10h, v10h, tshltr
      use masks, only: sm, sice
      use ctlblk_mod, only: jsta, jend, im, spval
!-------------------------------------------
      implicit none
      integer I, J
      real TSFC_C,TSHLTR_C,SST_C
      real, parameter :: C2K=273.15
      real, dimension(im,jsta:jend) :: pr, spd10
      real,intent(out) ::  ICEG(im,jsta:jend)

!      allocate (thsfc(im,jsta:jend),tsfc(im,jsta:jend))

      DO J=JSTA,JEND
        DO I=1,IM
!   CALCULATE SPEED
          SPD10(i,j)=SQRT(U10H(I,J)**2+V10H(I,J)**2)
          if (SPD10(i,j).gt.50) then
            iceg(i,j)=0.
            CYCLE
          endif

! Reverse of land mask use le instead of ge from original code
!!  MASK CHECK
          if((sice(i,j).ge.0.5).or.(sm(i,j).le.0.5)) then
            ICEG(i,j)=0.
            CYCLE
          endif

!!! CHANGE TEMP to FROM K to C
!!! TEMPERATURE CHECK
          SST_C=SST(I,J)-C2K  
          TSHLTR_C=TSHLTR(I,J)-C2K
          if((SST_C.lt.-1.7).OR. &
             (SST_C.gt.12.0)) then
            ICEG(I,j)=0.
            CYCLE
          endif

          if((TSHLTR_C.gt.0.).OR. &
             (TSHLTR_C.lt.-40.)) then
            ICEG(I,j)=0.
            CYCLE
          endif

!  CALCULATE ICE GROWTH
          PR(i,j)=SPD10(i,j)*(-1.7-TSHLTR_C)/(1.+.4*(SST_C+1.7))
          ICEG(i,j)=(2.73E-02)*PR(i,j)+(2.91E-04)*PR(i,j)*PR(i,j) &
                   +(1.84E-06)*PR(i,j)**3

!! ICE GROWTH CHECK
          if (ICEG(i,j).LT.0.) THEN
            ICEG(I,J)=0.
          else
! Convert to m/s from cm/hr
            ICEG(i,j)=(1./3.6E+05)*ICEG(I,J)
          endif

        ENDDO
      ENDDO

      END
