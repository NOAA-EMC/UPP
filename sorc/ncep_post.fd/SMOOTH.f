!> @file
!> @brief smooth() smooths a meteorological field using Shapiro smoother.
!>
!> @author Stan Benjamin FSL/PROFS @date 1990-06-15
!>
!> @note Reference: Shapiro, 1970: "Smoothing, filtering, and
!> boundary effects", REV. GEOPHYS. SP. PHYS., 359-387.
!> This filter is of the type 
!> @code
!>       Z(I) = (1-S)Z(I) + S(Z(I+1)+Z(I-1))/2
!> @endcode
!> For a filter which is supposed to damp 2DX waves completely
!> but leave 4DX and longer with little damping,
!> it should be run with 2 passes using SMTH (or s) of 0.5
!> and -0.5.
!>
!> @param[in] FIELD Real array FIELD(IX,IY) Meteorological field.
!> @param[in] HOLD Real array HOLD(IX,2) Holding the value for field.
!> @param[in] IX Integer X Coordinates of field.
!> @param[in] IY Integer Y Coordinates of field.
!> @param[in] SMTH Real.
!> @return FIELD Real array FIELD(IX,IY) Smoothed meteorological field.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1990-06-15 | S. Benjamin | Initial
!> 2014-03-03 | S. Moorthi  | Threading and slight cleanup
!>
!> @author Stan Benjamin FSL/PROFS @date 1990-06-15
!**********************************************************************
!**********************************************************************
!> smooth() smooths a meteorological field using Shapiro smoother.
!>
!> @param[inout] FIELD Real array FIELD(IX,IY) Meteorological field.
!> @note FIELD is converted to *smoothed* meteorological field in this function.
!> @param[in] HOLD Real array HOLD(IX,2) Holding the value for field.
!> @param[in] IX Integer X Coordinates of field.
!> @param[in] IY Integer Y Coordinates of field.
!> @param[in] SMTH Real.
!>
      SUBROUTINE SMOOTH (FIELD,HOLD,IX,IY,SMTH)

!**********************************************************************
!**********************************************************************

      implicit none

      integer :: i1, i2, j, it, i, ix, iy
      real    :: smth1, smth, smth2, smth3, smth4, smth5 
      real    :: sum1, sum2
      REAL      FIELD(IX,IY), HOLD (IX,2)
      SMTH1 = 0.25 * SMTH * SMTH
      SMTH2 = 0.5  * SMTH * (1.-SMTH)
      SMTH3 = (1.-SMTH) * (1.-SMTH)
      SMTH4 = (1.-SMTH)
      SMTH5 = 0.5 * SMTH
      I1 = 2
      I2 = 1
       DO J=2,IY-1
         IT = I1
         I1 = I2
         I2 = IT
!$omp parallel do private(i,sum1,sum2)
         DO I = 2,IX-1
           SUM1 = FIELD (I-1,J+1) + FIELD (I-1,J-1)                 &
                + FIELD (I+1,J+1) + FIELD (I+1,J-1)
           SUM2 = FIELD (I  ,J+1) + FIELD (I+1,J  )                 &
                + FIELD (I  ,J-1) + FIELD (I-1,J  )
           HOLD(I,I1) = SMTH1*SUM1 + SMTH2*SUM2 + SMTH3*FIELD(I,J)
         ENDDO
         IF (J > 2) then
!$omp parallel do private(i)
           DO I=2,IX-1
             FIELD(I,J-1) = HOLD(I,I2)
           ENDDO
         endif
       ENDDO

!$omp parallel do private(i)
       DO I = 2,IX-1
         FIELD (I,IY-1) = HOLD(I,I1)
       ENDDO

       DO I = 2,IX-1
         FIELD(I,1)  = SMTH4 * FIELD(I,1)                             &
                     + SMTH5 * (FIELD(I-1,1) + FIELD(I+1,1))
         FIELD(I,IY) = SMTH4 * FIELD(I,IY)                            &
                     + SMTH5 * (FIELD(I-1,IY) + FIELD(I+1,IY))
       ENDDO

       DO J = 2,IY-1
         FIELD(1,J)  = SMTH4 * FIELD(1,J)                             &
                     + SMTH5 * (FIELD(1,J-1) + FIELD(1,J+1))
         FIELD(IX,J) = SMTH4 * FIELD(IX,J)                            &
                     + SMTH5 * (FIELD(IX,J-1) + FIELD(IX,J+1))
       ENDDO

      RETURN
      END
!> @brief smoothc() smooths a meteorological field using Shapiro smoother.
!>
!> @author Stan Benjamin FSL/PROFS @date 1990-06-15
!> @note Reference: Shapiro, 1970: "Smoothing, filtering, and
!> boundary effects", REV. GEOPHYS. SP. PHYS., 359-387.
!> This filter is of the type 
!> @code
!>       Z(I) = (1-S)Z(I) + S(Z(I+1)+Z(I-1))/2
!> @endcode
!> For a filter which is supposed to damp 2DX waves completely
!> but leave 4DX and longer with little damping,
!> it should be run with 2 passes using SMTH (or s) of 0.5
!> and -0.5.
!>
!> @param[inout] FIELD Real array FIELD(IX,IY) Meteorological field.
!> @note FIELD is converted to *smoothed* meteorological field in this function.
!> @param[in] HOLD Real array HOLD(IX,2) Holding the value for field.
!> @param[in] IX Integer X Coordinates of field.
!> @param[in] IY Integer Y Coordinates of field.
!> @param[in] SMTH Real.
!>
!> ### Program history log:
!> Date | Programmer | Comments
!> -----|------------|---------
!> 1985-12-09 | S. Benjamin | Original version os smooth
!> 2014-03-03 | S. Moorthi  | Threading and slight cleanup
!> 2016-08-08 | S. Moorthi  | Modify for cyclic domain 
!>
!> @author Stan Benjamin FSL/PROFS @date 1990-06-15
!**********************************************************************
!**********************************************************************

      SUBROUTINE SMOOTHC (FIELD,HOLD,IX,IY,SMTH)

!**********************************************************************
!**********************************************************************

      implicit none

      integer :: i1, i2, j, it, i, ix, iy, im1, ip1
      real    :: smth1, smth, smth2, smth3, smth4, smth5 
      real    :: sum1, sum2
      REAL      FIELD(IX,IY), HOLD (IX,2)
      integer :: iw(ix), ie(ix)
!
      SMTH1 = 0.25 * SMTH * SMTH
      SMTH2 = 0.5  * SMTH * (1.-SMTH)
      SMTH3 = (1.-SMTH) * (1.-SMTH)
      SMTH4 = (1.-SMTH)
      SMTH5 = 0.5 * SMTH
!
      do i=2,ix-1
        ie(i) = i + 1
        iw(i) = i - 1
      enddo
      ie(ix) = 1
      iw(1)  =  ix
!
      I1 = 2
      I2 = 1
       DO J=2,IY-1
         IT = I1
         I1 = I2
         I2 = IT
!$omp parallel do private(i,sum1,sum2,ip1,im1)
         DO I = 1,IX
           ip1 = ie(i)
           im1 = iw(i)
           SUM1 = FIELD (Im1,J+1) + FIELD (Im1,J-1)                 &
                + FIELD (Ip1,J+1) + FIELD (Ip1,J-1)
           SUM2 = FIELD (I  ,J+1) + FIELD (Ip1,J  )                 &
                + FIELD (I  ,J-1) + FIELD (Im1,J  )
           HOLD(I,I1) = SMTH1*SUM1 + SMTH2*SUM2 + SMTH3*FIELD(I,J)
         ENDDO
         IF (J > 2) then
!$omp parallel do private(i)
           DO I=1,IX
             FIELD(I,J-1) = HOLD(I,I2)
           ENDDO
         endif
       ENDDO

!$omp parallel do private(i)
       DO I = 1,IX
         FIELD (I,IY-1) = HOLD(I,I1)
       ENDDO

       DO I = 1,IX
         ip1 = ie(i)
         im1 = iw(i)
         FIELD(I,1)  = SMTH4 * FIELD(I,1)                             &
                     + SMTH5 * (FIELD(Im1,1) + FIELD(Ip1,1))
         FIELD(I,IY) = SMTH4 * FIELD(I,IY)                            &
                     + SMTH5 * (FIELD(Im1,IY) + FIELD(Ip1,IY))
       ENDDO


      RETURN
      END
