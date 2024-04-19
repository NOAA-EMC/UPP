!> @file
!> @brief module: MACHINE_POST defines machine-dependent constants
      MODULE MACHINE_POST

      IMPLICIT NONE
      SAVE
      
      Integers 
      kind_io4(:,:) &       !< Integer parameter specifying the kind value for 4-byte integers
      ,kind_io8(:,:) &      !< Integer parameter specifying the kind value for 8-byte integers
      ,kind_phys(:,:) &     !< Integer parameter specifying the kind value for real numbers used in physics
      ,kind_rad(:,:) &      !< Integer parameter specifying the kind value for real numbers used in radiation calculations
      ,kint_mpi(:,:) &      !< Integer parameter specifying the kind value for MPI integer types
      
      parameter (kind_rad = selected_real_kind(13,60)) ! the '60' maps to 64-bit real
      parameter (kind_phys = selected_real_kind(13,60)) ! the '60' maps to 64-bit real
      parameter (kind_io4 = 4)
      parameter (kind_io8 = 8)
      parameter (kint_mpi = 4)

      END MODULE MACHINE_POST
