!> @file
!> @brief MACHINE_POST defines machine-dependent constants
      MODULE MACHINE_POST

      IMPLICIT NONE
      SAVE
      
      integer kind_io4 &    !< 4-byte I/O variables ?
      ,kind_io8 &           !< 8-byte I/O variables ?
      ,kind_phys &          !< Physics variables ?
      ,kind_rad &           !< Radiation variables ?
      ,kint_mpi             !< MPI variables ?
      
      parameter (kind_rad = selected_real_kind(13,60)) ! the '60' maps to 64-bit real
      parameter (kind_phys = selected_real_kind(13,60)) ! the '60' maps to 64-bit real
      parameter (kind_io4 = 4)
      parameter (kind_io8 = 8)
      parameter (kint_mpi = 4)

      END MODULE MACHINE_POST
