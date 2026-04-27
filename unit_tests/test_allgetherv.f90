! This is a test program for UPP.
!
! This program tests the ALLGETHERV() subroutine.
!
! Alyson Stahl, 4/2026
program test_allgetherv
    use mpi
    use ctlblk_mod, only: im, jm, num_procs, me, jsta, jend, ista, iend, mpi_comm_comp
    implicit none

    integer :: i, j, ierr, res
    real, allocatable :: GRID1(:,:), EXP_GRID1(:,:)

    interface
        subroutine ALLGETHERV(GRID1)
            use ctlblk_mod, only: im, jm
            REAL GRID1(im,jm)
        end subroutine ALLGETHERV
    end interface

    call MPI_INIT(ierr)
    mpi_comm_comp = MPI_COMM_WORLD
    call MPI_COMM_SIZE(mpi_comm_comp, num_procs, ierr)
    call MPI_COMM_RANK(mpi_comm_comp, me, ierr)

    im = 4
    jm = 4 * num_procs

    ista = 1
    iend = im
    jsta = me * (jm / num_procs) + 1
    jend = (me + 1) * (jm / num_procs)

    allocate(GRID1(im, jm))
    allocate(EXP_GRID1(im, jm))

    do i = ista, iend
        do j = jsta, jend
            GRID1(i,j) = real((j-1)*im + i)
            EXP_GRID1(i,j) = real((j-1)*im + i)
        end do
    end do

    call ALLGETHERV(GRID1)

    res = 0
    do i = ista, iend
        do j = jsta, jend
            if (GRID1(i,j) /= EXP_GRID1(i,j)) then
                print *, "ERROR: GRID1(", i, ",", j, ") = ", GRID1(i,j), &
                         " does not match expected value ", EXP_GRID1(i,j)
                res = 1
            end if
        end do
    end do

    call MPI_FINALIZE(ierr)

    deallocate(GRID1)
    deallocate(EXP_GRID1)

    if (res .ne. 0) stop 10

    print *, "SUCCESS!"

end program test_allgetherv
