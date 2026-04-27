! This is a test program for UPP.
!
! This program tests ASSIGNNEMSIOVAR() subroutine.
!
! Alyson Stahl, 2/2026
program test_assignnemsiovar
    use ctlblk_mod, only: me
    implicit none

    interface
        subroutine ASSIGNNEMSIOVAR(IM, JSTA, JEND, JSTA, JSTA_2L, JEND_2U, &
                                    L, NREC, FLDSIZE, SPVAL, TMP, RECNAME, &
                                    RECLEVTYP, RECLEV, VARNAME, VCOORDNAME, BUF)
            integer, intent(in) :: IM, JSTA, JEND, JSTA_2L, JEND_2U, L, NREC, FLDSIZE
            integer, intent(in) :: RECLEV(NREC)
            real, intent(in) :: SPVAL, TMP(FLDSIZE * NREC)
            character(*), intent(in) :: RECNAME(NREC), RECLEVTYP(NREC)
            character(*), intent(in) :: VARNAME, VCOORDNAME
            real, intent(out) :: BUF(IM, JSTA_2L:JEND_2U)
        end subroutine ASSIGNNEMSIOVAR
    end interface


    print *, "SUCCESS!"

end program test_assignnemsiovar
