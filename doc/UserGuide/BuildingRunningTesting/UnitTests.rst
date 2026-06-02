.. _unit-tests:

**********
Unit Tests
**********

Overview
========

Unit tests help verify that individual UPP routines behave as expected and
continue to produce correct results as the code evolves. They are an
important safeguard against regressions, especially when existing routines
are modified without an intended change in output.

New code should be written so that it can be fully tested. Ideally,
developers should be able to write unit tests that achieve 100% line and
branch coverage for the routine being tested.

By default, unit tests are compiled automatically with CMake when UPP is
built. The root CMake configuration includes the ``unit_tests``
directory, which contains its own ``CMakeLists.txt`` file. No additional
steps are required to build unit tests locally or in GitHub CI.

What Makes a Good Unit Test in UPP?
===================================

A good UPP unit test verifies the expected behavior of the full function
or subroutine being tested. It should exercise the routine thoroughly
enough to confirm that the implementation is correct, stable, and
protected against future regressions.

At a minimum, a unit test should verify:

* Every ``if``/``else`` branch
* Every error condition
* Important edge and corner cases
* Every output variable, including variables updated by the routine
* Expected error codes, when routines return them

Good unit tests should also be reliable and self-contained. They should:

* Initialize all required state explicitly
* Use deterministic inputs
* Define expected values independently from the implementation being
  tested
* Clean up allocated resources before exiting

Unit tests should return a non-zero stop code when a failure occurs. They
should also print descriptive error messages that make the failure easy to
diagnose. Useful error messages identify:

* The failing output variable
* The expected value and actual value
* The precise array index where the failure occurred, when applicable

Each test case should include a short comment that clearly describes its
purpose, such as the branch, error condition, or edge case being tested.

UPP code should also be written with testability in mind. When possible,
avoid adding:

* ``STOP`` statements, which Fortran cannot intercept, making it
  difficult to verify in a unit test that the ``STOP`` statement was
  called as expected
* Branches that cannot be reached with any input

Example: ``test_calicing.f90``
------------------------------

The ``test_calicing.f90`` unit test is a useful example for new tests. It
tests a small subroutine with minimal setup, but includes enough
``if``/``else`` branches, boundary conditions, and missing-value cases to
demonstrate how to write a thorough unit test.

.. code-block:: fortran

   ! This is a test program for UPP.
   !
   ! This program tests the CALICING() subroutine.
   !
   ! Alyson Stahl, 4/2026
   program test_calicing
       use ctlblk_mod, only: jsta, jend, spval, ista, iend
       implicit none

       real, parameter :: tol = 1.0e-8
       integer, parameter :: npts = 11
       integer :: i, res
       real :: T1(1, npts), RH(1, npts), OMGA(1, npts)
       real :: ICING(1, npts), EXP_ICING(1, npts)

       interface
           subroutine CALICING(T1,RH,OMGA, ICING)
               use ctlblk_mod, only: jsta, jend, ista, iend
               real, dimension(ista:iend,jsta:jend), intent(in) :: T1,RH,OMGA
               real, dimension(ista:iend,jsta:jend), intent(inout) :: ICING
           end subroutine CALICING
       end interface

       ! Grid parameters
       ista = 1
       iend = 1
       jsta = 1
       jend = npts
       spval = 9.9e10

       ! Test Case 1: OMGA < 0 & 251 < T1 < 273 & RH > 70 (expect ICING = 1)
       T1 = 260.0
       RH = 80.0
       OMGA = -0.1
       EXP_ICING(1,1) = 1.0

       ! Test Case 2: OMGA > 0 (expect ICING = 0)
       OMGA(1,2) = 0.1
       EXP_ICING(1,2) = 0.0

       ! Test Case 3: T1 < 251 (expect ICING = 0)
       T1(1,3) = 250.0
       EXP_ICING(1,3) = 0.0

       ! Test Case 4: T1 > 273 (expect ICING = 0)
       T1(1,4) = 274.0
       EXP_ICING(1,4) = 0.0

       ! Test Case 5: RH < 70 (expect ICING = 0)
       RH(1,5) = 60.0
       EXP_ICING(1,5) = 0.0

       ! Test Case 6: OMGA < 0 & T1 == 251 & RH > 70 (expect ICING = 1)
       T1(1,6) = 251.0
       EXP_ICING(1,6) = 1.0

       ! Test Case 7: OMGA < 0 & T1 == 273 & RH > 70 (expect ICING = 1)
       T1(1,7) = 273.0
       EXP_ICING(1,7) = 1.0

       ! Test Case 8: OMGA < 0 & 251 < T1 < 273 & RH == 70 (expect ICING = 1)
       RH(1,8) = 70.0
       EXP_ICING(1,8) = 1.0

       ! Test Case 9: OMGA > spval (expect ICING = spval)
       OMGA(1,9) = spval
       EXP_ICING(1,9) = spval

       ! Test Case 10: T1 > spval (expect ICING = spval)
       T1(1,10) = spval
       EXP_ICING(1,10) = spval

       ! Test Case 11: RH > spval (expect ICING = spval)
       RH(1,11) = spval
       EXP_ICING(1,11) = spval

       call CALICING(T1, RH, OMGA, ICING)

       res = 0
       do i = 1, npts
           if (abs(ICING(1,i) - EXP_ICING(1,i)) > tol) then
               print *, "Test Case ", i, " failed: ICING = ", ICING(1,i), &
                        " but expected ", EXP_ICING(1,i)
               res = 1
           end if
       end do

       if (res .ne. 0) stop 10

       print *, "SUCCESS!"
   end program test_calicing

This test is a good model because it demonstrates several practices that
are useful throughout UPP unit testing:

* It initializes required global state variables from ``ctlblk_mod`` so
  the subroutine can be tested outside of the full UPP workflow
* It creates an explicit interface for ``CALICING()``, which is needed
  because the subroutine is not part of a module
* It documents each test case with a comment that explains the condition
  being tested and the expected result
* It uses array inputs to evaluate multiple test cases with a single
  subroutine call
* It checks the computed value against the expected value for each test
  case using a defined floating-point tolerance
* It prints descriptive failure messages that identify the failing test
  case, the actual value, and the expected value
* It returns a non-zero stop code when any test case fails

Example: ``test_calgustconv.f90``
---------------------------------

The ``test_calgustconv.f90`` unit test is a useful template for routines
that depend on UPP global data arrays. It shows how to allocate,
initialize, and clean up arrays that are normally managed by the full UPP
workflow.

The following excerpt shows the global data setup:

.. code-block:: fortran

   program test_calgustconv
       use vrbls2d , only: u10, v10, ustar
       use ctlblk_mod, only: ista, iend, jsta, jend, ista_2l, iend_2u, &
                             jsta_2l, jend_2u, spval
       implicit none

       real, parameter :: tol = 1.0e-8
       integer, parameter :: npts = 5
       integer :: j, res
       real :: SPEED850(1, 1:npts), SPEED950(1, 1:npts)
       real :: GUSTCONV(1, 1:npts), EXP_GUSTCONV(1, 1:npts)

       ! Grid dimensions
       ista = 1
       iend = 1
       jsta = 1
       jend = npts
       ista_2l = 1
       iend_2u = 1
       jsta_2l = 1
       jend_2u = npts
       spval = 9.9e10

       allocate(u10(ista_2l:iend_2u,jsta_2l:jend_2u))
       allocate(v10(ista_2l:iend_2u,jsta_2l:jend_2u))
       allocate(ustar(ista_2l:iend_2u,jsta_2l:jend_2u))

       u10 = 12.0
       v10 = 16.0
       ustar = 0.7

       SPEED950 = 26.0
       SPEED850 = 34.0

       call CALGUSTCONV(SPEED850,SPEED950,GUSTCONV)

       deallocate(u10, v10, ustar)

This test is a good model because it:

* Imports the required global arrays from ``vrbls2d``
* Sets the related grid dimensions from ``ctlblk_mod`` before allocating
  the arrays
* Allocates and initializes the global arrays before calling
  ``CALGUSTCONV()``
* Deallocates the arrays after the subroutine call
* Covers typical input, negative wind-speed differences that are treated
  as zero, equal wind speeds, and ``spval`` handling

Example: ``test_select_channels.f90``
-------------------------------------

The ``test_select_channels.f90`` unit test is a useful example for tests
that call the same subroutine multiple times. Each call uses different
input values and checks the result before moving to the next case.

The following excerpt shows the repeated subroutine calls:

.. code-block:: fortran

   ! Test Case 1: All L = 0
   L = 0
   IGOT = 99

   call SELECT_CHANNELS_L(CHANNELINFO, NCHANNELS, CHANNELS, L, IGOT)

   if (IGOT .ne. 0) then
       print *, 'Test Case 1 Failed: Expected IGOT = 0, got ', IGOT
       res = 1
   end if

   do i = 1, NCHANNELS
       CHANNELINFO%Process_Channel(i) = .true.
   end do

   ! Test case 2: Mixed L values
   L    = (/ 1, 0, 1, 0, 1 /)
   IGOT = 5

   call SELECT_CHANNELS_L(CHANNELINFO, NCHANNELS, CHANNELS, L, IGOT)

   if (IGOT .ne. 5) then
       print *, 'Test Case 2 Failed: Expected IGOT = 5, got ', IGOT
       res = 1
   end if

This test is a good model because it:

* Calls ``SELECT_CHANNELS_L()`` multiple times with different inputs
* Resets ``CHANNELINFO%Process_Channel`` between test cases
* Checks both single-value output, such as ``IGOT``, and array values,
  such as ``CHANNELINFO%Process_Channel``
* Prints descriptive messages when a case fails