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
developers should be able to write unit tests that achieve 100% line
coverage and branch coverage for the routine being tested.

By default, unit tests are compiled automatically with :term:`CMake` when
UPP is built. The root CMake configuration includes the ``unit_tests``
directory, which contains its own ``CMakeLists.txt`` file. No additional
steps are required to build unit tests locally or in GitHub :term:`CI`.

.. _good-unit-test-upp:

What Makes a Good Unit Test in UPP?
===================================

A good UPP unit test verifies the expected behavior of a complete
function or subroutine. It should exercise the routine thoroughly enough
to confirm that the implementation is correct, stable, and protected
against future regressions.

Coverage
   At a minimum, a unit test should verify every ``if``/``else`` branch,
   every error condition, important edge and corner cases, every output
   variable, and any expected error codes returned by the routine.

Reliability
   Unit tests should be reliable and self-contained. They should
   initialize all required state explicitly, use deterministic inputs,
   define expected values independently of the implementation being
   tested, and clean up allocated resources before exiting.

Failure messages
   Unit tests should return a non-zero stop code when a failure occurs.
   They should also print descriptive error messages that identify the
   failing output variable, the expected value, the actual value, and the
   array index where the failure occurred, when applicable.

Test intent
   Each test case should include a short comment that clearly describes
   its purpose, such as the branch, error condition, or edge case being
   tested.

Testability
   UPP code should be written with testability in mind. When possible,
   avoid adding ``STOP`` statements, which Fortran cannot intercept, and
   avoid adding branches that cannot be reached with any input.

.. _example-test-calicing:

Example: ``test_calicing.f90``
------------------------------

The ``test_calicing.f90`` unit test is a useful example for new tests. It
tests a small subroutine with minimal setup and includes enough
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

       real, parameter :: tol = 1.0e-6
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

This example demonstrates several practices that are useful throughout UPP
unit testing:

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

.. _example-test-calgustconv:

Example: ``test_calgustconv.f90``
---------------------------------

The ``test_calgustconv.f90`` unit test is a useful template for routines
that depend on UPP global data arrays. It shows how to allocate,
initialize, and clean up arrays that the full UPP workflow normally
manages.

The following excerpt shows the global data setup:

.. code-block:: fortran

   program test_calgustconv
       use vrbls2d , only: u10, v10, ustar
       use ctlblk_mod, only: ista, iend, jsta, jend, ista_2l, iend_2u, &
                             jsta_2l, jend_2u, spval
       implicit none

       real, parameter :: tol = 1.0e-6
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

This example is useful for routines that depend on global arrays because
it:

* Imports the required global arrays from ``vrbls2d``
* Sets the related grid dimensions from ``ctlblk_mod`` before allocating
  the arrays
* Allocates and initializes the global arrays before calling
  ``CALGUSTCONV()``
* Deallocates the arrays after the subroutine call
* Covers typical input, negative wind-speed differences that are treated
  as zero, equal wind speeds, and :term:`spval` handling

.. _example-test-select-channels:

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

This example is useful for tests that call the same routine multiple times
because it:

* Calls ``SELECT_CHANNELS_L()`` multiple times with different inputs
* Resets ``CHANNELINFO%Process_Channel`` between test cases
* Checks both single-value output, such as ``IGOT``, and array values,
  such as ``CHANNELINFO%Process_Channel``
* Prints descriptive messages when a case fails

.. _add-unit-test-upp:

How to Add a Unit Test to UPP
=============================

Add new unit tests to the top-level ``unit_tests`` directory. Use the
standard UPP naming convention for test files:

.. code-block:: text

   test_<name>.f90

where ``<name>`` identifies the subroutine, file, or module under test.

If the routine under test is not part of a module, define an explicit
interface in the unit test:

.. code-block:: fortran

   interface
       subroutine ROUTINE_NAME(arg1, arg2)
           ! Argument declarations
       end subroutine ROUTINE_NAME
   end interface

The structure of each unit test depends on the routine being tested. In
general:

* Routines with single-value inputs should be called once for each test case
* Routines with array inputs can often test multiple cases in a single
  call by assigning each case to a different array index
* Each test case should include a comment that describes the condition
  being tested
* Failure messages should identify what failed and what value was expected
* Unique exit codes should be used where appropriate

After adding the test file, update ``unit_tests/CMakeLists.txt`` so the
test is included in the CMake build:

.. code-block:: cmake

   create_test(test_name)

Use ``create_test(test_name)`` for most unit tests. If the test requires
:term:`MPI`, use:

.. code-block:: cmake

   create_mpi_test(test_name nprocs)

If a test depends on an optional build setting, add logic so the test is
built only when that option is enabled.

.. note::

   New code should be written so that it can be fully tested. Developers
   should be able to write unit tests that achieve 100% line and branch
   coverage for the added code.

.. _update-existing-unit-test:

Updating an Existing Unit Test
==============================

When modifying an existing routine, first determine whether the expected
output is intended to change. If the expected output should remain the
same, the existing unit tests should continue to pass.

If a change adds new logic, such as a new branch in an ``if``/``else``
statement, add a corresponding test case. The exact update depends on the
structure of the existing unit test.

Many UPP unit tests store multiple test cases in input arrays. These tests
often define a parameter near the top of the file that controls the number
of cases. To add a case, increase that parameter, then add the new input
and expected output values.

.. _unit-test-edge-cases:

Common Edge Cases
=================

In addition to covering all branches and error conditions, developers
should consider test cases for common edge conditions, including:

* Threshold boundary values or values just outside thresholds
* Empty string inputs
* Missing files, corrupted formats, permission errors, or other I/O issues
* ``min`` and ``max`` outcomes, especially at boundaries

.. _unit-test-expected-values:

Determining Expected Values
===========================

Expected values usually need to be hardcoded to the appropriate precision.
How those values are determined depends on the output and the calculations
performed by the routine under test.

Possible sources for expected values include:

* Precalculated values
* Existing input or output data
* Proven reference implementations
* Known analytical solutions

Some edge or boundary cases also have predictable expected results.

.. _comparing-floating-point-values:

Comparing Floating-Point Values
===============================

The tolerance used to compare floating-point values depends on the
calculations performed by the routine under test. An absolute tolerance of
``1e-5`` is a reasonable starting point.

Many existing UPP unit tests use ``1e-6``. Use that tighter tolerance only
when the routine can reliably support that level of precision.

.. _unit-tests-ci:

How Unit Tests Run in CI
========================

Unit tests are compiled automatically with CMake when UPP is built. The
top-level CMake configuration adds the ``unit_tests`` directory to the
build, and ``unit_tests`` contains its own ``CMakeLists.txt`` file. No
additional steps are required to build unit tests locally or in GitHub CI.

The ``developer.yml`` :term:`GitHub Actions` workflow runs the unit tests
in the ``run-tests`` step. After all unit tests pass, this workflow
calculates :term:`code coverage` and uploads the coverage report.

Some tests may only compile or run when optional build settings are
enabled. If a new unit test requires a specific build option, add the
appropriate CMake option to the ``build`` step in ``developer.yml`` so the
test is included in the CI coverage results.

.. _view-code-coverage-reports:

Viewing Code Coverage Reports
=============================

Review the test coverage report to confirm that unit tests cover the
intended lines and branches.

The coverage report is available after the ``developer.yml`` workflow
completes successfully. To access the report from an open pull request:

#. Select the ``Checks`` tab under the pull request title.
#. Select the ``developer`` workflow to open the summary view.
#. Find ``UPP-test-coverage`` under ``Artifacts``. The artifact appears
   only after the workflow completes.
#. Download the compressed folder and extract it.
#. Open ``test-coverage.html`` from the extracted folder.
#. Navigate to the file or subroutine being tested.

Lines highlighted in red were not executed by any unit tests. Lines
highlighted in yellow typically indicate partial coverage.

:term:`GCOVR` can occasionally report yellow lines in unexpected ways, so
a yellow line does not always indicate missing coverage. Pay particular
attention to red blocks inside ``if``/``else`` statements. These blocks
may indicate that another test case is needed or that an existing test
case is not reaching the intended branch.

.. _run-unit-tests-locally:

Running Unit Tests Locally
==========================

After building UPP, use :term:`CTest` to run all unit tests from the build
directory:

.. code-block:: console

   ctest --verbose --output-on-failure --rerun-failed

This is the same command used by the GitHub CI workflow. The additional
options provide more output for debugging.

To run a single unit test, execute the test program from the
``<build-directory>/unit_tests/`` subdirectory.

When adding or updating a unit test, developers can usually rebuild it
quickly from the ``<build-directory>/unit_tests/`` subdirectory:

.. code-block:: console

   make

This rebuilds the unit tests without rebuilding all of UPP.

.. _debugging-unit-test-failures:

Debugging Unit Test Failures
============================

GitHub CI reports errors that occur while compiling UPP or running unit
tests. This output is typically similar to the output produced when the
same tests are run locally.

When a unit test fails, CI prints the name of the failing test and its
output. Well-written unit tests should include descriptive error messages
that identify the failed condition and help developers diagnose the issue.

.. _common-unit-test-challenges:

Common Challenges When Writing Unit Tests
=========================================

Writing unit tests for individual UPP routines can be challenging because
many routines were designed to be called by external drivers. Unit tests
run routines independently, so developers may need to recreate state that
is normally set elsewhere in the UPP workflow.

Some UPP routines depend on module variables for global state management,
such as variables from ``CTLBLK.f``. If a routine depends on global
variables, the unit test must set those values explicitly before calling
the routine.

Some routines require data files that are available during normal UPP
execution but not during an isolated unit test. In these cases, update
``unit_tests/CMakeLists.txt`` to copy the required files to the correct
location for the test. See the existing CMake file for examples.
