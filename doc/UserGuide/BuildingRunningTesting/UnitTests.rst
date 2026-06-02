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
