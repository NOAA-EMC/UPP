.. _unit-tests:

**********
Unit Tests
**********

Overview
========

Automated unit testing helps prevent regressions when existing code is
modified. If a routine is changed without an intentional change in
output, the existing unit tests should continue to pass.

New code should also be fully testable. Developers should be able to
write unit tests that achieve 100% line and branch coverage.

By default, unit tests are compiled automatically with CMake when UPP is
built. The root CMake configuration includes the ``unit_tests``
directory, which contains its own ``CMakeLists.txt`` file. No additional
steps are required to build unit tests locally or in GitHub CI.
