.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

.. _testing-upp:

*****************
Testing the UPP
*****************

Running UPP Regression Tests
=============================

To run the full regression test (RT) suite in preparation for opening a pull request (PR):

   #. Navigate to the local clone of your UPP fork containing the changes you would like to introduce, and run the included RT script within ``/ci``

     .. code-block:: console

        cd /path/to/UPP/ci
        ./rt.sh -a <my_account> -r $PWD/rundir -t $PWD/../ 
        # OR
        nohup ./rt.sh -a <my_account> -r $PWD/rundir -t $PWD/../ &

     where ``my_account`` is the name of an account where you have permissions to run jobs. 

     .. note:: 
        
        The time it takes for tests to run is queue-dependent. RTs can take as little as half an hour to run, but on machines with long queue times, it can take several hours to complete the full set of tests. 

   #. Then, check ``rt.log.<machine>/nohup.out`` for a short summary of any changes in results. 

      * The ``/work`` directory generated in ``UPP/ci`` will include output from each test in case of unexpected errors during runtime within ``out.post.<test_name>`` files. 
      * The ``/rundir`` directory generated within ``UPP/ci`` will include test case results, and ``.diff`` files located within each test's directory will outline changes in fields with the current baselines.
      * Confirm expected changes within the run directory ``.diff`` files if any are present.
      
         * Changes in the ``rap_pe_test`` case only consisting of field 708 Convective Cloud Layer may be ignored; this is a known bug and will always be present within the ``WRFPRS.diff`` file.

Additional Configuration
=========================
For repeated regression test runs, users can edit the ``rt.sh`` file and disable the specified test cases by changing their respective values to “no.” Users can disable the build step as well with the same value for the build variable above the tests. Please be sure to enable all test cases and build settings and conduct a full RT run in preparation for a pull request so that code managers (CMs) can confirm all changes in results are expected and consistent with the developer's results.

``rt.sh`` will allow for changing the configuration of the regression tests if users desire to do so with the following available options:

* ``w`` -- specify the work directory for test case job output
* ``r`` -- specify the run directory containing baselines and .diff files for comparison of changes in results
* ``t`` -- path to the UPP fork to be tested, as ``rt.sh`` is now included within ``UPP/ci``, this should usually be ``$PWD/../``

The following are legacy options for when ``rt.sh`` was not included within the UPP repository and may be ignored by developers: ``-b``, ``-u``, ``-c``