.. role:: underline
    :class: underline
.. role:: bolditalic
    :class: bolditalic

.. _regression-test:

****************
Regression Tests
****************

Running UPP Regression Tests
=============================

To run the full regression test (RT) suite in preparation for opening a pull request (PR):

   #. Navigate to the local clone of your UPP fork containing the changes you would like to introduce, and run the included RT script within ``ci``:

      .. code-block:: console

         cd /path/to/UPP/ci
         nohup ./rt.sh -a epic -C intel -r $PWD/rundir &

      where ``<my_account>`` is the name of an account where you have permissions to run jobs. The terminal will print a message like:

      .. code-block:: console
        
         nohup: ignoring input and appending output to ‘nohup.out’
     
      The user can continue to issue commands in the Terminal while the RTs run in the background. 

      .. note:: 
        
         The time it takes for tests to run is queue-dependent. RTs can take as little as half an hour to run, but on machines with long queue times, it can take a few hours to complete the full set of tests. 

   #. Check ``nohup.out`` or ``UPP/tests/logs/rt.log.<machine>_<compiler>`` for a short summary of any changes in results. The tests are finished when there are 17 timestamps and a final results summary (e.g., "No changes in test results detected."). 

      * The ``work`` directory generated in ``UPP/ci`` contains ``out.post.<test_name>`` files, which list output from each test, including any unexpected errors during runtime. 
      * The ``rundir`` directory generated within ``UPP/ci`` will include test case results, and ``.diff`` files located within each test's directory will outline changes in fields with the current baselines.
      * Confirm that any changes within the run directory ``.diff`` files are expected if any are present.
   
   #. Check for errors in the RT output directory (e.g., ``work-upp-<machine>``) using the following commands:

      .. code-block:: console

         cd work-upp-<machine>
         grep -ir "error" .
         grep -ir "fatal" .
   
   #. Push the test log into your local branch:

      .. code-block:: console

         cd /path/to/UPP/tests/logs
         git add rt.log.<machine>_<compiler>
         git commit -m "<machine> <compiler> rt log"
         git push origin <branch>

Additional Configuration
=========================
For repeated regression test runs, users can edit the ``rt.sh`` file and remove any tests that should not be run from ``test_list``. However, please be sure to enable all test cases and build settings and conduct a full RT run in preparation for a pull request so that code managers (CMs) can confirm all changes in results are expected and consistent with the developer's results.

``rt.sh`` will allow for changing the configuration of the regression tests if users desire to do so with the following available options:

* ``w`` -- specify the work directory for test case job output
* ``r`` -- specify the run directory containing baselines and ``.diff`` files for comparison of changes in results
* ``d`` -- disable IFI tests even if IFI is available
* ``e`` -- do not build the UPP executable

The following are legacy options from before ``rt.sh`` was included within the UPP repository; they may be ignored by developers: ``-b``, ``-u``, ``-c``, ``-t``.
