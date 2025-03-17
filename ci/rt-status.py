# check the UPP RT rundir for .diff files that would indicate a change in test results

'''
Updates
Aug/29/2023 - Fernando Andrade-Maldonado: Script creation
Feb/10/2025 - Sam Trahan: Detect if a test was unable to run.
'''

import os
import sys

# files used in result comparison
tests = [
    'nmmb',
    'gfs',
    'fv3gefs',
    'fv3r',
    'rap',
    'hrrr',
    'fv3hafs',
    'rtma'
]

# look for .diff files
# every case has its own directory in rundir
# loop through every test case sub directory and files, then match with the test name
def check_for_diff(tests):
    changed = False
    failed = False
    success = False
    rundir = os.getenv('rundir')
    for case_dir in os.listdir(rundir):
        full_case_dir=rundir+'/'+case_dir
        case_failed = False
        case_success = False
        diff_case = ''
        for file in os.listdir(full_case_dir):
            if file.endswith('.diff'):
                for test in tests:
                    if test in case_dir:
                        diff_case = test
        for file in os.listdir(full_case_dir):
            full_file = full_case_dir + '/' + file
            if file == 'TEST_ERROR':
                case_failed = True
                print('Error: Test case {} was unable to run! See {} for details.'
                      .format(diff_case, full_case_dir + '/TEST_ERROR'))
            elif file == 'SUCCESS':
                case_success = True
            elif file.endswith('.diff'):
                checked_file = full_case_dir + '/' + file.replace(".diff", "")
                if case_dir.endswith('pe_test'):
                    # the rap pe test currently has a false positive bug with WRFPRS
                    if 'rap' in case_dir and file == 'WRFPRS.GrbF16.diff':
                        with open(full_file) as f:
                            data = f.readlines()
                            if len(data) == 1 and 'CDCON:convective cloud layer:rpn_corr=-nan:rpn_rms=undefined' in data[0]:
                                continue
                    print('There are changes in results for case {}_pe_test in {}'.format(diff_case, checked_file))
                else:
                    print('There are changes in results for case {} in {}'.format(diff_case, checked_file))
                changed = True
        if not case_success and not case_failed:
            print('Error: Test case {} did not produce SUCCESS flag file, but provided no reason for failing. '
                  'Check {} for trouble.'.format(diff_case,rundir + '/' + case_dir))
        failed = case_failed or failed
        success = case_success and success
    if changed:
        print('Refer to .diff files in rundir: {} for details on differences in results for each case.'.format(rundir))
    if failed:
        print('Error: Some tests were unable to run!')
    if failed or changed:
        sys.exit(1)
    else:
        print('No changes in test results detected.')

def main():
    check_for_diff(tests)

if __name__ == "__main__":
    main()
