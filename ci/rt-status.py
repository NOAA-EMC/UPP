# check the UPP RT rundir for .diff files that would indicate a change in test results

'''
Updates
Aug/29/2023 - Fernando Andrade-Maldonado: Script creation

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
    rundir = os.getenv('rundir')
    for case_dir in os.listdir(rundir):
        for file in os.listdir(rundir+'/{}'.format(case_dir)):
            if file.endswith('.diff'):
                diff_case = ''
                for test in tests:
                    if test in case_dir:
                        diff_case = test
                if case_dir.endswith('pe_test'):
                    # the rap pe test currently has a false positive bug with WRFPRS
                    if 'rap' in case_dir and file == 'WRFPRS.GrbF16.diff':
                        with open('{}/{}/{}'.format(rundir,case_dir, file)) as f:
                            data = f.readlines()
                            if len(data) == 1 and 'CDCON:convective cloud layer:rpn_corr=-nan:rpn_rms=undefined' in data[0]:
                                continue
                    print('There are changes in results for case {}_pe_test in {}'.format(diff_case, file.replace(".diff", "")))
                else:
                    print('There are changes in results for case {} in {}'.format(diff_case, file.replace(".diff", "")))
                changed = True
    if changed:
        print('Refer to .diff files in rundir: {} for details on differences in results for each case.'.format(rundir))
        sys.exit(1)
    else:
        print('No changes in test results detected.')
                

def main():
    check_for_diff(tests)

if __name__ == "__main__":
    main()
