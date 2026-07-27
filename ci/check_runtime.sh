#!/bin/bash
##########################################################################
# This script is used to retrive runtime on R&D machines.
# Wen Meng, 05/2025, First version.
##########################################################################

#get runtime for each test
export some_failed="NO"
sleep 30
for job_id in $jobid_list; do
  ic=1
  sleep_loop_max=300
  while [ $ic -le $sleep_loop_max ]; do
     status=$(sacct --parsable -j "$job_id" --format=jobid,jobname,elapsed,state | awk -F"|" 'FNR == 2 {print $4}')

     if [ "$status" = "COMPLETED" ]; then
       break
     elif ( echo "$status" | grep -E 'FAIL|TIMEOUT|CANCEL|DEAD|SIGNAL|SPECIAL' > /dev/null ) ; then
       some_failed="YES"
       echo "? Job $job_id failed with status: $status"
       break
     else
      ((ic++))
      sleep 15
     fi
  done

  if [ $ic -lt $sleep_loop_max ]; then

     info=$(sacct --parsable -j "$job_id" --format=jobid,jobname,elapsed,state | awk -F"|" 'FNR == 2')
     runtime_fmt=$(echo "$info" | cut -d"|" -f3)
     jobname=$(echo "$info" | cut -d"|" -f2)

     runtime_b=$(grep -w "^${jobname}" "${runtime_log}" | awk '{print $2}')
     printf "%-10s %-16s %-10s %s\n" "$runtime_fmt" "$jobname" "baseline:" "$runtime_b"
     msg="Runtime: $jobname $runtime_fmt -- baseline ${runtime_b}"
     postmsg "$logfile" "$msg"
  fi
done

elapsed_time=$( printf '%02dh:%02dm:%02ds\n' $((SECONDS%86400/3600)) $((SECONDS%3600/60)) $((SECONDS%60)) )

python ${test_v}/ci/rt-status.py
test_results=$?

if [ "$some_failed" = "YES" ] ; then
  test_results=99
  echo WARNING: some tests exited with non-zero status.
fi

# Cleanup rt log
cd ${test_v}

UPP_HASH=$(git rev-parse HEAD)
SUBMODULE_HASHES=$(git submodule status --recursive)
DATE="$(date '+%Y%m%d %T')"

cd ${test_v}/ci

cat << EOF > ${rt_log}.temp
===== Start of UPP Regression Testing Log =====
UPP Hash Tested:
${UPP_HASH}

Submodule hashes:
${SUBMODULE_HASHES}

Run directory: ${rundir}
Baseline directory: ${homedir}

Total runtime: ${elapsed_time}
Test Date: ${DATE}
Summary Results:

EOF

if [ "$some_failed" = "YES" ] ; then
    echo "Warning: some tests exited with non-zero. status" >> ${rt_log}.temp
    echo >> ${rt_log}.temp
fi

cat ${rt_log} | grep "test:" >> ${rt_log}.temp
cat ${rt_log} | grep "baseline" >> ${rt_log}.temp
python ${test_v}/ci/rt-status.py >> ${rt_log}.temp
echo "===== End of UPP Regression Testing Log =====" >> ${rt_log}.temp
mv ${rt_log}.temp ${rt_log}
mv ${rt_log} ${test_v}/tests/logs
  
# should indicate failure to Jenkins
if [ $test_results -ne 0 ]; then
   python ${test_v}/ci/rt-status.py > changed_results.txt
   if [ "$some_failed" = "YES" ]; then
     echo "Warning: some tests exited with non-zero status." >> changed_results.txt
   fi
   exit 1
fi
