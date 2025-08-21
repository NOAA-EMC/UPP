#!/bin/bash

# Set machine-specific variables
source test.bash

case ${machine} in 
   hercules|orion)

      export QUEUE=batch
      export EXCLUSIVE=false
      ;;
   ursa)

      export QUEUE=batch
      export EXCLUSIVE=true
      ;;
   *)
      echo "Unknown machine. Exiting..."
      exit 2;;
esac

# Set test-specific variables
for $test in $test_list
   do
   case $test in
         3drtma|rtma)
            3drtma();;
         gefsv12|gefsv13)
            gefs();;
         gfs)
            gfs();;
         hafs)
            hafs();;
         hrrr)
            hrrr();;
         mpas)
            mpas();;
         mpas_hfip)
            mpas_hfip();;
         nmmb)
            nmmb();;
         rap)
            rap();;
         rrfs)
            rrfs();;
         rrfs_ifi_missing)
            rrfs_ifi_missing();;
         sfs)
            sfs();;
         *) exit 2;;
      esac 

   done
