#!/bin/bash

# Set machine-specific variables

if [[ ${machine} = "URSA" ]]; then
   export EXCLUSIVE='--exclusive'
   export STACK_SIZE=''
elif [[ ${machine} = "ORION" || ${machine} = "HERCULES" ]]; then
   ulimit -s unlimited
   export EXCLUSIVE=''
else
   echo "${machine} is not a supported machine. Cannot set machine-specific variables."
fi
