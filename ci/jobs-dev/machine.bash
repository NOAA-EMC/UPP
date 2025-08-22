#!/bin/bash

# Set machine-specific variables

if [[ ${machine} = "URSA" ]]; then
   export EXCLUSIVE='--exclusive'
   export STACK_SIZE=''
elif
   export STACK_SIZE='ulimit -s unlimited\nexport OMP_STACKSIZE=128M'
   export EXCLUSIVE=''
else
   echo "${machine} is not a supported machine. Cannot set machine-specific variables."
fi