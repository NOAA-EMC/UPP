#!/bin/bash

set_global() {
   
   export WTIME=00:30:00
   export QUEUE=batch

   if [[ ${machine} = "URSA" ]]; then
      export EXCLUSIVE='--exclusive'
   else
      export EXCLUSIVE=''
   fi
   
}

3drtma() {

   case $machine in
      ORION|HERCULES)
         export NODES=8
         export N_TASKS_PER_NODE=12
      ;;
      URSA)
         export WTIME=00:20:00
         export N_TASKS=128
         export TASKS_PER_NODE=32
      ;;
   esac
}

gefs() {
   case $machine in
      ORION|HERCULES)
         export NODES=3
         export N_TASKS_PER_NODE=12
      ;;
      URSA)
         export N_TASKS=48
         export TASKS_PER_NODE=24
      ;;
   esac

}

gfs() {
   case $machine in
      ORION|HERCULES)
         export NODES=6
         export N_TASKS_PER_NODE=40
      ;;
      URSA)
         export N_TASKS=400
         export TASKS_PER_NODE=40
      ;;
   esac

}

hafs() {

   export WTIME=00:20:00

   case $machine in
      ORION|HERCULES)
         export NODES=5
         export N_TASKS_PER_NODE=12
      ;;
      URSA)
         export N_TASKS=72
         export TASKS_PER_NODE=24
         export EXCLUSIVE=''
      ;;
   esac
   
}

hrrr() {

   # Same settings for all machines, unlike most tests
   export WTIME=00:20:00
   export NODES=2
   export N_TASKS_PER_NODE=24

}

mpas() {

   export WTIME=00:20:00

   case $machine in
      ORION|HERCULES)
         export N_TASKS=200
         export TASKS_PER_NODE=40
      ;;
      URSA)
         export N_TASKS=192
         export TASKS_PER_NODE=48
      ;;
   esac

}

mpas_hfip() {

   export EXCLUSIVE='--exclusive'
   export N_TASKS=256
   export CPUS_PER_TASK=4
   
   if [[ $machine = URSA ]]; then
      export MEM='--mem=0'
   fi

}

nmmb() {

   export WTIME=00:20:00

   case $machine in
      ORION|HERCULES)
         export NODES=2
         export N_TASKS_PER_NODE=8
      ;;
      URSA)
         export NODES=7
         export N_TASKS_PER_NODE=4
      ;;
   esac

}

rap() {

   export WTIME=00:20:00

   case $machine in
      ORION|HERCULES)
         export NODES=2
         export N_TASKS_PER_NODE=24
      ;;
      URSA)
         export NODES=4
         export N_TASKS_PER_NODE=12
      ;;
   esac

}

rrfs() {

   case $machine in
      ORION|HERCULES)
         export NODES=6
         export N_TASKS_PER_NODE=40
      ;;
      URSA)
         export N_TASKS=240
         export TASKS_PER_NODE=48
      ;;
   esac

}

rrfs_ifi_missing() {

   case $machine in
      ORION|HERCULES)
         export NODES=8
         export N_TASKS_PER_NODE=12
      ;;
      URSA)
         export N_TASKS=240
         export TASKS_PER_NODE=48
      ;;
   esac

}

sfs() {

   case $machine in
      ORION|HERCULES)
         export NODES=3
         export N_TASKS_PER_NODE=12
      ;;
      URSA)
         export N_TASKS=48
         export TASKS_PER_NODE=24
      ;;
   esac

}
