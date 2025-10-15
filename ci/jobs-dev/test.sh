#!/bin/bash

set_global() {
   
   export WTIME=00:30:00
   export QUEUE=batch

   # Reset variables between tests
   export NODES=''
   export N_TASKS_PER_NODE=''
   export N_TASKS=''
   export TASKS_PER_NODE=''
   export CPUS_PER_TASK=''
   export NUMX=''
   export MEM=''

}

3drtma() {

   case $machine in
      ORION|HERCULES)
         export NODES='-N 8'
         export N_TASKS_PER_NODE='--ntasks-per-node=12'
      ;;
      URSA)
         export WTIME=00:20:00
         export N_TASKS='--ntasks 128'
         export TASKS_PER_NODE='--tasks-per-node 32'
      ;;
   esac
}

gefsv12() {
   case $machine in
      ORION|HERCULES)
         export NODES='-N 3'
         export N_TASKS_PER_NODE='--ntasks-per-node=12'
      ;;
      URSA)
         export N_TASKS='--ntasks 48'
         export TASKS_PER_NODE='--tasks-per-node 24'
      ;;
   esac

}

gefsv13() {
   case $machine in
      ORION|HERCULES)
         export NODES='-N 3'
         export N_TASKS_PER_NODE='--ntasks-per-node=12'
      ;;
      URSA)
         export N_TASKS='--ntasks 48'
         export TASKS_PER_NODE='--tasks-per-node 24'
      ;;
   esac

}

gfs() {
   case $machine in
      ORION|HERCULES)
         export NODES='-N 3'
         export N_TASKS_PER_NODE='--ntasks-per-node=40'
      ;;
      URSA)
         export N_TASKS='--ntasks 120'
         export TASKS_PER_NODE='--tasks-per-node 40'
      ;;
   esac

}

hafs() {

   export WTIME=00:20:00

   case $machine in
      ORION|HERCULES)
         export NODES='-N 5'
         export N_TASKS_PER_NODE='--ntasks-per-node=12'
      ;;
      URSA)
         export N_TASKS='--ntasks 72'
         export TASKS_PER_NODE='--tasks-per-node 24'
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
         export NUMX=',numx=4'
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
         export NODES='-N 6'
         export N_TASKS_PER_NODE='--ntasks-per-node=40'
      ;;
      URSA)
         export N_TASKS='--ntasks 240'
         export TASKS_PER_NODE='--tasks-per-node 48'
      ;;
   esac

}

rrfs_ifi_missing() {

   case $machine in
      ORION|HERCULES)
         export NODES='-N 8'
         export N_TASKS_PER_NODE='--ntasks-per-node=12'
      ;;
      URSA)
         export N_TASKS='--ntasks 240'
         export TASKS_PER_NODE='--tasks-per-node 48'
      ;;
   esac

}

sfs() {

   case $machine in
      ORION|HERCULES)
         export NODES='-N 3'
         export N_TASKS_PER_NODE='--ntasks-per-node=12'
         export STACK_SIZE='export OMP_STACKSIZE=128M'
      ;;
      URSA)
         export N_TASKS='--ntasks 48'
         export TASKS_PER_NODE='--tasks-per-node 24'
      ;;
   esac

}
