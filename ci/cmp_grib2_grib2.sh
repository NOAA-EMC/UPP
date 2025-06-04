#!/bin/sh
#
# this routine compares two grib2 files
# the files must have the same fields in the same order
#
#set -x

  # for wgrib2 v2.0.7beta2 5/2017
  wgrib2 $2 -var -lev -rpn "sto_1" -import_grib $1 -rpn "rcl_1:print_corr:print_rms" | \
  egrep -v "rpn_corr=1"
#  egrep -v "rpn_corr=1|rpn_rms=undefined"

  # use rpn_corr=0 for all fields
