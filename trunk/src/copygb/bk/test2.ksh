#!/bin/ksh

set -x

#NEW=/global/save/wx20gg/bgrids/copygb.fd/copygb
NEW=/meso/save/wx20er/nextjif/nwtest/util/exec/copygb

efile="/nwprod/fix/nam_mxsnoalb.grb"

grid="255 205 954 835 -7491 -144134 136 54000 -106000 126 108 64"
$NEW -g "${grid}" -i3 -x $efile  bfile1.grb

grid2="255 205 954 835 -7491 -144134 136 54000 -106000 126 108 64 44539 14802"
$NEW -g "${grid2}" -i3 -x $efile  bfile2.grb

exit 0
