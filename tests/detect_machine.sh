#!/bin/bash

case $(hostname -f) in

  alogin0[12].acorn.wcoss2.ncep.noaa.gov)    MACHINE_ID=wcoss2_a ;; ### acorn
  adecflow0[12].acorn.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2_a ;; ### acorn
  dlogin0[1-9].dogwood.wcoss2.ncep.noaa.gov) MACHINE_ID=wcoss2 ;; ### dodwood
  clogin0[1-9].cactus.wcoss2.ncep.noaa.gov)  MACHINE_ID=wcoss2 ;; ### cactus

  gaea5[1-8])                    MACHINE_ID=gaea ;; ### gaea51-58
  gaea5[1-8].ncrc.gov)           MACHINE_ID=gaea ;; ### gaea51-58

  hfe0[1-8])                   MACHINE_ID=hera_c ;; ### hera01-08
  hfe09)                       MACHINE_ID=hera ;; ### hera09
  hfe1[0-2])                   MACHINE_ID=hera ;; ### hera10-12

  fe[1-4])                     MACHINE_ID=jet_c ;; ### jet1-4
  fe[5-8])                     MACHINE_ID=jet ;; ### jet5-8
  tfe[12])                     MACHINE_ID=jet ;; ### tjet1-2

  Orion-login-[1-4].HPC.MsState.Edu) MACHINE_ID=orion ;; ### orion1-4

  [Hh]ercules-login-[1-4].[Hh][Pp][Cc].[Mm]s[Ss]tate.[Ee]du) MACHINE_ID=hercules ;; ### hercules1-4

  login[1-4].stampede2.tacc.utexas.edu) MACHINE_ID=stampede ;; ### stampede1-4

  s4-submit.ssec.wisc.edu) MACHINE_ID=s4 ;; ### S4

  *) MACHINE_ID=unknown
esac

# Overwrite auto-detect with RT_MACHINE if set
MACHINE_ID=${RT_MACHINE:-${MACHINE_ID}}

# Append compiler
#if [ $MACHINE_ID = orion ] || [ $MACHINE_ID = hera ] || [ $MACHINE_ID = cheyenne ] || [ $MACHINE_ID = jet ] || [ $MACHINE_ID = gaea ] || [ $MACHINE_ID = stampede ] || [ $MACHINE_ID = s4 ]; then
#    MACHINE_ID=${MACHINE_ID}.${RT_COMPILER}
#fi

echo "Machine: " $MACHINE_ID 
