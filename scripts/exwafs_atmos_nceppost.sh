#####################################################################
echo "-----------------------------------------------------"
echo " exwafs_nceppost.sh" 
echo " Jun 23 - Mao - script for wafs of a simpler version"
echo "-----------------------------------------------------"
#####################################################################
set -x

cd $DATA

# specify model output format type: 4 for nemsio, 3 for sigio
msg="HAS BEGUN on `hostname`"
postmsg "$msg"

export POSTGPSH=${POSTGPSH:-$USHgfs/wafs_nceppost.sh}

export PREFIX=${PREFIX:-${RUN}.t${cyc}z.}

export OUTTYP=${OUTTYP:-4}

if [ $OUTTYP -eq 4 ] ; then
  if [ $OUTPUT_FILE = "netcdf" ]; then
      export SUFFIX=".nc"
      export MODEL_OUT_FORM=${MODEL_OUT_FORM:-netcdfpara}
  else
      export SUFFIX=".nemsio"
      export MODEL_OUT_FORM=${MODEL_OUT_FORM:-binarynemsiompiio}
  fi
else
  export SUFFIX=
fi

export PGBOUT=wafsfile
export PGIOUT=wafsifile

##########################################################
# Specify variable to directly output pgrb2 files for GDAS/GFS
##########################################################
export IDRT=${IDRT:-0} # IDRT=0 is setting for outputting grib files on lat/lon grid


export stime=`echo $post_times | cut -c1-3`
############################################################
# Post Analysis Files before starting the Forecast Post
############################################################
if [ ${stime} = "anl" ]; then
#----------------------------------
    export VDATE=${PDY}${cyc}
    export GFSOUT=${PREFIX}gfsioanl

   if [ $OUTTYP -eq 4 ] ; then
      export loganl=$COMIN/${PREFIX}atmanl${SUFFIX}
   else
      export loganl=$COMIN/${PREFIX}sanl
   fi

   if test -f $loganl ; then

      [[ -f flxfile ]] && rm flxfile ; [[ -f nemsfile ]] && rm nemsfile
      if [ $OUTTYP -eq 4 ] ; then
	  ln -fs $COMIN/${PREFIX}atmanl${SUFFIX} nemsfile
	  export NEMSINP=nemsfile
	  ln -fs $COMIN/${PREFIX}sfcanl${SUFFIX} flxfile
	  export FLXINP=flxfile
      fi

##########################  WAFS U/V/T analysis start ##########################
# U/V/T on ICAO pressure levels for WAFS verification
      if [[ $RUN = gfs && $GRIBVERSION = 'grib2' ]] ; then

	 #For MDL2P.f, WAFS pressure levels are different from master file
	 export POSTGPVARS="KPO=56,PO=84310.,81200.,78190.,75260.,72430.,69680.,67020.,64440.,61940.,59520.,57180.,54920.,52720.,50600.,48550.,46560.,44650.,42790.,41000.,39270.,37600.,35990.,34430.,32930.,31490.,30090.,28740.,27450.,26200.,25000.,23840.,22730.,21660.,20650.,19680.,18750.,17870.,17040.,16240.,15470.,14750.,14060.,13400.,12770.,12170.,11600.,11050.,10530.,10040.,9570.,9120.,8700.,8280.,7900.,7520.,7170.,popascal=.true., numx=1"

	 export PostFlatFile=$PARMpost/postxconfig-NT-GFS-WAFS-ANL.txt
	 export CTLFILE=$PARMpost/postcntrl_gfs_wafs_anl.xml

	 $POSTGPSH=
	 export err=$?

	 if [ $err -ne 0 ] ; then
	     echo " *** GFS POST WARNING: WAFS output failed for analysis, err=$err"
	 else

	    # Need to be saved for WAFS U/V/T verification, 
	    # resolution higher than WAFS 1.25 deg for future compatibility
	    wafsgrid="latlon 0:1440:0.25 90:721:-0.25"
	    $WGRIB2 $PGBOUT -set_grib_type same -new_grid_winds earth \
		    -new_grid_interpolation bilinear -set_bitmap 1 \
		    -new_grid $wafsgrid ${PGBOUT}.tmp

	    if test $SENDCOM = "YES"
	    then
		cp ${PGBOUT}.tmp $COMOUT/${PREFIX}wafs.0p25.anl
		$WGRIB2 -s ${PGBOUT}.tmp > $COMOUT/${PREFIX}wafs.0p25.anl.idx

#         if [ $SENDDBN = YES ]; then
#            $DBNROOT/bin/dbn_alert MODEL GFS_WAFS_GB2 $job $COMOUT/${PREFIX}wafs.0p25.anl
#            $DBNROOT/bin/dbn_alert MODEL GFS_WAFS_GB2__WIDX $job $COMOUT/${PREFIX}wafs.0p25.anl.idx
#         fi
	    fi
	    rm $PGBOUT ${PGBOUT}.tmp
	 fi
      fi
   fi
##########################  WAFS U/V/T analysis end  ##########################
else
##########################  WAFS forecast  start ##########################
   SLEEP_LOOP_MAX=`expr $SLEEP_TIME / $SLEEP_INT`

   for fhr in $post_times ; do

      # Start Looping for the existence of the restart files
      echo 'Start processing fhr='$fhr
      set -x
      ic=1
      while [ $ic -le $SLEEP_LOOP_MAX ]
      do
	  if [  -f $restart_file${fhr}.txt ] ; then
              break
	  else
              ic=`expr $ic + 1`
              sleep $SLEEP_INT
	  fi
	  if [ $ic -eq $SLEEP_LOOP_MAX ] ; then
              echo " *** FATAL ERROR: No model output in nemsio for f${fhr} "
              export err=9
              err_chk
	  fi
      done
      set -x

      export VDATE=`${NDATE} +${fhr} ${PDY}${cyc}`
      [[ -f flxfile ]] && rm flxfile ; [[ -f nemsfile ]] && rm nemsfile
      if [ $OUTTYP -eq 4 ] ; then
	  ln -fs $COMIN/${PREFIX}atmf${fhr}${SUFFIX} nemsfile
	  export NEMSINP=nemsfile
	  ln -fs $COMIN/${PREFIX}sfcf${fhr}${SUFFIX} flxfile
	  export FLXINP=flxfile
      fi
    
      # Generate WAFS products on ICAO standard level.
      # Do not need to be sent out to public, WAFS package will process the data.
      if [[ $fhr -le 120 ]] ; then
	 if [[ ( $RUN = gfs  || $RUN = gefs ) && $GRIBVERSION = 'grib2' ]] ; then

	    #For MDL2P.f, WAFS pressure levels are different from master file
	    export POSTGPVARS="KPO=56,PO=84310.,81200.,78190.,75260.,72430.,69680.,67020.,64440.,61940.,59520.,57180.,54920.,52720.,50600.,48550.,46560.,44650.,42790.,41000.,39270.,37600.,35990.,34430.,32930.,31490.,30090.,28740.,27450.,26200.,25000.,23840.,22730.,21660.,20650.,19680.,18750.,17870.,17040.,16240.,15470.,14750.,14060.,13400.,12770.,12170.,11600.,11050.,10530.,10040.,9570.,9120.,8700.,8280.,7900.,7520.,7170.,gtg_on=.true.,popascal=.true., numx=1"

	    run=`echo $RUN | tr '[a-z]' '[A-Z]'`
	    # Extend WAFS icing and gtg up to 120 hours
            if [  $fhr -le 48  ] ; then
		export PostFlatFile=$PARMpost/postxconfig-NT-${run}-WAFS.txt
		export CTLFILE=$PARMpost/postcntrl_${RUN}_wafs.xml
	    else
		export PostFlatFile=$PARMpost/postxconfig-NT-GFS-WAFS-EXT.txt
		export CTLFILE=$PARMpost/postcntrl_gfs_wafs_ext.xml
	    fi

            # gtg has its own configurations
            cp $HOMEgfs/sorc/ncep_post.fd/post_gtg.fd/gtg.config.$RUN .
            cp $HOMEgfs/sorc/ncep_post.fd/post_gtg.fd/imprintings.gtg_${RUN}.txt .
	    cp $HOMEgfs/sorc/ncep_post.fd/post_gtg.fd/gtg.input.$RUN .

            # WAFS data is processed:
            #   hourly if fhr<=24
            #   every 3 forecast hour if 24<fhr<=48
            #   every 6 forecast hour if 48<fhr<=120
	    if [  $fhr -le 24  ] ; then
		$POSTGPSH
            elif [  $fhr -le 48  ] ; then
		if [  $((10#$fhr%3)) -eq 0  ] ; then
		    $POSTGPSH
		fi
            elif [  $((10#$fhr%6)) -eq 0  ] ; then
		$POSTGPSH
	    fi

            export err=$?

	    if [ $err -ne 0 ] ; then
	       echo " *** GFS POST WARNING: WAFS output failed for f${fhr}, err=$err"
	    else
	       if [ -e $PGBOUT ] ; then
		  if [  $SENDCOM = "YES" ] ; then
		     cp $PGBOUT $COMOUT/${PREFIX}wafs.grb2f$fhr
		     $WGRIB2 -s $PGBOUT > $PGIOUT # WAFS products exist from ush/gfs_nceppost.sh before running anything else
		     cp $PGIOUT $COMOUT/${PREFIX}wafs.grb2f$fhr.idx
		  fi
               fi
	    fi
	 fi
	 [[ -f wafsfile ]] && rm wafsfile ; [[ -f wafsifile ]] && rm wafsifile
      fi
   done
fi
###########################  WAFS forecast end ###########################
echo "PROGRAM IS COMPLETE!!!!!"
