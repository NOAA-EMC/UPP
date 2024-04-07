set(PERLXML_SCRIPT "${CMAKE_CURRENT_SOURCE_DIR}/PostXMLPreprocessor.pl")

# Function to split a string into a list based on a delimiter
function(split_delimited_string_into_list INPUT_STRING DELIMITER OUTPUT_LIST)
  string(REPLACE "${DELIMITER}" ";" TEMP_LIST "${INPUT_STRING}")
  set(${OUTPUT_LIST} ${TEMP_LIST} PARENT_SCOPE)
endfunction()

# Function to get filename from path, and replace '-' with '_'
# The replacing is done because cmake targets cannot have '-'
# The filename is used to name the target e.g. postxconfig_NT_gfs
function(convert_filename FILE_PATH OUTPUT_VAR)
  get_filename_component(FILE_NAME_WE ${FILE_PATH} NAME_WE)
  string(REPLACE "-" "_" NEW_FILE_NAME ${FILE_NAME_WE})
  set(${OUTPUT_VAR} ${NEW_FILE_NAME} PARENT_SCOPE)
endfunction()

# List of output, input and available files
list(APPEND ALL_TXT
  "./gfs/postxconfig-NT-gfs.txt,./gfs/postcntrl_gfs.xml,post_avblflds.xml"
  "./gfs/postxconfig-NT-gfs-goes.txt,./gfs/postcntrl_gfs_goes.xml,post_avblflds.xml"
  "./gfs/postxconfig-NT-gfs-anl.txt,./gfs/postcntrl_gfs_anl.xml,post_avblflds.xml"
  "./gfs/postxconfig-NT-gfs-f00.txt,./gfs/postcntrl_gfs_f00.xml,post_avblflds.xml"
  "./gfs/postxconfig-NT-gfs-flux.txt,./gfs/postcntrl_gfs_flux.xml,post_avblflds.xml"
  "./gfs/postxconfig-NT-gfs-flux-f00.txt,./gfs/postcntrl_gfs_flux_f00.xml,post_avblflds.xml"
  "./gfs/postxconfig-NT-gfs-two.txt,./gfs/postcntrl_gfs_two.xml,post_avblflds.xml"
  "./gfs/postxconfig-NT-gfs-f00-two.txt,./gfs/postcntrl_gfs_f00_two.xml,post_avblflds.xml"
  "./gefs/postxconfig-NT-gefs.txt,./gefs/postcntrl_gefs.xml,post_avblflds.xml"
  "./gefs/postxconfig-NT-gefs-f00.txt,./gefs/postcntrl_gefs_f00.xml,post_avblflds.xml"
  "./gefs/postxconfig-NT-gefs-aerosol.txt,./gefs/postcntrl_gefs_aerosol.xml,post_avblflds.xml"
  "./gefs/postxconfig-NT-gefs-f00-aerosol.txt,./gefs/postcntrl_gefs_aerosol_f00.xml,post_avblflds.xml"
  "./sfs/postxconfig-NT-sfs.txt,./sfs/postcntrl_sfs.xml,post_avblflds.xml"
  "./postxconfig-NT-NMM.txt,./nam_cntrl_cmaq.xml,nam_post_avblflds.xml"
  "./postxconfig-NT-NGAC.txt,./ngac_postcntrl.xml,ngac_post_avblflds.xml"
  "./postxconfig-NT-hafs_sat.txt,./postcntrl_hafs_sat.xml,post_avblflds.xml"
  "./postxconfig-NT-hafs_nosat.txt,./postcntrl_hafs_nosat.xml,post_avblflds.xml"
  "./postxconfig-NT-hafs.txt,./postcntrl_hafs.xml,post_avblflds.xml"
  "./postxconfig-NT-hrrr.txt,./hrrr_postcntrl.xml,post_avblflds_raphrrr.xml"
  "./postxconfig-NT-rap.txt,./rap_postcntrl.xml,post_avblflds_raphrrr.xml"
  "./postxconfig-NT-fv3lam_rrfs.txt,./fv3lam_rrfs.xml,post_avblflds.xml"
  "./postxconfig-NT-UFS-aerosol.txt,./postcntrl_ufs_aerosol.xml,post_avblflds.xml"
  "./postxconfig-NT-UFS-aerosol-F00.txt,./postcntrl_ufs_aerosol_f00.xml,post_avblflds.xml"
  "./postxconfig-NT-AQM.txt,./aqm.xml,post_avblflds.xml"
)

# Loop over the list of lists and create a target for each line
foreach(TXT_STRING IN LISTS ALL_TXT)

  # Split each line into a list and get individual values
  split_delimited_string_into_list(${TXT_STRING} "," LIST_STRING)
  list(GET LIST_STRING 0 target_txt)
  list(GET LIST_STRING 1 control_xml)
  list(GET LIST_STRING 2 avblflds_xml)

  # Name the target based on the target txt file (remove path and replace - with _)
  convert_filename(${target_txt} target_name)

  # Using CMAKE_CURRENT_SOURCE_DIR as WORKING_DIRECTORY is not a CMake best practice
  # To eliminate this, the perl script will need modifications
  add_custom_command(
    OUTPUT ${target_txt}
    COMMAND ${PERL_EXECUTABLE} ${PERLXML_SCRIPT} "${control_xml}" "${avblflds_xml}" "${target_txt}"
    DEPENDS ${control_xml} ${avblflds_xml}
    WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
    COMMENT "Processing ${target_txt}"
    VERBATIM
  )

  # Create a make target to run the above custom command
  add_custom_target(${target_name} DEPENDS ${target_txt})  # Do not add to "ALL" targets

endforeach()

