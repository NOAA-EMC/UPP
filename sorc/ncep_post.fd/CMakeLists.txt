cmake_minimum_required(VERSION 2.6)

set(EXENAME ncep_post )
set(LIBNAME nceppost )
set(LIB_SRC 
AllGETHERV_GSD.f
ALLOCATE_ALL.f
ASSIGNNEMSIOVAR.f
AVIATION.f
blockIO.c
BNDLYR.f
BOUND.f
CALCAPE.f
CALDRG.f
CALDWP.f
CALGUST.f
CALHEL.f
CALLCL.f
CALMCVG.f
CALMICT.f
CALPBL.f
CALPBLREGIME.f
CALPOT.f
CALPW.f
CALRAD_WCLOUD_newcrtm.f
CALRCH.f
CALRH.f
CALRH_GFS.f
CALRH_GSD.f
CALRH_PW.f
CALSTRM.f
CALTAU.f
CALTHTE.f
CALUPDHEL.f
CALVIS.f
CALVIS_GSD.f
CALVOR.f
CALWXT_BOURG.f
CALWXT_DOMINANT.f
CALWXT_EXPLICIT.f
CALWXT.f
CALWXT_RAMER.f
CALWXT_REVISED.f
CANRES.f
CLDFRAC_ZHAO.f
CLDRAD.f
CLMAX.f
CMASSI.f
COLLECT.f
COLLECT_LOC.f
CTLBLK.f
cuparm.f
DEALLOCATE.f
DEWPOINT.f
ETCALC.f
EXCH2.f
EXCH.f
FDLVL.f
FGAMMA.f
FILL_PSETFLD.f
FIXED.f
FRZLVL2.f
FRZLVL.f
GEO_ZENITH_ANGLE.f
GET_BITS.f
GETGBANDSCATTER.f
getIVariableN.f
GETNEMSNDSCATTER.f
get_postfilename.f
getVariable.f
GFIP3.f
GFSPOST.F
GPVS.f
grib2_module.f
GRIBIT.F
GRIDAVG.f
GRIDSPEC.f
gtg_algo.f90
gtg_compute.f90
gtg_config.f90
gtg_ctlblk.f90
gtg_filter.f90
gtg_indices.f90
gtg_smoothseams.f90
ICAOHEIGHT.f
INITPOST.F
INITPOST_GFS.f
INITPOST_GFS_NEMS.f
INITPOST_GFS_NEMS_MPIIO.f
INITPOST_GFS_SIGIO.f
INITPOST_NEMS.f
INITPOST_NEMS_MPIIO.f
INITPOST_NETCDF.f
INITPOST_NMM.f
kinds_mod.F
LFMFLD.f
LFMFLD_GFS.f
LOOKUP.f
machine.f
map_routines.f90
MAPSSLP.f
MASKS_mod.f
MDL2AGL.f
MDL2P.f
MDL2SIGMA2.f
MDL2SIGMA.f
MDL2THANDPV.f
MDLFLD.f
MICROINIT.F
MISCLN.f
MIXLEN.f
MPI_FIRST.f
MPI_LAST.f
MSFPS.f
native_endianness.f
NGMFLD.f
NGMSLP.f
OTLFT.f
OTLIFT.f
PARAMR.f
params.F
PARA_RANGE.f
physcons.f
PMICRPH.f
POLEAVG.f
PROCESS.f
READCNTRL.f
READ_xml.f
retrieve_index.f
RHGRD.f
RQSTFLD.f
SCLFLD.f
SELECT_CHANNELS.f
SERVER.f
SET_LVLSXML.f
SET_OUTFLDS.f
SETUP_SERVERS.f
SLP_new.f
SLP_NMM.f
SMOOTH.f
SNFRAC.f
SNFRAC_GFS.f
SOIL_mod.f
SPLINE.f
SURFCE.f
svptbl.f
TABLE.f
TABLEQ.f
TRPAUS.f
TTBLEX.f
VRBLS2D_mod.f
VRBLS3D_mod.f
VRBLS4D_mod.f
WETBULB.f
WETFRZLVL.f
wrf_io_flags.f
wrf_io_flags.h
WRFPOST.f
xml_perl_data.f
ZENSUN.f)
file(GLOB EXE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.f ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 ${CMAKE_CURRENT_SOURCE_DIR}/*.F)
file(GLOB EXE_C_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
list( REMOVE_ITEM EXE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/TRPAUS_NAM.f)
list( REMOVE_ITEM EXE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/WRF_STUBS.f)
if(IntelComp)
    message("setting intel flags")
    set(CMAKE_Fortran_8_FLAGS " -free -O3 -convert big_endian -traceback -g -fp-model source ${OpenMP_Fortran_FLAGS}" CACHE INTERNAL "" )
    set(CMAKE_C_FLAGS "-DLINUX -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'" CACHE INTERNAL "")
elseif(GNUComp)
    message("setting gnu flags")
    set(CMAKE_Fortran_8_FLAGS " -ffree-form -O3 -fconvert=big-endian -fbacktrace -g ${OpenMP_Fortran_FLAGS}" CACHE INTERNAL "" )
    if(APPLE)
      set(CMAKE_C_FLAGS "-DAPPLE -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'" CACHE INTERNAL "")
    elseif(UNIX)
      set(CMAKE_C_FLAGS "-DLINUX -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'" CACHE INTERNAL "")
    endif(APPLE)
else()
    set(CMAKE_Fortran_8_FLAGS " -ffree-form -O3 -fconvert=big-endian -fbacktrace -g ${OpenMP_Fortran_FLAGS}" CACHE INTERNAL "" )
    set(CMAKE_C_FLAGS "-DLINUX -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'" CACHE INTERNAL "")
endif()

add_executable(${EXENAME} ${EXE_SRC} ${EXE_C_SRC} )
add_library(${LIBNAME} ${EXE_SRC} ${EXE_C_SRC} )
set_source_files_properties( ${EXE_SRC} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_8_FLAGS} )
set_source_files_properties( ${LIB_SRC} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_8_FLAGS} )
set_source_files_properties( ${EXE_C_SRC} PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS} )
set_target_properties( ${EXENAME} PROPERTIES LINK_FLAGS ${OpenMP_Fortran_FLAGS} )
message("hey, g2 inc is ${G2TMPL_INCd}")
include_directories( ${MPI_Fortran_MODULE_DIR} ${MPI_Fortran_INCLUDE_PATH} ${NETCDF_INCLUDES} ${SIGIO_INC} ${NEMSIO_INC} ${SFCIO_INC} ${W3EMC_INC4} ${G2_INC4} ${G2TMPL_INCd} ${GFSIO_INC} ${CRTM_INC} ${CMAKE_INSTALL_PREFIX}/include_4 )
target_link_libraries(${EXENAME} ${BACIO_LIB4} ${SP_LIB4} ${W3NCO_LIB4} ${SIGIO_LIB} ${NEMSIO_LIB} ${IP_LIB4} ${SFCIO_LIB} ${W3EMC_LIB4} ${ZLIB} ${G2TMPL_LIBd} ${G2_LIB4}  ${PNG_LIBRARIES} ${JASPER_LIBRARIES} ${GFSIO_LIB4} ${CRTM_LIB} ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARIES_F90} ${NETCDF_LIBRARIES} ${BACIO_LIB4}  )


install(TARGETS ${EXENAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 

install(TARGETS ${LIBNAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 

install(DIRECTORY ${CMAKE_INCLUDE_OUTPUT_DIRECTORY}/ DESTINATION include_4 )
