list(APPEND LIB_SRC
    AllGETHERV_GSD.f
    ALLOCATE_ALL.f
    AVIATION.f
    BNDLYR.f
    BOUND.f
    CALDRG.f
    CALDWP.f
    CALGUST.f
    CALHEL.f
    CALHEL2.f
    CALLCL.f
    CALMCVG.f
    CALMICT.f
    CALPBL.f
    CALPBLREGIME.f
    CALPOT.f
    CALPW.f
    CALRAD_WCLOUD_newcrtm.f
    CALRCH.f
    CALSTRM.f
    CALTAU.f
    CALTHTE.f
    CALUPDHEL.f
    CALVESSEL.f
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
    ETAMP_Q2F.f
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
    get_postfilename.f
    GFIP3.f
    GFSPOST.F
    GPVS.f
    grib2_module.f
    GRIDSPEC.f
    gtg_algo.f90
    gtg_compute.f90
    gtg_config.f90
    gtg_ctlblk.f90
    gtg_filter.f90
    gtg_indices.f90
    gtg_smoothseams.f90
    ICAOHEIGHT.f
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
    MDL2STD_P.f
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
    READCNTRL.F
    READ_xml.f
    retrieve_index.f
    RHGRD.f
    RQSTFLD.F
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
    UPP_MATH.f
    UPP_PHYSICS.f
    VRBLS2D_mod.f
    VRBLS3D_mod.f
    VRBLS4D_mod.f
    WETBULB.f
    WETFRZLVL.f
    xml_perl_data.f
    ZENSUN.f)

list(APPEND EXE_SRC
    ASSIGNNEMSIOVAR.f
    GETNEMSNDSCATTER.f
    GFSPOSTSIG.F
    INITPOST.F
    INITPOST_GFS.f
    INITPOST_GFS_NEMS.f
    INITPOST_GFS_NEMS_MPIIO.f
    INITPOST_GFS_NETCDF.f
    INITPOST_GFS_NETCDF_PARA.f
    INITPOST_GFS_SIGIO.f
    INITPOST_NEMS.f
    INITPOST_NEMS_MPIIO.f
    INITPOST_NETCDF.f
    INITPOST_NMM.f
    TIMEF.f
    WRFPOST.f
    getIVariableN.f
    getVariable.f
    getlvls.f
    intio_tags.f
    native_endianness.f
    retrieve_index.f
    wrf_io_flags.f
    wrf_io_flags.h)

# Use IO stubs in place of WRFLIB
if(NOT wrf_io_FOUND)
  list(APPEND EXE_SRC
    io_int_stubs.f)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS
      "-g -traceback -fp-model source -free -convert big_endian")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS
      "-g -fbacktrace -ffree-form -ffree-line-length-none -fconvert=big-endian")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
  if(CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL 10)
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch")
  endif()
endif()

set(LIBNAME "upp")
set(EXENAME "upp.x")

set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")

add_library(${LIBNAME} STATIC ${LIB_SRC})
add_library(${PROJECT_NAME}::${LIBNAME} ALIAS ${LIBNAME})

set_target_properties(${LIBNAME} PROPERTIES Fortran_MODULE_DIRECTORY
                                            ${module_dir})

target_include_directories(${LIBNAME} INTERFACE
  $<BUILD_INTERFACE:${module_dir}>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(${LIBNAME} PUBLIC
  bacio::bacio_4
  crtm::crtm
  g2::g2_4
  g2tmpl::g2tmpl
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${LIBNAME} PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(BUILD_POSTEXEC)
  add_executable(${EXENAME} ${EXE_SRC})
  target_link_libraries(${EXENAME} PRIVATE
    ${LIBNAME}
    gfsio::gfsio
    w3nco::w3nco_4
    nemsio::nemsio
    ip::ip_4
    sp::sp_4
    sfcio::sfcio
    sigio::sigio)
  if(wrf_io_FOUND)
    target_link_libraries(${EXENAME} PRIVATE
      wrf_io::wrf_io)
  endif()
  install(TARGETS ${EXENAME} RUNTIME DESTINATION bin)
endif()

install(DIRECTORY ${module_dir} DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS ${LIBNAME}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)

### Package config
include(CMakePackageConfigHelpers)
set(CONFIG_INSTALL_DESTINATION lib/cmake/${PROJECT_NAME})

export(EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake)

configure_package_config_file(
  ${CMAKE_SOURCE_DIR}/cmake/PackageConfig.cmake.in ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  INSTALL_DESTINATION ${CONFIG_INSTALL_DESTINATION})
install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

write_basic_package_version_file(
  ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion)
install(FILES ${CMAKE_BINARY_DIR}/${PROJECT_NAME}-config-version.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})

install(EXPORT ${PROJECT_NAME}Exports
  NAMESPACE ${PROJECT_NAME}::
  FILE ${PROJECT_NAME}-targets.cmake
  DESTINATION ${CONFIG_INSTALL_DESTINATION})
