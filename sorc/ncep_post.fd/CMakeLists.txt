list(APPEND LIB_SRC
    AllGETHERV_GSD.f
    ALLOCATE_ALL.f
    AVIATION.f
    BNDLYR.f
    BOUND.f
    CALCHIPSI.f
    CALDRG.f
    CALDWP.f
    CALGUST.f
    CALHEL.f
    CALHEL2.f
    CALHEL3.f
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
    EXCH.f
    EXCH_c_float.f
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
    ICAOHEIGHT.f
    IFI.F
    kinds_mod.F
    LFMFLD.f
    LFMFLD_GFS.f
    LOOKUP.f
    machine.f
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

list(APPEND GTG_LIB_SRC
  gtg_algo.F90
  gtg_compute.F90
  gtg_cdf_mod.F90
  gtg_ReadSetupParams.F90
  gtg_read_config.F90
  gtg_mlmodel.F90
  gtg_prob.F90
  gtg_interp.F90
  gtg_ctlblk.F90
  gtg_filter.F90
  gtg_indices.F90
  gtg_smoothseams.F90
  map_routines.F90)

if(BUILD_WITH_GTG)
  list(TRANSFORM GTG_LIB_SRC PREPEND post_gtg.fd/)
endif()

list(APPEND LIB_SRC ${GTG_LIB_SRC})

list(APPEND EXE_SRC
    INITPOST.F
    INITPOST_NETCDF.f
    WRFPOST.F
    getIVariableN.f
    getVariable.f
    getlvls.f
    intio_tags.f
    retrieve_index.f
    wrf_io_flags.f
    wrf_io_flags.h)

# Use IO stubs in place of WRFLIB
if(NOT wrf_io_FOUND)
  list(APPEND EXE_SRC
    io_int_stubs.f)
endif()

if(BUILD_WITH_NEMSIO)
  list(APPEND EXE_SRC
    ASSIGNNEMSIOVAR.f
    GETNEMSNDSCATTER.f
    GFSPOSTSIG.F
    INITPOST_GFS_NEMS_MPIIO.f
    INITPOST_NEMS.f)
endif()

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS
      "-g -traceback -fp-model source -free -convert big_endian")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ftrapuv -check all -check noarg_temp_created")
  set_source_files_properties(INITPOST_GFS_NEMS_MPIIO.f INITPOST_NETCDF.f INITPOST_NEMS.f PROPERTIES COMPILE_FLAGS -qoverride-limits)
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

target_include_directories(${LIBNAME} PUBLIC
  $<BUILD_INTERFACE:${module_dir}>
  $<INSTALL_INTERFACE:include>)

target_link_libraries(${LIBNAME} PUBLIC
  bacio::bacio_4
  crtm::crtm
  g2::g2_4
  g2tmpl::g2tmpl
  MPI::MPI_Fortran
  NetCDF::NetCDF_Fortran)

if(IFI_FOUND)
  if(INTERNAL_IFI)
    target_link_libraries(${LIBNAME} PRIVATE IFI)
    target_link_libraries(${LIBNAME} PRIVATE tdrp)
    target_include_directories(${LIBNAME} PUBLIC $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>/include $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>)
  else()
    target_link_libraries(${LIBNAME} PUBLIC IFI)
    target_link_libraries(${LIBNAME} PUBLIC tdrp)
  endif()
endif()

if(BUILD_WITH_GTG)
  target_link_libraries(${LIBNAME} PUBLIC
  ip::ip_4)
endif()

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${LIBNAME} PUBLIC OpenMP::OpenMP_Fortran)
endif()

if(BUILD_POSTEXEC)
  add_executable(${EXENAME} ${EXE_SRC})
  target_link_libraries(${EXENAME} PRIVATE
    ${LIBNAME}
    w3emc::w3emc_4
    sp::sp_4)
  set_target_properties(${EXENAME} PROPERTIES LINKER_LANGUAGE Fortran)
  if(IFI_FOUND)
    target_link_libraries(${EXENAME} PRIVATE IFI)
    target_link_libraries(${EXENAME} PRIVATE tdrp)
  endif()
  if(BUILD_WITH_NEMSIO)
    set_target_properties(nemsio::nemsio PROPERTIES
      INTERFACE_LINK_LIBRARIES "w3emc::w3emc_4;bacio::bacio_4")
    target_link_libraries(${EXENAME} PRIVATE
      nemsio::nemsio sigio::sigio)
    target_compile_definitions(${EXENAME} PRIVATE BUILD_WITH_NEMSIO)
  endif()
  if(wrf_io_FOUND)
    target_link_libraries(${EXENAME} PRIVATE
      wrf_io::wrf_io)
  endif()
  install(TARGETS ${EXENAME} RUNTIME DESTINATION ${exec_dir})
endif()

install(DIRECTORY ${module_dir} DESTINATION ${CMAKE_INSTALL_PREFIX})

install(
  TARGETS ${LIBNAME}
  EXPORT ${PROJECT_NAME}Exports
  RUNTIME DESTINATION ${exec_dir}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
