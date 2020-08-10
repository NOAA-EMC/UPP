set(LIB_SRC
    AllGETHERV_GSD.f
    ALLOCATE_ALL.f
    ASSIGNNEMSIOVAR.f
    AVIATION.f
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
    VRBLS2D_mod.f
    VRBLS3D_mod.f
    VRBLS4D_mod.f
    WETBULB.f
    WETFRZLVL.f
    wrf_io_flags.f
    wrf_io_flags.h
    xml_perl_data.f
    ZENSUN.f)

set(EXE_SRC
    CALCAPE2.f
    CALHEL2.f
    CALVESSEL.f
    GFSPOSTSIG.F
    INITPOST_GFS_NETCDF.f
    MDL2STD_P.f
    TIMEF.f
    WRFPOST.f
    getlvls.f
    intio_tags.f
    io_int_stubs.f)

if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel)$")
  set(CMAKE_Fortran_FLAGS
      "-g -traceback -fp-model source -free -convert big_endian")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
  set(CMAKE_Fortran_FLAGS
      "-g -fbacktrace -ffree-form -ffree-line-length-none -fconvert=big-endian")
  set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
  set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ggdb -fno-unsafe-math-optimizations -frounding-math -fsignaling-nans -ffpe-trap=invalid,zero,overflow -fbounds-check")
endif()

if(APPLE)
  add_compile_definitions(APPLE)
elseif(UNIX)
  add_compile_definitions(LINUX)
endif()

set(LIBNAME nceppost)
set(EXENAME ncep_post)

set(module_dir "${CMAKE_CURRENT_BINARY_DIR}/include")

add_library(${LIBNAME} ${LIB_SRC})
set_target_properties(${LIBNAME} PROPERTIES Fortran_MODULE_DIRECTORY
                                            ${module_dir})
target_include_directories(${LIBNAME} PUBLIC ${NETCDF_INCLUDE_DIRS} ${module_dir})
target_link_libraries(
  ${LIBNAME}
  sp_4
  sigio_4
  sfcio_4
  gfsio_4
  nemsio
  ip_4
  w3emc_4
  w3nco_4
  g2tmpl
  g2_4
  crtm
  bacio_4
  MPI::MPI_Fortran
  PNG::PNG
  ${JASPER_LIBRARY})

if(OpenMP_Fortran_FOUND)
  target_link_libraries(${LIBNAME} OpenMP::OpenMP_Fortran)
endif()

install(
  TARGETS ${LIBNAME}
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib)
install(DIRECTORY ${module_dir} DESTINATION include)

add_executable(${EXENAME} ${EXE_SRC})
target_link_libraries(
  ${EXENAME}
  ${LIBNAME}
  ${NETCDF_LIBRARIES}
  ${CMAKE_DL_LIBS})

install(TARGETS ${EXENAME} RUNTIME DESTINATION bin)
