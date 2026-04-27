# This is the cmake build file for the tests directory of the UPP project.
#
# Alyson Stahl, 10/2025

message(STATUS "Building test directory...")

# This function builds, links, and runs a Fortran test program.
function(create_test name)
  # optional source files for testing files that are only part of the EXE source list.
  cmake_parse_arguments(ARG "" "" "SRC_FILES" ${ARGN})
  if(ARG_SRC_FILES)
    if(CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|IntelLLVM)$")
      set_source_files_properties(${ARG_SRC_FILES} PROPERTIES COMPILE_FLAGS "-free")
    elseif(CMAKE_Fortran_COMPILER_ID MATCHES "^(GNU|Clang|AppleClang)$")
      set_source_files_properties(${ARG_SRC_FILES} PROPERTIES COMPILE_FLAGS "-ffree-form")
    endif()
  endif()
  add_executable(${name} ${name}.f90 ${ARG_SRC_FILES})
  target_link_libraries(${name} PRIVATE upp::upp)
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE Fortran)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

# This function builds, links, and runs a Fortran test program using MPI with a specified number of processes.
function(create_mpi_test name num_procs)
  add_executable(${name}_np_${num_procs} ${name}.f90)
  target_link_libraries(${name}_np_${num_procs} PRIVATE upp::upp MPI::MPI_Fortran)
  set_target_properties(${name}_np_${num_procs} PROPERTIES LINKER_LANGUAGE Fortran)

  add_test(
    NAME ${name}_np_${num_procs}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${num_procs}
            ${MPIEXEC_PREFLAGS} $<TARGET_FILE:${name}_np_${num_procs}> ${MPIEXEC_POSTFLAGS}
  )
  set_tests_properties(${name}_np_${num_procs} PROPERTIES PROCESSORS ${num_procs})
endfunction()

create_test(test_dewpoint)
create_test(test_bound)
create_test(test_caldrg)
create_test(test_calstrm)
create_test(test_callcl)
create_test(test_calpot)
create_test(test_calvis)
create_test(test_calvis_gsd)
create_test(test_calwxt_dominant)
create_test(test_calthte)
create_test(test_calvessel)
create_test(test_calwxt_explicit)
create_test(test_ttblex)
create_test(test_wetfrzlvl)
create_test(test_zensun)
create_test(test_geo_zenith_angle)
create_test(test_icaoheight)
create_test(test_gpvs)
create_test(test_etamp_q2f)
create_test(test_etcalc)
create_test(test_fgamma)
create_test(test_canres)
create_test(test_clmax)
create_test(test_cldfrac_zhao)
create_test(test_mdl2sigma2)
create_test(test_mixlen)
create_test(test_msfps)
create_test(test_para_range)
create_test(test_poleavg)
create_test(test_smooth)
create_test(test_sclfld)
create_test(test_spline)
create_test(test_snfrac)
create_test(test_snfrac_gfs)
create_test(test_set_lvlsxml)
create_test(test_select_channels)
create_test(test_caldwp)

if(BUILD_WITH_NEMSIO)
  create_test(test_assignnemsiovar SRC_FILES ${CMAKE_SOURCE_DIR}/sorc/ncep_post.fd/ASSIGNNEMSIOVAR.f)
endif()

# The following tests require MPI to run.
create_mpi_test(test_allgetherv 1)
create_mpi_test(test_allgetherv 2)