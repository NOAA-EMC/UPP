# This is the cmake build file for the tests directory of the UPP project.
#
# Alyson Stahl, 10/2025

message(STATUS "Building test directory...")

# This function builds, links, and runs a Fortran test program.
function(create_test name)
  add_executable(${name} ${name}.f90)
  target_link_libraries(${name} PRIVATE upp::upp)
  set_target_properties(${name} PROPERTIES LINKER_LANGUAGE Fortran)
  add_test(NAME ${name} COMMAND ${name})
endfunction()

create_test(test_dewpoint)