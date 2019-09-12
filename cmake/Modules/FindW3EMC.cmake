# This module looks for environment variables detailing where W3EMC lib is
# If variables are not set, W3EMC will be built from external source 
include(ExternalProject)
if(DEFINED ENV{W3EMC_LIBd} )
  message("HEY!! setting W3EMC library via environment variable")
  set(W3EMC_LIBRARY $ENV{W3EMC_LIBd} CACHE STRING "SP Library Location" )
  set(W3EMC_4_LIBRARY $ENV{W3EMC_LIB4} CACHE STRING "W3EMC4 Library Location" )
  set(W3EMCINC $ENV{W3EMC_INC} CACHE STRING "W3EMC d Include Location" )
  set(W3EMC4INC $ENV{W3EMC_INC4} CACHE STRING "W3EMC4 Include Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-w3emc 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-w3emc
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DEXTRA_INCLUDE=${PROJECT_BINARY_DIR}/include
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-w3emc 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  add_dependencies(NCEPLIBS-w3emc NCEPLIBS-sigio)
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-w3emc OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("w3emc version is ${LIBVERSION}")
  set( W3EMC_LIBRARY ${PROJECT_BINARY_DIR}/lib/libw3emc_${LIBVERSION}_d.a )
  set( W3EMC_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libw3emc_${LIBVERSION}_4.a )
  set(W3EMCINC ${PROJECT_BINARY_DIR}/include )
  set(W3EMC4INC ${PROJECT_BINARY_DIR}/include_4 )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${W3EMC_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-w3emc )
  else()
      set( CORE_BUILT ${W3EMC_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-w3emc )
  endif()
endif()

set( W3EMC_LIBRARY_PATH ${W3EMC_LIBRARY} CACHE STRING "W3EMC Library Location" )
set( W3EMC_4_LIBRARY_PATH ${W3EMC_4_LIBRARY} CACHE STRING "W3EMC_4 Library Location" )
set( W3EMC_INCLUDE_PATH ${W3EMCINC} CACHE STRING "W3EMC Include Location" )
set( W3EMC_INCLUDE_4_PATH ${W3EMC4INC} CACHE STRING "W3EMC_4 Include Location" )

