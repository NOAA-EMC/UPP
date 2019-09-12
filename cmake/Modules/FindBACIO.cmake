# This module looks for environment variables detailing where BACIO lib is
# If variables are not set, BACIO will be built from external source 
include(ExternalProject)
if(DEFINED ENV{BACIO_LIB4} )
  set(BACIO_LIBRARY $ENV{BACIO_LIB4} CACHE STRING "BACIO Library Location" )
  set(BACIO_VER $ENV{BACIO_VER} CACHE STRING "BACIO Version")
  set(BACIO_8_LIBRARY $ENV{BACIO_LIBd} CACHE STRING "BACIO_8 Library Location")
  set(BACIO_d_LIBRARY $ENV{BACIO_LIBd} CACHE STRING "BACIO_8 Library Location")
else()  
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-bacio 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-bacio
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-bacio 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-bacio OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("bacio version is ${LIBVERSION}")
  set( BACIO_LIBRARY ${PROJECT_BINARY_DIR}/lib/libbacio_${LIBVERSION}_4.a )
  set( BACIO_8_LIBRARY ${PROJECT_BINARY_DIR}/lib/libbacio_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${BACIO_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-bacio )
  else()
      set( CORE_BUILT ${BACIO_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-bacio )
  endif()

  set( BACIO_LIBRARY_PATH ${BACIO_LIBRARY} CACHE STRING "BACIO Library Location" )
endif()
