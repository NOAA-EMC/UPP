# This module looks for environment variables detailing where CRTM lib is
# If variables are not set, CRTM will be built from external source 
include(ExternalProject)
if(DEFINED ENV{CRTM_LIB} )
  message("HEY!! setting CRTM library via environment variable")
  set(CRTM_LIBRARY $ENV{CRTM_LIB} CACHE STRING "CRTM Library Location" )
  set(CRTM_INC $ENV{CRTM_INC} CACHE STRING "CRTM Include Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-crtm 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-crtm
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-crtm 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-crtm OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("crtm version is ${LIBVERSION}")
  set( CRTM_LIBRARY ${PROJECT_BINARY_DIR}/lib/libcrtm_${LIBVERSION}_d.a )
  set( CRTM_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libcrtm_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${CRTM_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-crtm )
  else()
      set( CORE_BUILT ${CRTM_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-crtm )
  endif()
  
endif()

set( CRTM_LIBRARY_PATH ${CRTM_LIBRARY} CACHE STRING "CRTM Library Location" )
set( CRTM_INCLUDE_PATH ${CRTM_INC} CACHE STRING "CRTM Include Location" )
