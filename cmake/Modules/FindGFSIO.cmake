# This module looks for environment variables detailing where GFSIO lib is
# If variables are not set, GFSIO will be built from external source 
include(ExternalProject)
if(DEFINED ENV{GFSIO_DIR} )
  message("HEY!! setting GFSIO library via environment variable")
  set(GFSIO_LIBRARY $ENV{GFSIO_LIB4} CACHE STRING "GFSIO Library Location" )
  set(GFSIO_INC $ENV{GFSIO_INC4} CACHE STRING "GFSIO_4 Include Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-gfsio 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-gfsio
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-gfsio 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-gfsio OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("gfsio version is ${LIBVERSION}")
  set( GFSIO_LIBRARY ${PROJECT_BINARY_DIR}/lib/libgfsio_${LIBVERSION}_d.a )
  set( GFSIO_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libgfsio_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${GFSIO_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-gfsio )
  else()
      set( CORE_BUILT ${GFSIO_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-gfsio )
  endif()
  
  set( GFSIO_LIBRARY_PATH ${GFSIO_LIBRARY} CACHE STRING "GFSIO Library Location" )
  set( GFSIO_INCLUDE_PATH ${GFSIO_INC} CACHE STRING "GFSIO Include Location" )
  
endif()
