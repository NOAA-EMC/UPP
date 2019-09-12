# This module looks for environment variables detailing where SFCIO lib is
# If variables are not set, SFCIO will be built from external source 
include(ExternalProject)
if(DEFINED ENV{SFCIO_DIR} )
  message("HEY!! setting SFCIO library via environment variable")
  set(SFCIO_LIBRARY $ENV{SFCIO_LIB} CACHE STRING "SP Library Location" )
  set(SFCIO_4_LIBRARY $ENV{SFCIO_LIB4} CACHE STRING "SFCIO4 Library Location" )
  set(SFCIO_INC $ENV{SFCIO_INC} CACHE STRING "SFCIO d Include Location" )
  set(SFCIO_4_INC $ENV{SFCIO_INC4} CACHE STRING "SFCIO4 Include Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-sfcio 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-sfcio
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-sfcio 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-sfcio OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("sfcio version is ${LIBVERSION}")
  set( SFCIO_LIBRARY ${PROJECT_BINARY_DIR}/lib/libsfcio_${LIBVERSION}.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${SFCIO_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-sfcio )
  else()
      set( CORE_BUILT ${SFCIO_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-sfcio )
  endif()
endif()

set( SFCIO_LIBRARY_PATH ${SFCIO_LIBRARY} CACHE STRING "SFCIO Library Location" )
set( SFCIO_INCLUDE_PATH ${SFCIO_LIBRARY} CACHE STRING "SFCIO Include Location" )

