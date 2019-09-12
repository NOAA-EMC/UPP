# This module looks for environment variables detailing where SIGIO lib is
# If variables are not set, SIGIO will be built from external source 
include(ExternalProject)
if(DEFINED ENV{SIGIO_DIR} )
  message("HEY!! setting SIGIO library via environment variable")
  set(SIGIO_LIBRARY $ENV{SIGIO_LIB} CACHE STRING "SP Library Location" )
  set(SIGIO_4_LIBRARY $ENV{SIGIO_LIB4} CACHE STRING "SIGIO4 Library Location" )
  set(SIGIO_INC $ENV{SIGIO_INC} CACHE STRING "SIGIO d Include Location" )
  set(SIGIO_4_INC $ENV{SIGIO_INC4} CACHE STRING "SIGIO4 Include Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-sigio 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-sigio
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-sigio 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-sigio OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("sigio version is ${LIBVERSION}")
  set( SIGIO_LIBRARY ${PROJECT_BINARY_DIR}/lib/libsigio_${LIBVERSION}.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${SIGIO_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-sigio )
  else()
      set( CORE_BUILT ${SIGIO_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-sigio )
  endif()

endif()
set( SIGIO_LIBRARY_PATH ${SIGIO_LIBRARY} CACHE STRING "SIGIO Library Location" )
set( SIGIO_INCLUDE_PATH ${SIGIO_LIBRARY} CACHE STRING "SIGIO include Location" )

