# This module looks for environment variables detailing where G2 lib is
# If variables are not set, G2 will be built from external source 
include(ExternalProject)
if(DEFINED ENV{G2_DIR} )
  message("HEY!! setting G2 library via environment variable")
  set(G2_LIBRARY $ENV{G2_LIBd} CACHE STRING "G2 Library Location" )
  set(G2_4_LIBRARY $ENV{G2_LIB4} CACHE STRING "G2_4 Library Location" )
  set(G2_INC4 $ENV{G2_INC4} CACHE STRING "G2_4 Include Location" )
  set(G2_INCd $ENV{G2_INCd} CACHE STRING "G2_d Include Location" )
  set(G2_LIBRARY_PATH ${G2_LIBRARY} CACHE STRING "G2 Library Location" )
  set(G2_4_LIBRARY_PATH ${G2_4_LIBRARY} CACHE STRING "G2_4 Library Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-g2 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-g2
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-g2 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-g2 OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("g2 version is ${LIBVERSION}")
  set( G2_LIBRARY ${PROJECT_BINARY_DIR}/lib/libg2_${LIBVERSION}_d.a )
  set( G2_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libg2_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${G2_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-g2 )
  else()
      set( CORE_BUILT ${G2_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-g2 )
  endif()
  
  set( G2_LIBRARY_PATH ${G2_LIBRARY} CACHE STRING "G2 Library Location" )
  set( G2_4_LIBRARY_PATH ${G2_4_LIBRARY} CACHE STRING "G2_4 Library Location" )
  
endif()
