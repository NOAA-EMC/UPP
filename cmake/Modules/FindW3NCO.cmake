# This module looks for environment variables detailing where W3NCO lib is
# If variables are not set, W3NCO will be built from external source 
include(ExternalProject)
if(DEFINED ENV{W3NCO_LIB4} )
  set(W3NCO_LIBRARY $ENV{W3NCO_LIBd} CACHE STRING "W3NCO Library Location" )
  set(W3NCO_4_LIBRARY $ENV{W3NCO_LIB4} CACHE STRING "W3NCO_4 Library Location" )
  set(W3NCO_8_LIBRARY $ENV{W3NCO_LIB8} CACHE STRING "W3NCO_4 Library Location" )
else()  
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  set(SRC_DIR ${PROJECT_SOURCE_DIR}/../NCEPLIBS-w3nco)
  ExternalProject_Add(NCEPLIBS-w3nco 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-w3nco
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${SRC_DIR}
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${SRC_DIR} OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("w3nco version is ${LIBVERSION}")
  set( W3NCO_LIBRARY ${PROJECT_BINARY_DIR}/lib/libw3nco_${LIBVERSION}_d.a )
  set( W3NCO_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libw3nco_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${W3NCO_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-w3nco )
  else()
      set( CORE_BUILT ${W3NCO_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-w3nco )
  endif()

endif()
set( W3NCO_LIBRARY_PATH ${W3NCO_LIBRARY} CACHE STRING "W3NCO Library Location" )
set( W3NCO_4_LIBRARY_PATH ${W3NCO_4_LIBRARY} CACHE STRING "W3NCO_4 Library Location" )
set( W3NCO_8_LIBRARY_PATH ${W3NCO_8_LIBRARY} CACHE STRING "W3NCO_4 Library Location" )


