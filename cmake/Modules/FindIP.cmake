# This module looks for environment variables detailing where IP lib is
# If variables are not set, IP will be built from external source 

include(ExternalProject)
if(DEFINED ENV{IP_LIBd} )
  message("HEY!! setting IP library via environment variable")
  set(IP_LIBRARY $ENV{IP_LIBd} CACHE STRING "IP Library Location" )
  set(IP_4_LIBRARY $ENV{IP_LIB4} CACHE STRING "IP_4 Library Location" )
  set(IP_8_LIBRARY $ENV{IP_LIB4} CACHE STRING "IP_8 Library Location" )
  set(IP_INC4 $ENV{IP_INC4} CACHE STRING "IP_4 Include Location" )
  set(IP_INC8 $ENV{IP_INC8} CACHE STRING "IP_8 Include Location" )
  set(IP_INCd $ENV{IP_INCd} CACHE STRING "IP_8 Include Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-ip 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-ip 
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-ip 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-ip OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("ip version is ${LIBVERSION}")
  set( IP_LIBRARY ${PROJECT_BINARY_DIR}/lib/libip_${LIBVERSION}_d.a )
  set( IP_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libip_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${IP_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-ip )
  else()
      set( CORE_BUILT ${IP_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-ip )
  endif()
endif()

set( IP_LIBRARY_PATH ${IP_LIBRARY} CACHE STRING "IP Library Location" )
set( IP_4_LIBRARY_PATH ${IP_4_LIBRARY} CACHE STRING "IP_4 Library Location" )

