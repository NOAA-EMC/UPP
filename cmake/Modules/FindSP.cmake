# This module looks for environment variables detailing where SP lib is
# If variables are not set, SP will be built from external source 

include(ExternalProject)
if(DEFINED ENV{SP_LIBd} )
  message("HEY!! setting SP library via environment variable")
  set(SP_LIBRARY $ENV{SP_LIBd} CACHE STRING "SP Library Location" )
  set(SP_4_LIBRARY $ENV{SP_LIB4} CACHE STRING "SP_4 Library Location" )
  set(SP_8_LIBRARY $ENV{SP_LIB4} CACHE STRING "SP_8 Library Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-sp 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-sp 
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-sp 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-sp OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("sp version is ${LIBVERSION}")
  set( SP_LIBRARY ${PROJECT_BINARY_DIR}/lib/libsp_${LIBVERSION}_d.a )
  set( SP_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libsp_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${SP_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-sp )
  else()
      set( CORE_BUILT ${SP_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-sp )
  endif()
endif()

set( SP_LIBRARY_PATH ${SP_LIBRARY} CACHE STRING "SP Library Location" )
set( SP_4_LIBRARY_PATH ${SP_4_LIBRARY} CACHE STRING "SP_4 Library Location" )

