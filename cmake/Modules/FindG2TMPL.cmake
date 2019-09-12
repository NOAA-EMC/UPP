# This module looks for environment variables detailing where G2TMPL lib is
# If variables are not set, G2TMPL will be built from external source 
include(ExternalProject)
if(DEFINED ENV{G2TMPL_DIR} )
  message("HEY!! setting G2TMPL library via environment variable")
  set(G2TMPL_LIBRARY $ENV{G2TMPL_LIBd} CACHE STRING "G2TMPL Library Location" )
  set(G2TMPL_INC $ENV{G2TMPL_INCd} CACHE STRING "G2TMPL_d Include Location" )
else()
  set(CMAKE_INSTALL_PREFIX ${PROJECT_BINARY_DIR})
  ExternalProject_Add(NCEPLIBS-g2tmpl 
    PREFIX ${PROJECT_BINARY_DIR}/NCEPLIBS-g2tmpl
    CMAKE_ARGS
      -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
      -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
      -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
      -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      -DCMAKE_BUILD_TYPE=RELEASE
    SOURCE_DIR ${PROJECT_SOURCE_DIR}/NCEPLIBS-g2tmpl 
    INSTALL_DIR ${CMAKE_INSTALL_PREFIX}
    BUILD_COMMAND make
    INSTALL_COMMAND make install
  )
  execute_process(COMMAND grep "set(VERSION" CMakeLists.txt WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}/NCEPLIBS-g2tmpl OUTPUT_VARIABLE LIBVERSION)
  string(REPLACE "set(VERSION " "" LIBVERSION ${LIBVERSION})
  string(REPLACE ")" "" LIBVERSION ${LIBVERSION})
  string(REPLACE "\n" "" LIBVERSION ${LIBVERSION})
  message("g2tmpl version is ${LIBVERSION}")
  set( G2TMPL_LIBRARY ${PROJECT_BINARY_DIR}/lib/libg2tmpl_${LIBVERSION}_d.a )
  set( G2TMPL_4_LIBRARY ${PROJECT_BINARY_DIR}/lib/libg2tmpl_${LIBVERSION}_4.a )
  if( CORE_BUILT )
      list( APPEND CORE_BUILT ${G2TMPL_LIBRARY} )
      list( APPEND EXT_BUILT NCEPLIBS-g2tmpl )
  else()
      set( CORE_BUILT ${G2TMPL_LIBRARY} )
      set( EXT_BUILT NCEPLIBS-g2tmpl )
  endif()
  
  set( G2TMPL_LIBRARY_PATH ${G2TMPL_LIBRARY} CACHE STRING "G2TMPL Library Location" )
  set( G2TMPL_INCLUDE_PATH ${G2TMPL_INC} CACHE STRING "G2TMPL Include Location" )
  
endif()
