cmake_minimum_required(VERSION 3.15)
project(ncep_post VERSION 1.1.1)
set(${PROJECT_NAME}_VERSION ${PROJECT_VERSION} CACHE INTERNAL "${PROJECT_NAME} version number")

enable_language (Fortran)

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE RELEASE CACHE STRING
    "Choose the type of build, options are: PRODUCTION Debug Release."
    FORCE)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
  set(IntelComp true )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU*" OR CMAKE_CXX_COMPILER_ID MATCHES "Clang*")
  set(GNUComp true )
elseif(CMAKE_CXX_COMPILER_ID MATCHES "pgc*")
  set(PGIComp true )
endif()

STRING(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "RelWithDebInfo" BUILD_RELEASE)
if(NOT BUILD_RELEASE )
  STRING(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "RELEASE" BUILD_RELEASE)
endif()
STRING(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "PRODUCTION" BUILD_PRODUCTION)
STRING(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "DEBUG" BUILD_DEBUG)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()

find_package(MPI REQUIRED)
find_package(OpenMP)

find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(HDF5 COMPONENTS C Fortran HL REQUIRED)
find_package(NetCDF REQUIRED)
find_package(Jasper REQUIRED)

if(NOT TARGET sigio)
  find_package(sigio REQUIRED)
endif()

if(NOT TARGET w3nco_4)
  find_package(w3nco REQUIRED)
endif()

if(NOT TARGET g2_4)
  find_package(g2 REQUIRED)
endif()

if(NOT TARGET g2tmpl_d)
  find_package(g2tmpl REQUIRED)
endif()

if(NOT TARGET bacio_4)
  find_package(bacio REQUIRED)
endif()

if(NOT TARGET ip_4)
  find_package(ip REQUIRED)
endif()

if(NOT TARGET sp_4)
  find_package(sp REQUIRED)
endif()

if(NOT TARGET sfcio)
  find_package(sfcio REQUIRED)
endif()

if(NOT TARGET nemsio)  
  find_package(nemsio REQUIRED)
endif()

if(NOT TARGET gfsio)
  find_package(gfsio REQUIRED)
endif()

if(NOT TARGET w3emc_4)
  find_package(w3emc REQUIRED)
endif()

if(NOT TARGET crtm)
  find_package(CRTM REQUIRED)
endif()


add_subdirectory(sorc/ncep_post.fd)
