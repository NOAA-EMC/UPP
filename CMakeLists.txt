cmake_minimum_required(VERSION 3.15)

file(STRINGS "VERSION" pVersion)

project(
  ncep_post
  VERSION ${pVersion}
  LANGUAGES Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/")

if(POLICY CMP0054)
  cmake_policy(SET CMP0054 NEW)
endif()
if(POLICY CMP0012)
  cmake_policy(SET CMP0012 NEW)
endif()

option(OPENMP "use OpenMP threading" OFF)

if(NOT CMAKE_BUILD_TYPE MATCHES "^(Debug|Release|RelWithDebInfo|MinSizeRel)$")
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE
      "Release"
      CACHE STRING "Choose the type of build." FORCE)
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
                                               "MinSizeRel" "RelWithDebInfo")
endif()

if(NOT CMAKE_Fortran_COMPILER_ID MATCHES "^(Intel|GNU|Clang|AppleClang)$")
  message(WARNING "Compiler not officially supported: ${CMAKE_Fortran_COMPILER_ID}")
endif()

find_package(MPI REQUIRED)
find_package(PNG REQUIRED)
find_package(ZLIB REQUIRED)
find_package(Jasper REQUIRED)
find_package(HDF5 COMPONENTS Fortran HL REQUIRED)
find_package(NetCDF REQUIRED)

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

set(CMAKE_FIND_PACKAGE_PREFER_CONFIG true)

if(NOT TARGET sigio_4)
  find_package(sigio REQUIRED)
endif()

if(NOT TARGET w3nco_4)
  find_package(w3nco REQUIRED)
endif()

if(NOT TARGET g2_4)
  find_package(g2 REQUIRED)
endif()

if(NOT TARGET g2tmpl)
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

if(NOT TARGET sfcio_4)
  find_package(sfcio REQUIRED)
endif()

if(NOT TARGET nemsio)
  find_package(nemsio REQUIRED)
endif()

if(NOT TARGET gfsio_4)
  find_package(gfsio REQUIRED)
endif()

if(NOT TARGET w3emc_4)
  find_package(w3emc REQUIRED)
endif()

if(NOT TARGET crtm)
  find_package(crtm REQUIRED)
endif()

add_subdirectory(sorc)

install(FILES parm/postxconfig-NT-GFS.txt parm/postxconfig-NT-GFS-F00.txt DESTINATION share)
