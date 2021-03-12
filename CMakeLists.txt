# This is the CMake build file for EMC_post (also known as UPP).
#
cmake_minimum_required(VERSION 3.15)

file(STRINGS "VERSION" pVersion LIMIT_COUNT 1)

project(
  upp
  VERSION ${pVersion}
  LANGUAGES Fortran)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules/Modules")

# Handle user options.
option(OPENMP "use OpenMP threading" OFF)
option(BUILD_POSTEXEC "Build NCEPpost executable" ON)
option(BUILD_WITH_WRFIO "Build NCEPpost with WRF-IO library" OFF)
option(ENABLE_DOCS "Enable generation of doxygen-based documentation." OFF)

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
find_package(NetCDF REQUIRED COMPONENTS Fortran)

if(OPENMP)
  find_package(OpenMP REQUIRED COMPONENTS Fortran)
endif()

find_package(bacio REQUIRED)
find_package(crtm REQUIRED)
find_package(g2 REQUIRED)
find_package(g2tmpl REQUIRED)

if(BUILD_POSTEXEC)
  find_package(gfsio REQUIRED)
  find_package(ip REQUIRED)
  find_package(nemsio REQUIRED)
  find_package(sfcio REQUIRED)
  find_package(sigio REQUIRED)
  find_package(sp REQUIRED)
  find_package(w3nco REQUIRED)
  if(BUILD_WITH_WRFIO)
    find_package(wrf_io REQUIRED)
  endif()
endif()

add_subdirectory(sorc)
add_subdirectory(parm)

# If desired, build the doxygen docs.
if(ENABLE_DOCS)
  find_package(Doxygen REQUIRED)
  add_subdirectory(docs)
endif()

