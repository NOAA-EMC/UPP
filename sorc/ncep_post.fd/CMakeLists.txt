cmake_minimum_required(VERSION 2.6)

set(EXENAME ncep_post )
file(GLOB EXE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.f ${CMAKE_CURRENT_SOURCE_DIR}/*.f90 ${CMAKE_CURRENT_SOURCE_DIR}/*.F)
file(GLOB EXE_C_SRC ${CMAKE_CURRENT_SOURCE_DIR}/*.c)
list( REMOVE_ITEM EXE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/TRPAUS_NAM.f)
list( REMOVE_ITEM EXE_SRC ${CMAKE_CURRENT_SOURCE_DIR}/WRF_STUBS.f)
if(IntelComp)
    message("setting intel flags")
    set(CMAKE_Fortran_8_FLAGS " -free -O3 -convert big_endian -traceback -g -fp-model source ${OpenMP_Fortran_FLAGS}" CACHE INTERNAL "" )
    set(CMAKE_C_FLAGS "-DLINUX -Dfunder -DFortranByte=char -DFortranInt=int -DFortranLlong='long long'" CACHE INTERNAL "")
else()
    message("unkown compiler!")
    exit()
endif()

set_source_files_properties( ${EXE_C_SRC} PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS} )
add_executable(${EXENAME} ${EXE_SRC} ${EXE_C_SRC} )
set_source_files_properties( ${EXE_SRC} PROPERTIES COMPILE_FLAGS ${CMAKE_Fortran_8_FLAGS} )
set_source_files_properties( ${EXE_C_SRC} PROPERTIES COMPILE_FLAGS ${CMAKE_C_FLAGS} )
set_target_properties( ${EXENAME} PROPERTIES LINK_FLAGS -qopenmp ) 
message("hey, g2 inc is ${G2TMPL_INCd}")
include_directories( ${MPI_Fortran_INCLUDE_PATH} ${NETCDF_INCLUDES} ${SIGIO_INC} ${NEMSIO_INC} ${SFCIO_INC} ${W3EMC_INC4} ${G2_INC4} ${G2TMPL_INCd} ${GFSIO_INC} ${CRTM_INC} )
target_link_libraries(${EXENAME} ${BACIO_LIB4} ${SP_LIB4} ${W3NCO_LIB4} ${SIGIO_LIB} ${NEMSIO_LIB} ${IP_LIB4} ${SFCIO_LIB} ${W3EMC_LIB4} ${ZLIB} ${G2TMPL_LIBd} ${G2_LIB4}  ${PNG_LIBRARIES} ${JASPER_LIBRARIES} ${GFSIO_LIB4} ${CRTM_LIB} ${MPI_Fortran_LIBRARIES} ${NETCDF_LIBRARIES_F90} ${NETCDF_LIBRARIES} ${BACIO_LIB4}  )


install(TARGETS ${EXENAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
        LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
        ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib) 
install(DIRECTORY ${CMAKE_INCLUDE_OUTPUT_DIRECTORY}/ DESTINATION include_4 )

