# - Find GEMPAK 
if(DEFINED ENV{GEMPAK})
    set(GEMINC $ENV{GEMPAK}/include CACHE STRING "GEMPAK Include Path" )
    set(OS_INC $ENV{OS_INC} CACHE STRING "GEMPAK OS Include Path" )
else ()
    message(FATAL_ERROR "Could not find GEMPAK module")
endif()
