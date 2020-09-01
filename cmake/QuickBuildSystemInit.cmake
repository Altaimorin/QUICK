if(SHARED)
	message(STATUS "Build shared libs")
	set(BUILD_SHARED_LIBS TRUE)
endif()

if(DEBUG)
	message(STATUS "Build type is debug")
	set(CMAKE_BUILD_TYPE "Debug")
	add_definitions(-DDEBUG)
endif()

add_compile_options($<$<CONFIG:Debug>:-g> $<$<CONFIG:Debug>:-DDEBUG>)

set(OPT_FFLAGS "-O3")
set(OPT_CFLAGS "-O3")
set(OPT_CXXFLAGS "-O3")
set(CMAKE_Fortran_FLAGS_DEBUG "-O0")
set(CMAKE_C_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")

if(${COMPILER} STREQUAL GNU)
	add_definitions(-DGNU)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -mtune=native -ffree-form -cpp")  
elseif(${COMPILER} STREQUAL INTEL)
	set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ip -cpp -diag-disable 8291")
	set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} -traceback")
	set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -traceback")
	set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -traceback")
endif()

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR})
include(QuickCompilerConfig)
include(CopyTarget)
include(ConfigModuleDirs)

