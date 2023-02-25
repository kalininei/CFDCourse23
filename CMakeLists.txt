cmake_minimum_required(VERSION 3.0)

Project(cfdapp)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)

# build all binaries to build/bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ------ boost
include_directories(${CMAKE_SOURCE_DIR}/..)

# ------ amgcl
FetchContent_Declare(
	amgcl
	GIT_REPOSITORY https://github.com/ddemidov/amgcl
	GIT_TAG 1.4.2
)
FetchContent_GetProperties(amgcl)
if(NOT amgcl_POPULATED)
	FetchContent_Populate(amgcl)
endif()
include_directories(${amgcl_SOURCE_DIR})

# ------ cfdlib
include_directories(${CMAKE_SOURCE_DIR}/lib)
set(CFDLIB_TARGET cfdlib)

# default input directory to "input_data/"
add_definitions(-DCFDLIB_INPUT_DIR_PATH="${CMAKE_SOURCE_DIR}/input_data/")

# default output directory
if (MSVC)
	# build/output for visual studio
	set(OUTPUT_DIR_PATH ${CMAKE_BINARY_DIR}/output/)
	add_definitions(-DCFDLIB_OUTPUT_DIR_PATH="${OUTPUT_DIR_PATH}")
	file(MAKE_DIRECTORY ${OUTPUT_DIR_PATH})
else(MSVC)
	# ./ for linux
	add_definitions(-DCFDLIB_OUTPUT_DIR_PATH="")
endif(MSVC)

add_subdirectory(lib)
add_subdirectory(progs)
add_subdirectory(doc)
