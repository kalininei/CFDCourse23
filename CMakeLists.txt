cmake_minimum_required(VERSION 3.0)

Project(cfdapp)
include(FetchContent)

set(CMAKE_CXX_STANDARD 17)

# build all binaries to build/bin
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# ------ boost
SET(BOOST_MODULES
	"libs/core;libs/static_assert;libs/type_traits;libs/property_tree;libs/config;libs/optional;libs/throw_exception;"
	"libs/assert;libs/move;libs/utility;libs/any;libs/type_index;libs/container_hash;libs/mpl;libs/preprocessor;libs/multi_index;"
	"libs/serialization;libs/tuple;libs/bind;libs/foreach;libs/iterator;libs/detail;libs/math;libs/geometry;"
	"libs/multiprecision;libs/range;libs/array;libs/lexical_cast;libs/variant;libs/algorithm;"
	"libs/function;")
FetchContent_Declare(
	boost
	GIT_REPOSITORY https://github.com/boostorg/boost.git
	GIT_SUBMODULES "tools/cmake;" ${BOOST_MODULES}
)
FetchContent_GetProperties(boost)
if(NOT boost_POPULATED)
	FetchContent_Populate(boost)
endif()
set(BOOST_INCLUDE)
foreach (m in ${BOOST_MODULES})
	LIST(APPEND BOOST_INCLUDE ${boost_SOURCE_DIR}/${m}/include)
endforeach()
message(${BOOST_INCLUDE})
include_directories(${BOOST_INCLUDE})

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
