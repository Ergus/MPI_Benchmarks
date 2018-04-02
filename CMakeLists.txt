cmake_minimum_required (VERSION 3.0)

project(mpi_benchmarks)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CMAKE_VERBOSE_MAKEFILE ON)

enable_testing()
set(Python_ADDITIONAL_VERSIONS 3.7 3.6 3.5 3.4 3.3 3.2 3.1)
find_package(PythonInterp MODULE)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  set(OMPFLAG "-qopenmp")
  set(FLAGS "-link_mpi=opt_mt -O2")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  set(OMPFLAG "-fopenmp")
  set(FLAGS "-O2 -lrt")
else()
  message(FATAL_ERROR "Unknown compiler: ${CMAKE_CXX_COMPILER_ID}")
endif ()

set(CMAKE_C_COMPILER mpicc)
set(CMAKE_CXX_COMPILER mpicxx)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${FLAGS} ${OMPFLAG}")
set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS} ${FLAGS} ${OMPFLAG}")

add_definitions(-DEXECUTABLE_TYPE="mpi")

file(GLOB dirlist ${PROJECT_SOURCE_DIR}/*_mpi*)

include_directories(${PROJECT_SOURCE_DIR}/common)
add_subdirectory(${PROJECT_SOURCE_DIR}/common)

foreach (thedir ${dirlist})
  add_subdirectory(${thedir})
endforeach ()

set(CMAKE_INSTALL_PREFIX ${PROJECT_SOURCE_DIR}/build/install)

configure_file(extrae_template.xml extrae_template.xml COPYONLY)
