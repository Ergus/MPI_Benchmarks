# Copyright (C) 2021-2022  Jimmy Aguilar Mena

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

# 3.9 required to use FindOpenMP
# 3.12.4 required to use FindPython
cmake_minimum_required(VERSION 3.10)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(CCOMPILER gnu CACHE STRING "C compiler to use gnu or intel")
set(CCOMPILERVALUES gnu intel)
set_property(CACHE CCOMPILER PROPERTY STRINGS ${CCOMPILERVALUES})

list(FIND CCOMPILERVALUES ${CCOMPILER} INDEX)
if(INDEX EQUAL -1)
  message(FATAL_ERROR "Invalid CCOMPILER value: ${CCOMPILER}")
endif()

if ("${CCOMPILER}" STREQUAL "intel")
  find_package(INTEL REQUIRED)
  add_compile_options(--Wn,-diag-disable=858)
endif ()

project(mpi_benchmarks)

option(WITH_ASAN "Build tests with address sanitizer" false)
option(WITH_EXTRAE "Build with Extrae traces support" false)
option(WITH_SERIAL "Build serial version of the code." false)
option(WITH_MKL "Use MKL as BLAS library." false)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)

add_definitions(_GNU_SOURCE)

# Enable address sanitizer
if (WITH_ASAN)
  message("-- Using ASAN")
  add_compile_options(-fno-omit-frame-pointer -fsanitize=address)

  find_library(LIBASAN
    NAMES asan
    HINTS ENV LIBRARY_PATH LD_LIBRARY_PATH)

  set (ASAN_OPTIONS
    "\"detect_leaks=0:detect_stack_use_after_return=1:check_initialization_order=1:strict_init_order=1\"")

  set (CMAKE_EXE_LINKER_FLAGS
    "${CMAKE_EXE_LINKER_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
endif ()

# Enable Extrae
if (WITH_EXTRAE)
  find_package(Extrae REQUIRED COMPONENTS ompitrace)
  if (Extrae_FOUND)
    add_compile_definitions(__WITH_EXTRAE=1)
    include_directories(${EXTRAE_INCLUDES})
    link_libraries(${EXTRAE_LIBRARIES})
  endif ()
endif()

if (CMAKE_BUILD_TYPE MATCHES ".*Debug.*")
  set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
endif()

find_package(MPI REQUIRED)
link_libraries(MPI::MPI_C)

find_package(OpenMP REQUIRED COMPONENTS C)
link_libraries(OpenMP::OpenMP_C)

find_package(Python MODULE COMPONENTS Interpreter NumPy)

# Lapack - MKL
if (WITH_MKL)
  set(BLA_VENDOR Intel10_64lp_seq)
endif ()

find_package(LAPACK REQUIRED)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
link_libraries(${LAPACK_LIBRARIES})

# add the subdirectories here
add_subdirectory(utils/ArgParserC)
link_libraries(argparser)

include_directories(utils)
enable_testing()

configure_file(${PROJECT_SOURCE_DIR}/utils/ArgParserBash/argparse.sh
  ${PROJECT_BINARY_DIR}/argparse.sh COPYONLY)

# Lists of include files
file(GLOB children RELATIVE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/*_mpi)
# Include all subdirectories with _ompss2 suffix.
foreach (child ${children})
  message ("-- Adding subdirectory?: ${child}")
  if ((IS_DIRECTORY ${CMAKE_CURRENT_LIST_DIR}/${child})
      AND (EXISTS ${CMAKE_CURRENT_LIST_DIR}/${child}/CMakeLists.txt))
    message ("-- Adding subdirectory: ${child}")
    add_subdirectory(${child})
  endif ()
endforeach ()
