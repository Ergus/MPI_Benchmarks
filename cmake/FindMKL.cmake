# - Find Intel MKL
# Find the MKL libraries
#
# Options:
#
#   MKL_STATIC        :   use static linking
#   MKL_MULTI_THREADED:   use multi-threading
#   CMAKE_OMP_LIB   :     omp library to use (imp5 | gomp)
#   CMAKE_MPI_LIB   :     mpi library to use (openmpi | impi | mpich2)
#   MKL_WITH_SCALAPACK:   use scalapack?
#
# This module defines the following variables:
#
#   MKL_FOUND        : True if MKL_INCLUDE_DIR are found
#   MKL_INCLUDES     : set when MKL_INCLUDE_DIR found
#   MKL_LIBRARIES    : the library to link against.
#   MKL_DEFINITIONS  : compiler definitions
#   MKL_EXTRA_LINKS  : extra libraries needed in link time

#  The commands are here
#  https://software.intel.com/en-us/articles/intel-mkl-link-line-advisor

include(FindPackageHandleStandardArgs)

# GEt environment MKLROOT and test it
set(MKLROOT $ENV{MKLROOT})

set(MKL_VALID_OMP iomp5 gomp)
set(MKL_VALID_MPI openmpi impi mpich2)

if (NOT MKLROOT)
  message(FATAL_ERROR "Environment variable MKLROOT is empty. Load MKL")
elseif (NOT EXISTS ${MKLROOT}) #if set, test it
  message(FATAL_ERROR "Directory ${MKLROOT} doesn't exist")
endif()

message(STATUS "MKLROOT: ${MKLROOT}")
message(STATUS "MKL_STATIC: ${MKL_STATIC}")
message(STATUS "MKL_MULTI_THREADED: ${MKL_MULTI_THREADED}")
message(STATUS "MKL_WITH_SCALAPACK: ${MKL_WITH_SCALAPACK}")
message(STATUS "CMAKE_OMP_LIB: ${CMAKE_OMP_LIB}")
message(STATUS "CMAKE_MPI_LIB: ${CMAKE_MPI_LIB}")

# Find headers
find_path(MKL_INCLUDE_DIR mkl.h ${MKLROOT}/include)

# Find libraries
# Save suffix temporally to set desired
set(_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES ${CMAKE_FIND_LIBRARY_SUFFIXES})

if (MKL_STATIC)  # Set the static or dynamic
  set(CMAKE_FIND_LIBRARY_SUFFIXES .a)
else ()
  set(CMAKE_FIND_LIBRARY_SUFFIXES .so)
endif ()

# Libraries according to architecture
set(MKL_LIBS_LIST mkl_core)
set(MKL_LIBS_FOUND "")
set(MKL_EXTRA_LINKS "-lpthread -lm -ldl")

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  set(MKL_EXTRA_LINKS "-Wl,--no-as-needed ${MKL_EXTRA_LINKS}")
endif ()

#Threads
if (MKL_MULTI_THREADED OR MKL_WITH_SCALAPACK)  # Search the omp library if multi-threaded
  # Find lib for openmp
  if (CMAKE_OMP_LIB IN_LIST MKL_VALID_OMP) # if declared test it
	set(OMP_NAME ${CMAKE_OMP_LIB})
  elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	set(OMP_NAME gomp)
  else ()
	set(OMP_NAME iomp5)
  endif ()
  message(STATUS "OMP_NAME: ${OMP_NAME}")

  # After searched or set, test and append properly
  list(APPEND MKL_EXTRA_LINKS "-l${OMP_NAME}") # adds -liomp5 or -lgomp
  if (OMP_NAME STREQUAL "gomp")
	list(APPEND MKL_LIBS_LIST mkl_gnu_thread)
  elseif (OMP_NAME STREQUAL "iomp5")
	list(APPEND MKL_LIBS_LIST mkl_intel_thread)
  else ()
	message(FATAL_ERROR "The OMP library set is wrong for MKL")
  endif ()
else ()
  list(APPEND MKL_LIBS_LIST mkl_sequential)
endif ()

# Set parameters dependent of the architecture
if (CMAKE_SIZEOF_VOID_P EQUAL 8) # 64 bit
  set(MKL_LIB_PATH "${MKLROOT}/lib/intel64")
  list(APPEND MKL_LIBS_LIST mkl_intel_ilp64)
  set(MKL_DEFINITIONS "-DMKL_ILP64 -m64")
  if (MKL_WITH_SCALAPACK) # Scalapack only makes sense in 64 bits
	list(APPEND MKL_LIBS_LIST mkl_scalapack_lp64)

	if (CMAKE_MPI_LIB)
	  if (CMAKE_MPI_LIB IN_LIST MKL_VALID_MPI)
		set(MPI_NAME ${CMAKE_MPI_LIB})
	  else ()
		message(FATAL_ERROR "The MPI library set is unknown for MKL")
	  endif ()
	elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
	  set(MPI_NAME openmpi)
	elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
	  set(MPI_NAME iomp5)
	else()
	  set(MPI_NAME mpich2)
	endif ()

	if (MPI_NAME STREQUAL "openmpi")
	  list(APPEND MKL_LIBS_LIST mkl_blacs_openmpi_ilp64)
	else ()
	  list(APPEND MKL_LIBS_LIST mkl_blacs_intelmpi_ilp64)
	endif ()
  endif ()
else () # 32 bits
  set(MKL_LIB_PATH "${MKLROOT}/lib/ia32")
  list(APPEND MKL_LIBS_LIST mkl_intel)
  set(MKL_DEFINITIONS "-m32")
endif ()

# Check the validity of lib path
if (NOT EXISTS ${MKL_LIB_PATH})
  message(FATAL_ERROR "The MKL_LIB_PATH: ${MKL_LIB_PATH} does not exist")
endif ()

# Search all the libraries in the list
foreach (thelib ${MKL_LIBS_LIST})
  SET(FOUND_LIB "FOUND_LIB-NOTFOUND")
  find_library(FOUND_LIB ${thelib} ${MKL_LIB_PATH})
  set(${thelib} ${FOUND_LIB})
  list(APPEND MKL_LIBS_FOUND ${FOUND_LIB})
endforeach ()

set(CMAKE_FIND_LIBRARY_SUFFIXES ${_MKL_ORIG_CMAKE_FIND_LIBRARY_SUFFIXES})

find_package_handle_standard_args(MKL DEFAULT_MSG
  MKL_INCLUDE_DIR ${MKL_LIBS_LIST})

if(MKL_FOUND)
  set(MKL_INCLUDES ${MKL_INCLUDE_DIR})
  set(MKL_LIBRARIES "${MKL_LIBS_FOUND}")
  message(STATUS "MKL_LIBRARIES: ${MKL_LIBRARIES}")
endif()
