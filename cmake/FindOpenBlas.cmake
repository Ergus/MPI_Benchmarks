###############################################################################
# Find OpenBlas
#
# This sets the following variables:
# OPENBLAS_FOUND        - True if OpenBLAS was found.
# OPENBLAS_INCLUDES     - Directories containing the OpenBLAS include files.
# OPENBLAS_LIBRARIES    - Libraries needed to use OpenBLAS.
###############################################################################

INCLUDE(FindPackageHandleStandardArgs)

MESSAGE(STATUS "CMAKE_OPENBLAS_ROOT: ${CMAKE_OPENBLAS_ROOT}")
MESSAGE(STATUS "OPENBLAS_ROOT: $ENV{OPENBLAS_ROOT}")

# Set extended paths
SET(INC_PATHS
  /usr/include
  /usr/local/include
  /usr/include )

SET(LIB_PATHS
  /usr/lib
  /usr/lib64
  /usr/local/lib
  /usr/lib/x86_64-linux-gnu )

# Find where are the header files.
FIND_PATH(OPENBLAS_INCLUDE
  NAMES cblas.h
  HINTS ${CMAKE_OPENBLAS_ROOT} ENV OPENBLAS_ROOT
  PATHS ${INC_PATHS} ENV CPLUS_INCLUDE_PATH
  DOC "OpenBlas include path"
  PATH_SUFFIXES include
  )

FIND_LIBRARY(OPENBLAS_LIBRARY
  NAMES openblas
  HINTS ${CMAKE_OPENBLAS_ROOT} ENV OPENBLAS_ROOT
  PATHS ${LIB_PATHS} ENV LD_LIBRARY_PATH
  DOC "OpenBlas library"
  PATH_SUFFIXES lib lib64 lib32)

FIND_PACKAGE_HANDLE_STANDARD_ARGS(OpenBlas
  REQUIRED_VARS OPENBLAS_INCLUDE OPENBLAS_LIBRARY)

MARK_AS_ADVANCED(OPENBLAS_INCLUDES OPENBLAS_LIBRARIES)

SET(OPENBLAS_INCLUDES  ${OPENBLAS_INCLUDE})
SET(OPENBLAS_LIBRARIES ${OPENBLAS_LIBRARY})
