# Compilation


# Copyright (C) 2019  Jimmy Aguilar Mena

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

set(BLA_VENDOR Intel10_64ilp)
find_package(BLAS REQUIRED)

find_path(
  MKL_INCLUDE_DIR mkl.h
  HINTS ENV MKLROOT
  PATH_SUFFIXES include)

# Add the executables and link
include_directories(${MKL_INCLUDE_DIR})

add_executable(cholesky cholesky.c util.c)
target_link_options(cholesky PUBLIC ${BLAS_LINKER_FLAGS})
target_link_libraries(cholesky
  PUBLIC OpenMP::OpenMP_C
  PUBLIC MPI::MPI_C
  ${BLAS_LIBRARIES}
  argparser)
