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

set(TEST "cholesky")

# Add the executables and link
include_directories(${MKL_INCLUDE_DIR})

set(executable "cholesky_omp_mpi")

message("- Add benchmark: ${EXECUTABLE}")

add_executable(${executable} cholesky_omp_mpi.c)

add_test(NAME "Test_${executable}"
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${executable} 512 16 0)

# submiter
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_dim.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_extrae.sh" submit_extrae.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)
