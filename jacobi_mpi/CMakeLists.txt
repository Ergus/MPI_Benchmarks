# Copyright (C) 2021  Jimmy Aguilar Mena

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

set (CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set (CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

add_executable(jacobi_omp_mpi jacobi_omp_mpi.c)
target_compile_definitions(jacobi_omp_mpi PUBLIC -DISMATVEC=0)
target_link_libraries(jacobi_omp_mpi
  PUBLIC OpenMP::OpenMP_C
  PUBLIC MPI::MPI_C
  argparser)

add_test(NAME Test_jacobi
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./jacobi_omp_mpi 128 8 3)

# # Submitter scripts
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/utils/process.py" process.py COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/create_trace.py" create_paraver_trace.py COPYONLY)
