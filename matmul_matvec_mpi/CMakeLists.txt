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

add_executable(matmul_omp_mpi matmul_omp_mpi.c)
target_compile_definitions(matmul_omp_mpi PUBLIC -DISMATVEC=0)
target_link_libraries(matmul_omp_mpi
  PUBLIC OpenMP::OpenMP_C
  PUBLIC MPI::MPI_C
  argparser)
add_test(NAME Test_matmul
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./matmul_omp_mpi 128 8 3)

add_executable(matvec_omp_mpi matmul_omp_mpi.c)
target_compile_definitions(matvec_omp_mpi PUBLIC -DISMATVEC=1)
target_link_libraries(matvec_omp_mpi
  PUBLIC OpenMP::OpenMP_C
  PUBLIC MPI::MPI_C
  argparser)
add_test(NAME Test_matvec
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./matvec_omp_mpi 128 16 1)

# Run test

# python test
configure_file(prove.py prove.py COPYONLY)

if (Python_FOUND AND Python_Interpreter_FOUND AND Python_NumPy_FOUND)
  CONFIGURE_FILE(prove.py prove.py)
  ADD_TEST(NAME Test_matmul_results
    COMMAND ${Python_EXECUTABLE} prove.py matmul_A.mat matmul_B.mat matmul_C.mat)
  ADD_TEST(NAME Test_matvec_results
    COMMAND ${Python_EXECUTABLE} prove.py matvec_A.mat matvec_B.mat matvec_C.mat)
else ()
  MESSAGE("Python interpreter not found")
endif () # PYTHONINTERP_FOUND

# Submitter scripts
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)
CONFIGURE_FILE("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

configure_file("${PROJECT_SOURCE_DIR}/utils/process.py" process.py COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/create_trace.py" create_paraver_trace.py COPYONLY)

