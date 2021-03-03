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

# Run test

# python test
# if (Python_FOUND AND Python_Interpreter_FOUND AND Python_NumPy_FOUND)
#   CONFIGURE_FILE(prove.py prove.py)
#   ADD_TEST(NAME Test_matmul_results
#     COMMAND ${Python_EXECUTABLE} prove.py matmul_A.mat matmul_B.mat matmul_C.mat)
#   ADD_TEST(NAME Test_matvec_results
#     COMMAND ${Python_EXECUTABLE} prove.py matvec_A.mat matvec_B.mat matvec_C.mat)
# else ()
#   MESSAGE("Python interpreter not found")
# endif () # PYTHONINTERP_FOUND

# # Submitter scripts
CONFIGURE_FILE(submiter_mn.sh submiter_mn.sh @ONLY)
CONFIGURE_FILE(submit_mn.sh submit_mn.sh @ONLY)