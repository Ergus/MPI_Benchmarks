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
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ./matmul_omp_mpi 128 1 1 1)

add_executable(matvec_omp_mpi matmul_omp_mpi.c)
target_compile_definitions(matvec_omp_mpi PUBLIC -DISMATVEC=1)
target_link_libraries(matvec_omp_mpi
  PUBLIC OpenMP::OpenMP_C
  PUBLIC MPI::MPI_C
  argparser)
add_test(NAME Test_matvec
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ./matvec_omp_mpi 128 1 1 1)

# Run test

# # Submitter scripts
# CONFIGURE_FILE(submit.sh submit.sh COPYONLY)
# CONFIGURE_FILE(profile-matmul.sh profile-matmul.sh COPYONLY)
# CONFIGURE_FILE(gnuplot.plt gnuplot.plt COPYONLY)

# python test
IF(PYTHONINTERP_FOUND)
  CONFIGURE_FILE(prove.py prove.py)
  ADD_TEST(NAME Test_matmul_results
    COMMAND python3 prove.py matmul_A.mat matmul_B.mat matmul_C.mat)
  ADD_TEST(NAME Test_matvec_results
    COMMAND python3 prove.py matvec_A.mat matvec_B.mat matvec_C.mat)
ELSE()
  MESSAGE("Python interpreter not found")
ENDIF() # PYTHONINTERP_FOUND
