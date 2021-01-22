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


add_executable(matvec_omp_mpi matvec_omp_mpi.cpp matvec.c)
target_link_libraries(matvec_strong argparser)

# # Run test
# ADD_TEST(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 64 8 1 1 matrix)

# # python test
# IF(PYTHONINTERP_FOUND)
#   CONFIGURE_FILE(prove.py prove.py)
#   ADD_TEST(NAME Test_${execut}_results COMMAND
#     python3 prove.py matrix_A.mat matrix_B.mat matrix_C.mat)
# ELSE()
#   MESSAGE("Python interpreter not found")
# ENDIF(PYTHONINTERP_FOUND)

# # Submitter scripts
# #CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/matmul_mpi/submit.sh submit.sh COPYONLY)
# #CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/matmul_mpi/profile-matmul.sh profile-matmul.sh COPYONLY)
# #CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/matmul_mpi/gnuplot.plt gnuplot.plt COPYONLY)

# configure_file(run_matvec.sh run_matvec.sh @ONLY)
# configure_file(profile_matvec.sh profile_matvec.sh @ONLY)
# configure_file(extrae_template.xml extrae_template.xml COPYONLY)
