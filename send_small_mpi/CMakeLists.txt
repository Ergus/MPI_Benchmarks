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

add_executable(normal_mpi normal_mpi.c)
target_link_libraries(normal_mpi
  PUBLIC MPI::MPI_C
  argparser)
add_test(NAME Test_normal
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./normal_mpi 128)

add_executable(type_mpi type_mpi.c)
target_link_libraries(type_mpi
  PUBLIC MPI::MPI_C
  argparser)
add_test(NAME Test_type
  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./type_mpi 128)

