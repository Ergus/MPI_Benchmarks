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

set(TEST "jacobi")

# CMAKE_C_FLAGS_RELEASE is a space separated string.
# But add_custom_command needs a list. In cmake lists are strings with ;
string(REGEX REPLACE "[ \t]+" ";" LIBBASE_FLAGS ${CMAKE_C_FLAGS_RELEASE})

# libbase.a
add_custom_command(
  OUTPUT libbase.o
  COMMAND ${ICC} -c ${CMAKE_CURRENT_SOURCE_DIR}/jacobi_base.c ${LIBBASE_FLAGS} -o libbase.o
  MAIN_DEPENDENCY ${CMAKE_CURRENT_SOURCE_DIR}/jacobi_base.c
  COMMENT "Building base object")

add_library(base STATIC libbase.o)

foreach(EXECUTABLE jacobi_omp_mpi jacobi_task_mpi jacobi_serial)

  message("- Add benchmark: ${EXECUTABLE}")

  # jacobi with parallel for
  add_executable(${EXECUTABLE} "${EXECUTABLE}.c")
  target_link_libraries(${EXECUTABLE} argparser base)

  add_test(NAME "Test_${EXECUTABLE}"
    COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${EXECUTABLE} 128 8 3)
endforeach()

# submiter
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_${TEST}_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_${TEST}_dim.sh @ONLY)

# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_extrae.sh" submit_extrae.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

# Process
configure_file("${PROJECT_SOURCE_DIR}/utils/process_mn.py" process_mn.py COPYONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/process_dim.py" process_dim.py COPYONLY)

