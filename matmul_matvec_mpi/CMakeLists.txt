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

if (Python_FOUND AND Python_Interpreter_FOUND AND Python_NumPy_FOUND)
  configure_file(validate.py validate.py @ONLY)
else ()
  message("Python interpreter not found")
endif () # PYTHONINTERP_FOUND


set (TESTLIST matmul matvec)
foreach (ISMATVEC 0;1)
  list(GET TESTLIST ${ISMATVEC} TEST)

  set(executable "${TEST}_omp_mpi")
  add_executable(${executable} matmul_omp_mpi.c)
  target_compile_definitions(${executable} PUBLIC -DISMATVEC=${ISMATVEC})
  target_link_libraries(${executable} argparser)

  add_test(NAME Test_${executable}
    COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 4 ./${executable} 256 16 1 1)
  # Create FIXTURES_SETUP for Test_${executable}
  set_tests_properties(Test_${executable} PROPERTIES
    ENVIRONMENT "NANOS6_CONFIG=${CMAKE_BINARY_DIR}/nanos6.toml"
    FIXTURES_SETUP ST_${executable})

  # Add result test
  if (Python_FOUND AND Python_Interpreter_FOUND AND Python_NumPy_FOUND)
    add_test(NAME Test_${executable}_results
      COMMAND ${Python_EXECUTABLE} validate.py
      ${executable}_A.mat ${executable}_B.mat ${executable}_C.mat)
    set_tests_properties(Test_${executable}_results PROPERTIES
      FIXTURES_REQUIRED ST_${executable})
  endif () # Python_FOUND

  configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_${TEST}_dim.sh @ONLY)
  configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_${TEST}_dim.sh @ONLY)

endforeach ()

# Submit
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_mn.sh" submit_mn.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_mn.sh" submiter_mn.sh @ONLY)

# Extrae
configure_file("${PROJECT_SOURCE_DIR}/utils/trace.sh" trace.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/extrae_template.xml" extrae.xml @ONLY)

# Process
# configure_file("${PROJECT_SOURCE_DIR}/utils/process_mn.py" process_mn.py COPYONLY)
# configure_file("${PROJECT_SOURCE_DIR}/utils/process_dim.py" process_dim.py COPYONLY)
