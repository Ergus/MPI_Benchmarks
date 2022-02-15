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

set(TASKLIST serial parallelfor task)
foreach(TASKTYPE 0 1 2)
  if ((NOT WITH_SERIAL) AND (TASKTYPE EQUAL 0)) # No serial unless WITH_SERIAL
    continue()
  endif()

  list(GET TASKLIST ${TASKTYPE} TASKSUFIX)

    set(P2PLIST nop2p p2p)
    foreach(P2P 0 1)
      list(GET P2PLIST ${P2P} P2PSUFIX)

      if ((TASKTYPE EQUAL 0) AND (P2P GREATER 0))
	continue()
      endif()

      set(BLASLIST noblas blas)
      foreach(BLAS 0 1)
	list(GET BLASLIST ${BLAS} BLASSUFIX)

	set(EXECUTABLE "jacobi_${TASKSUFIX}_${P2PSUFIX}_${BLASSUFIX}_mpi")

	message("- Add benchmark: ${EXECUTABLE}")

	add_executable(${EXECUTABLE} jacobi_mpi.c)
	target_compile_definitions(${EXECUTABLE}
	  PUBLIC P2P=${P2P}
	  PUBLIC TASKTYPE=${TASKTYPE}
	  PUBLIC BLAS=${BLAS})

	add_test(NAME "Test_${EXECUTABLE}"
	  COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ./${EXECUTABLE} 128 8 3)
    endforeach()
  endforeach()
endforeach()

# submiter
configure_file("${PROJECT_SOURCE_DIR}/utils/submit_dim.sh" submit_dim.sh @ONLY)
configure_file("${PROJECT_SOURCE_DIR}/utils/submiter_dim.sh" submiter_dim.sh @ONLY)

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

