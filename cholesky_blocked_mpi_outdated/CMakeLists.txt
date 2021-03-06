# Compilation

add_definitions(-DUSE_MPI -DOMPI_SKIP_MPICXX)

# Search modules
find_package(MKL REQUIRED)
include_directories(${MKL_INCLUDES})
add_definitions("-DHAVE_MKL_H")

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MKL_EXTRA_LINKS}")

# Add all the files with the name pattern mpi*.c
file(GLOB MAIN_SOURCES mpi*.c)

foreach (mainfile ${MAIN_SOURCES})
  get_filename_component(execut ${mainfile} NAME_WE)  # delete path

  add_executable(${execut} ${mainfile} util.c)
  target_link_libraries(${execut} ${MKL_LIBRARIES})

  add_test(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 16 4 1)
  add_test(NAME "Test_${execut}_bigger" COMMAND mpirun -np 4 ./${execut} 128 64 1)
endforeach (mainfile ${MAIN_SOURCES})

# MKL benchmarks (with scalapack)
set(MKL_SCALAPACK ON)
find_package(MKL REQUIRED)

file(GLOB MAIN_SOURCES mkl*.c)

# Sequential
foreach (mainfile ${MAIN_SOURCES})
  get_filename_component(execut ${mainfile} NAME_WE)  # delete path
  set(execut "${execut}_sec")

  add_executable(${execut} ${mainfile} util.c)
  target_link_libraries(${execut} ${MKL_LIBRARIES})

  add_test(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 16 4 1)
  add_test(NAME "Test_${execut}_bigger" COMMAND mpirun -np 4 ./${execut} 128 64 1)
endforeach (mainfile ${MAIN_SOURCES})

# Multithreaded
set(MKL_MULTI_THREADED ON)
find_package(MKL REQUIRED)

set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MKL_EXTRA_LINKS}")

foreach (mainfile ${MAIN_SOURCES})
  get_filename_component(execut ${mainfile} NAME_WE)  # delete path
  set(execut "${execut}_par")

  add_executable(${execut} ${mainfile} util.c)
  target_link_libraries(${execut} ${MKL_LIBRARIES})

  add_test(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 16 4 1)
  add_test(NAME "Test_${execut}_bigger" COMMAND mpirun -np 4 ./${execut} 128 64 1)
endforeach (mainfile ${MAIN_SOURCES})

# Submitter scripts
file(GLOB CONF_SOURCES submit*.sh)
foreach (conffile ${CONF_SOURCES})
  get_filename_component(conffile ${conffile} NAME)  # delete path
  CONFIGURE_FILE(${conffile} ${conffile} COPYONLY)
endforeach()

configure_file(cholesky.awk cholesky.awk COPYONLY)
configure_file(cholesky_babis.awk cholesky_babis.awk COPYONLY)
