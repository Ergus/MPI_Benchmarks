# Compilation

add_definitions(-DUSE_MPI -DOMPI_SKIP_MPICXX)

file(GLOB MAIN_SOURCES mpi*.c)

find_package(MKL REQUIRED)

include_directories(${MKL_INCLUDES})
add_definitions("-DHAVE_MKL_H")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${MKL_EXTRA_LINKS}")

foreach (mainfile ${MAIN_SOURCES})

  get_filename_component(execut ${mainfile} NAME_WE)  # delete path

  add_executable(${execut} ${mainfile} util.c)
  target_link_libraries(${execut} common ${MKL_LIBRARIES})

  add_test(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 16 4 1)

endforeach (mainfile ${MAIN_SOURCES})

# Submitter scripts
# CONFIGURE_FILE(submit.sh submit.sh COPYONLY)
# CONFIGURE_FILE(profile-matmul.sh profile-matmul.sh COPYONLY)
# CONFIGURE_FILE(gnuplot.plt gnuplot.plt COPYONLY)

# python test
# IF(PYTHONINTERP_FOUND)
#   CONFIGURE_FILE(prove.py prove.py)  
#   ADD_TEST(NAME Test_${execut}_results COMMAND 
#     python3 prove.py matrix_A.mat matrix_B.mat matrix_C.mat)
# ELSE()
#   MESSAGE("Python interpreter not found")
# ENDIF(PYTHONINTERP_FOUND)
