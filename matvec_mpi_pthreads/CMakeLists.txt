# Compilation
SET(execut matvec-mpi-pthreads)

ADD_EXECUTABLE(${execut} driver.cxx ${execut}.c ${execut}.h)
TARGET_LINK_LIBRARIES(${execut} common pthread)

# Run test
ADD_TEST(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 64 8 1 1 matrix)

# python test
IF(PYTHONINTERP_FOUND)
  CONFIGURE_FILE(prove.py prove.py)
  ADD_TEST(NAME Test_${execut}_results COMMAND
    python3 prove.py matrix_A.mat matrix_B.mat matrix_C.mat)
ELSE()
  MESSAGE("Python interpreter not found")
ENDIF(PYTHONINTERP_FOUND)

# Submitter scripts
#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/matmul_mpi/submit.sh submit.sh COPYONLY)
#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/matmul_mpi/profile-matmul.sh profile-matmul.sh COPYONLY)
#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/matmul_mpi/gnuplot.plt gnuplot.plt COPYONLY)

#CONFIGURE_FILE(run_matvec.sh run_matvec.sh @ONLY)
#CONFIGURE_FILE(profile_matvec.sh profile_matvec.sh @ONLY)
