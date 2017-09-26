# Compilation
SET(execut matvec-allgather)

ADD_EXECUTABLE(${execut} driver.cxx matvec-mpi.c)
TARGET_LINK_LIBRARIES(${execut} common)

# Run test
ADD_TEST(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 16 8 matrix)  

# Submitter scripts
#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/mpi_matmul/submit.sh submit.sh COPYONLY)
#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/mpi_matmul/profile-matmul.sh profile-matmul.sh COPYONLY)
#CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/mpi_matmul/gnuplot.plt gnuplot.plt COPYONLY)

CONFIGURE_FILE(run_matvec.sh run_matvec.sh)
CONFIGURE_FILE(profile_matvec.sh profile_matvec.sh @ONLY)

# python test
IF(PYTHONINTERP_FOUND)
  CONFIGURE_FILE(prove.py prove.py)  
  ADD_TEST(NAME Test_${execut}_results COMMAND 
    python3 prove.py matrix_A.mat matrix_B.mat matrix_C.mat)
ELSE()
  MESSAGE("Python interpreter not found")
ENDIF(PYTHONINTERP_FOUND)

