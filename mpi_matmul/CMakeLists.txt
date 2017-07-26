# Compilation

set(execut matmul-allgather)

ADD_EXECUTABLE(${execut} driver.cxx matmul-mpi.c)
TARGET_LINK_LIBRARIES(${execut} common)

# Run test
ADD_TEST(NAME Test_${execut} COMMAND mpirun -np 4 ./${execut} 16 8 matrix)  

# Submitter scripts
CONFIGURE_FILE(submit.sh submit.sh COPYONLY)
CONFIGURE_FILE(profile-matmul.sh profile-matmul.sh COPYONLY)
CONFIGURE_FILE(gnuplot.plt gnuplot.plt COPYONLY)

# python test
IF(PYTHONINTERP_FOUND)
  CONFIGURE_FILE(prove.py prove.py)  
  ADD_TEST(NAME Test_${execut}_results COMMAND 
    python3 prove.py matrix_A.mat matrix_B.mat matrix_C.mat)
ELSE()
  MESSAGE("Python interpreter not found")
ENDIF(PYTHONINTERP_FOUND)

