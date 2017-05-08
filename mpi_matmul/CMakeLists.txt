# Compilation
ADD_EXECUTABLE(matmul_mpi driver.cxx matmul-mpi.c)
TARGET_LINK_LIBRARIES(matmul_mpi common)

# Run test
ADD_TEST(NAME Test_Matmul_run COMMAND mpirun -np 4 matmul_mpi 16 8 true)  

# python test
IF(PYTHONINTERP_FOUND)
  CONFIGURE_FILE(prove.py prove.py)  
  ADD_TEST(NAME Test_Matmul_result COMMAND 
    prove.py matrix_A.mat matrix_B.mat matrix_C.mat)
ELSE()
  MESSAGE("Python interpreter not found")
ENDIF(PYTHONINTERP_FOUND)
