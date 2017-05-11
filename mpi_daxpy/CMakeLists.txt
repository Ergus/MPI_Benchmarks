# Compilation
ADD_EXECUTABLE(daxpy-gather driver.cxx daxpy-mpi.c)
TARGET_LINK_LIBRARIES(daxpy-gather common)

# Run test
ADD_TEST(NAME Test_Daxpy_Gather COMMAND mpirun -np 4 daxpy-gather 16 8 3 vector)  

# Submitter scripts
CONFIGURE_FILE(submit.sh submit.sh COPYONLY)
CONFIGURE_FILE(profile-daxpy.sh profile-daxpy.sh COPYONLY)
CONFIGURE_FILE(gnuplot.plt gnuplot.plt COPYONLY)

# python test
IF(PYTHONINTERP_FOUND)
  CONFIGURE_FILE(prove.py prove.py)  
  ADD_TEST(NAME Test_Daxpy_Gather_results COMMAND 
    python3 prove.py 3 vector_X.mat vector_in_Y.mat vector_out_Y.mat)
ELSE()
  MESSAGE("Python interpreter not found")
ENDIF(PYTHONINTERP_FOUND)

