
# SET(execut k-mean-mpi)
# 
# ADD_EXECUTABLE(${execut} driver.cxx k-mean-mpi.c)
# TARGET_LINK_LIBRARIES(${execut} common)
# 
# # ADD tests
# ADD_TEST(NAME Test_${execut} COMMAND mpirun -np 4 ${execut} 16 8 3)
# 
# # python test
# IF(PYTHONINTERP_FOUND)
#   CONFIGURE_FILE(prove.py prove.py)  
#   ADD_TEST(NAME Test_Daxpy_Gather_omp_results COMMAND 
#     python3 prove.py 3 omp_X.mat omp_in_Y.mat omp_out_Y.mat)
#   ADD_TEST(NAME Test_Daxpy_Gather_sec_results COMMAND 
#     python3 prove.py 3 omp_X.mat omp_in_Y.mat omp_out_Y.mat)
# ELSE()
#   MESSAGE("Python interpreter not found")
# ENDIF(PYTHONINTERP_FOUND)
# 
# # Submitter scripts and filters
# CONFIGURE_FILE(submit.sh submit.sh COPYONLY)
# CONFIGURE_FILE(profile-daxpy.sh profile-daxpy.sh COPYONLY)
# CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/mpi_matmul/gnuplot.plt gnuplot.plt COPYONLY)
