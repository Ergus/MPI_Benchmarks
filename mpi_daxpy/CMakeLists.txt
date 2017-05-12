#To be save as -fopenmp is default option
STRING(REPLACE "-fopenmp " "" CMAKE_C_FLAGS ${CMAKE_C_FLAGS})
STRING(REPLACE "-fopenmp " "" CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})

# Selective compilation with openMP
ADD_LIBRARY(daxpy-mpi-omp OBJECT daxpy-mpi.c)
TARGET_COMPILE_OPTIONS(daxpy-mpi-omp PUBLIC "-fopenmp")
ADD_EXECUTABLE(daxpy-gather-omp driver.cxx $<TARGET_OBJECTS:daxpy-mpi-omp> )
TARGET_LINK_LIBRARIES(daxpy-gather-omp common )
SET_TARGET_PROPERTIES(daxpy-gather-omp PROPERTIES LINK_FLAGS "-fopenmp")

# COmpile secuential version
ADD_EXECUTABLE(daxpy-gather-sec driver.cxx daxpy-mpi.c)
TARGET_LINK_LIBRARIES(daxpy-gather-sec common)

# ADD tests
ADD_TEST(NAME Test_Daxpy_Gather-opm COMMAND mpirun -np 4 daxpy-gather-omp 16 8 3 omp)
ADD_TEST(NAME Test_Daxpy_Gather-sec COMMAND mpirun -np 4 daxpy-gather-sec 16 8 3 sec)

# python test
IF(PYTHONINTERP_FOUND)
  CONFIGURE_FILE(prove.py prove.py)  
  ADD_TEST(NAME Test_Daxpy_Gather_omp_results COMMAND 
    python3 prove.py 3 omp_X.mat omp_in_Y.mat omp_out_Y.mat)
  ADD_TEST(NAME Test_Daxpy_Gather_sec_results COMMAND 
    python3 prove.py 3 omp_X.mat omp_in_Y.mat omp_out_Y.mat)
ELSE()
  MESSAGE("Python interpreter not found")
ENDIF(PYTHONINTERP_FOUND)

# Submitter scripts and filters
CONFIGURE_FILE(submit.sh submit.sh COPYONLY)
CONFIGURE_FILE(profile-daxpy.sh profile-daxpy.sh COPYONLY)
CONFIGURE_FILE(${PROJECT_SOURCE_DIR}/mpi_matmul/gnuplot.plt gnuplot.plt COPYONLY)

