# Compilation
SET(execut dynamicOOP)

FILE(GLOB_RECURSE SOURCES "*.cxx")

ADD_EXECUTABLE(${execut} main.cxx ${SOURCES})
TARGET_LINK_LIBRARIES(${execut} readline common)

# Run test
ADD_TEST(NAME Test_${execut} COMMAND ${execut} 5 2)  

# Submitter scripts and filters
CONFIGURE_FILE(test_spawn.sh test_spawn.sh COPYONLY)
CONFIGURE_FILE(test_shrink.sh test_shrink.sh COPYONLY)

CONFIGURE_FILE(gnuplot.plt gnuplot.plt COPYONLY)
