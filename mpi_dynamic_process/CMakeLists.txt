# Compilation
SET(execut dynamicOOP)

FILE(GLOB_RECURSE SOURCES "*.cxx")

ADD_EXECUTABLE(${execut} main.cxx ${SOURCES})
TARGET_LINK_LIBRARIES(${execut} readline common)

# Run test
ADD_TEST(NAME Test_${execut} COMMAND ${execut} 5 2)  

# Submitter scripts and filters
CONFIGURE_FILE(test1.sh test1.sh COPYONLY)
