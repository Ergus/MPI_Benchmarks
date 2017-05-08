
FILE(GLOB libsources ${PROJECT_SOURCE_DIR}/common/*.cpp)

ADD_LIBRARY(common STATIC ${libsources})
