cmake_minimum_required(VERSION 3.5)
project(long_predict)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES src/main.c src/suffixArray.c src/configFileParser.c src/configFileParser.h src/svm.h src/svm.cpp src/kseq.h)
add_executable(long_predict ${SOURCE_FILES})
target_link_libraries(long_predict z)