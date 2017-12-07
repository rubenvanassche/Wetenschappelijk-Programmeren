cmake_minimum_required(VERSION 3.8)
project(Orthogonal_Basisfunctions)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/Extra")

find_package(GSL)

INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIRS})

set(SOURCE_FILES main.cpp PointsWriter.cpp PointsWriter.h)
add_executable(main ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(main ${GSL_LIBRARIES})