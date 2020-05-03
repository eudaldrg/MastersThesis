cmake_minimum_required(VERSION 3.15)

set(PROJECT_NAME SWIFT)
project(${PROJECT_NAME})

set(CMAKE_CXX_STANDARD 17)

set(BOOST_ROOT C:\\Program Files\\boost\\boost_1_71_0)
message(${BOOST_ROOT})
find_package(Boost 1.71.0 REQUIRED #[[COMPONENTS system]])

#add_library(${PROJECT_NAME} library.cpp library.h)
target_include_directories(${PROJECT_NAME} PRIVATE ${Boost_INCLUDE_DIRS})
add_executable(${PROJECT_NAME} library.cpp library.h)
