CMAKE_MINIMUM_REQUIRED ( VERSION 3.1.0 )

project(sfc CXX)

set (CMAKE_CXX_STANDARD 17)

include_directories(include)
include_directories(ext/zorder_knn/zorder_knn/include)

add_executable(test1 src/test1.cpp)
add_executable(test2 src/test2.cpp)
