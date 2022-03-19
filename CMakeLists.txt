cmake_minimum_required(VERSION 3.5)
project(gmm)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(ament_cmake REQUIRED)

add_library(${PROJECT_NAME}
  GMM.h
  GMM.cpp
  Kmeans.h
  Kmeans.cpp
  Matrix.h
  Matrix.cpp
)

add_executable(${PROJECT_NAME}_test
  main.cpp
  GMM.h
  GMM.cpp
  Kmeans.h
  Kmeans.cpp
  Matrix.h
  Matrix.cpp
)

install(
  TARGETS ${PROJECT_NAME}
  DESTINATION lib/${PROJECT_NAME}
)

install(
  TARGETS ${PROJECT_NAME}_test
  DESTINATION bin/${PROJECT_NAME}
)

ament_package()
