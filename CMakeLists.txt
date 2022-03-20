cmake_minimum_required(VERSION 3.5)
project(gmm)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(ament_cmake REQUIRED)

add_library(${PROJECT_NAME}
  src/GMM.cpp
  src/Kmeans.cpp
  src/Matrix.cpp
)
target_include_directories(${PROJECT_NAME}
  PRIVATE
  include
)

add_executable(${PROJECT_NAME}_test
  src/main.cpp
  src/GMM.cpp
  src/Kmeans.cpp
  src/Matrix.cpp
)
target_include_directories(${PROJECT_NAME}_test
  PRIVATE
  include
)

install(
  TARGETS ${PROJECT_NAME}
  DESTINATION lib/${PROJECT_NAME}
)
install(DIRECTORY include/
  DESTINATION include
)
install(
  TARGETS ${PROJECT_NAME}_test
  DESTINATION bin/${PROJECT_NAME}
)

ament_package()
