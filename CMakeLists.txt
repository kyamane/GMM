cmake_minimum_required(VERSION 3.5)
project(gmm)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_BUILD_TYPE Release)

find_package(ament_cmake REQUIRED)

add_library(${PROJECT_NAME} SHARED
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

ament_export_include_directories(
  include
)
ament_export_libraries(
  ${PROJECT_NAME}
)

install(TARGETS ${PROJECT_NAME} ${PROJECT_NAME}_test
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib
)
install(DIRECTORY include/
  DESTINATION include
)

ament_package()
