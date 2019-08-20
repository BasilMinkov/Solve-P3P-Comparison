cmake_minimum_required(VERSION 2.8)
project( compareSloveP3P )
find_package( OpenCV REQUIRED )
find_package (Eigen3 3.3 REQUIRED NO_MODULE)  #
include_directories( ${OpenCV_INCLUDE_DIRS} )
add_executable( compareSloveP3P compareSloveP3P.cpp )
target_link_libraries( compareSloveP3P ${OpenCV_LIBS} Eigen3::Eigen)