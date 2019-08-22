cmake_minimum_required(VERSION 2.8) # Check CMake version. If the version of the installed program older than specified, an emergency exit will occur.

project( compareSloveP3P ) # Project name

find_package( OpenCV REQUIRED ) 
find_package (Eigen3 3.3 REQUIRED NO_MODULE)

include_directories(include)
include_directories( ${OpenCV_INCLUDE_DIRS}) # Location of header files

add_executable(
	compareSloveP3P 
	scr/compareSloveP3P.cpp 
	scr/TdS_Math/TdS_MathBasic.cpp 
	scr/TdS_Math/poly34_tsC.cpp 
	scr/P3P/P3P_MinkovSawada.cpp
	)

target_link_libraries( compareSloveP3P ${OpenCV_LIBS} Eigen3::Eigen )