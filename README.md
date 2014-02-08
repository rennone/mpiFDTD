mpiFDTD
=======

FDTD simulator by MPI Library

以下 CMakeLists.txt のサンプル

#要求するcmakeのバージョン
cmake_minimum_required(VERSION 2.8)

#project name
PROJECT(mpi_fdtd2D)

set(CMAKE_C_FLAGS "-Wall -std=c99")

# find MPI
find_package(MPI)
find_package(OPENGL)
find_package(GLUT)
find_package(GLEW)

include_directories(${MPI_C_INCLUDE_PATH})
include_directories(~/glew-1.10.0/include)  #include_directories(${GLEW_INCLUDE_PATH})
include_directories(/usr/local/include/)    #include_directories(${GLUT_INCLUDE_DIR})
include_directories(/opt/local/include)

link_directories(${MPI_LIBRARY_PATH})
link_directories(~/glew-1.10.0/lib)     #link_directories(${GLEW_LIBRARY_PATH})
#link_directories(${GLUT_LIBRARY_PATH})

set(CMAKE_BUILD_TYPE release)
add_executable(mpi_fdtd2D main.c field.c simulator.c mpiTM_UPML.c mpiTE_UPML.c models.c noModel.c circleModel.c drawer.c)
target_link_libraries(mpi_fdtd2D ${MPI_LIBRARIES})
#target_link_libraries(mpi_fdtd2D ${GLUT_LIBRARY} ${OPENGL_LIBRARY} m)
target_link_libraries(mpi_fdtd2D /System/Library/Frameworks/GLUT.framework ${OPENGL_LIBRARY} m)
target_link_libraries(mpi_fdtd2D GLEW.a)
