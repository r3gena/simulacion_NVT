cmake_minimum_required(VERSION 3.13)

project(NVT)

include_directories(${PROJECT_SOURCE_DIR})
add_library(applib progress_bar.cpp potenciales.cpp Sistema_NVT.cpp Sistema.cpp)
add_executable(NVT simulacion_NVT.cpp)
target_link_libraries(NVT applib)
