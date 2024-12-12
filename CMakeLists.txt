cmake_minimum_required(VERSION 3.29)
project(labelling)

set(CMAKE_CXX_STANDARD 20)
include_directories("lib/libBV/include")
add_executable(labelling main.cpp)
add_library(libBV STATIC lib/libBV/lib/Grafik.cxx
        lib/libBV/lib/SimpleFont.cxx
        lib/libBV/lib/AviRead.cxx
        lib/libBV/lib/AviWrite.cxx
        lib/libBV/lib/BmpRead.cxx
        lib/libBV/lib/BmpWrite.cxx
        lib/libBV/lib/FFT.cxx
        lib/libBV/lib/FileIO.cxx
        lib/libBV/lib/Grafik.cxx
        lib/libBV/lib/ConvImg.cxx
)
target_link_libraries(labelling libBV)