cmake_minimum_required(VERSION 2.8)
project( mgranul )
add_definitions("-std=c++11")
find_package( OpenCV REQUIRED )
include_directories(${OpenCV_INCLUDE_DIRS})

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -lpng -lz -ltiff -ljasper -ljpeg -lgtk-x11-2.0 -lgdk-x11-2.0 -latk-1.0 -lgio-2.0 -lpangoft2-1.0 -lpangocairo-1.0 -lgdk_pixbuf-2.0 -lcairo -lpango-1.0 -lfontconfig -lgobject-2.0 -lfreetype -lgthread-2.0 -lglib-2.0 -ldc1394 -lavcodec -lavformat -lavutil -lswscale -ldl -lm -lpthread -lrt -fopenmp")

add_library( lib_granul

             STATIC

             lib/granul.cpp
             lib/granul.h
             lib/meikon.h)

add_executable(mgranul main.cpp)

target_link_libraries( mgranul ${OpenCV_LIBS} lib_granul)

