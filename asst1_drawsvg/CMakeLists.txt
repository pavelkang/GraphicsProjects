cmake_minimum_required(VERSION 3.3)
project(asst1_drawsvg)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")

set(SOURCE_FILES
    build/CMakeFiles/3.2.1/CompilerIdC/CMakeCCompilerId.c
    build/CMakeFiles/3.2.1/CompilerIdCXX/CMakeCXXCompilerId.cpp
    build/CMakeFiles/feature_tests.c
    build/CMakeFiles/feature_tests.cxx
    include/GLEW/glew.h
    include/GLEW/glxew.h
    include/GLEW/wglew.h
    include/GLFW/glfw.h
    include/GLFW/glfwnative.h
    include/base64.h
    include/CMU462.h
    include/color.h
    include/complex.h
    include/matrix3x3.h
    include/osdtext.h
    include/renderer.h
    include/tinyxml2.h
    include/vector2D.h
    include/vector3D.h
    include/viewer.h
    src/drawsvg.cpp
    src/drawsvg.h
    src/hardware_renderer.cpp
    src/hardware_renderer.h
    src/illustrator.cpp
    src/illustrator.h
    src/main.cpp
    src/png.cpp
    src/png.h
    src/software_renderer.cpp
    src/software_renderer.h
    src/svg.cpp
    src/svg.h
    src/svg_renderer.h
    src/texture.cpp
    src/texture.h
    src/triangulation.cpp
    src/triangulation.h
    src/viewport.cpp
    src/viewport.h)

add_executable(asst1_drawsvg ${SOURCE_FILES})