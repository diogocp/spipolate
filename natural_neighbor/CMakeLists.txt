cmake_minimum_required(VERSION 3.1)
project(spipolate_natural_neighbor)

find_package(CGAL REQUIRED COMPONENTS Core)

if(APPLE)
    set(STATA_SYSTEM APPLEMAC)
elseif(UNIX)
    set(STATA_SYSTEM OPUNIX)
elseif(WIN32)
    set(STATA_SYSTEM STWIN32)
endif()

add_library(${PROJECT_NAME} MODULE stplugin.c natural_neighbor.cpp)
target_compile_definitions(${PROJECT_NAME} PRIVATE SYSTEM=${STATA_SYSTEM})
set_target_properties(${PROJECT_NAME} PROPERTIES
        PREFIX ""
        SUFFIX ".plugin"
        C_STANDARD 11
        C_EXTENSIONS OFF
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF)
