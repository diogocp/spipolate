cmake_minimum_required(VERSION 3.4)
project(spipolate_natural_neighbor)

option(USE_STATIC_LIBS "Build using static libraries" ON)
if(USE_STATIC_LIBS)
    set(CGAL_Boost_USE_STATIC_LIBS ON)
endif()

find_package(CGAL CONFIG REQUIRED)

if(APPLE)
    set(STATA_SYSTEM APPLEMAC)
elseif(UNIX)
    set(STATA_SYSTEM OPUNIX)
elseif(WIN32)
    set(STATA_SYSTEM STWIN32)
endif()

add_library(${PROJECT_NAME} MODULE stplugin.c natural_neighbor.cpp)
target_compile_definitions(${PROJECT_NAME} PRIVATE SYSTEM=${STATA_SYSTEM})
target_compile_definitions(${PROJECT_NAME} PRIVATE CGAL_NO_AUTOLINK)
target_include_directories(${PROJECT_NAME} PRIVATE ${CGAL_INCLUDE_DIRS})
set_target_properties(${PROJECT_NAME} PROPERTIES
        PREFIX ""
        SUFFIX ".plugin"
        C_STANDARD 11
        C_EXTENSIONS OFF
        CXX_STANDARD 11
        CXX_STANDARD_REQUIRED ON
        CXX_EXTENSIONS OFF)

if(NOT CMAKE_COMPILER_IS_GNUCXX)
    # For some reason, GCC does not require this
    target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL)
endif()

if(MSVC AND USE_STATIC_LIBS)
    # Link to the static CLR
    target_compile_options(${PROJECT_NAME} PRIVATE "/MT")
endif()
