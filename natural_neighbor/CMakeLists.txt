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
    # vcpkg renames libraries, breaking auto-linking when using a *-windows-static triplet
    target_compile_definitions(${PROJECT_NAME} PRIVATE CGAL_NO_AUTOLINK)

    # Link to the static CLR
    foreach(config ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE})
        string(TOUPPER "${config}" config)
        foreach(lang C CXX)
            string(APPEND CMAKE_${lang}_FLAGS_${config} " /MT")
        endforeach(lang)
    endforeach(config)
endif()
