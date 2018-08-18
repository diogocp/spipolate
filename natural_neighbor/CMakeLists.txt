cmake_minimum_required(VERSION 3.4)
project(spipolate_natural_neighbor)

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


if(MSVC)
    # This seems to only be required on Windows, for some reason
    target_link_libraries(${PROJECT_NAME} PRIVATE CGAL::CGAL)

    # Link to the static CLR when using a *-windows-static vcpkg triplet
    if(VCPKG_TARGET_TRIPLET MATCHES "-windows-static$")
        foreach(config ${CMAKE_CONFIGURATION_TYPES} ${CMAKE_BUILD_TYPE})
            string(TOUPPER "${config}" config)
            foreach(lang C CXX)
                string(APPEND CMAKE_${lang}_FLAGS_${config} " /MT")
            endforeach(lang)
        endforeach(config)
    endif(VCPKG_TARGET_TRIPLET MATCHES "-windows-static$")
endif(MSVC)
