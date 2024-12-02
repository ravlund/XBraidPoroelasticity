cmake_minimum_required(VERSION 3.5)

set(pluginName XBraidPoroelasticity)
project(UG_PLUGIN_${pluginName})

set(UG4_XBraidPoroelasticity_plugin_VERSION_MAJOR 0.1)
set(UG4_XBraidPoroelasticity_plugin_VERSION_MINOR 0.0)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
#set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_INIT} -fPIC")
#set(CMAKE_CXX_FLAGS " ${CMAKE_CXX_FLAGS_INIT} -Wpedantic")
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(SOURCES
        src/xbraidporoelasticity_plugin.cpp
        src/braid_biot_estimator.h
        src/braid_biot_control.h
        src/braid_biot_precomputed.h
        src/braid_heat_check.h
        src/biot_error_data.h
        src/biot_braid_displacement_norm.h
)


if(NOT UG_ROOT_CMAKE_PATH)
    set(UG_ROOT_CMAKE_PATH /storage/pro/ug/ugcore/cmake)
endif(NOT UG_ROOT_CMAKE_PATH)

include(${UG_ROOT_CMAKE_PATH}/ug_plugin_includes.cmake)

if (buildEmbeddedPlugins)
    EXPORTSOURCES(${CMAKE_CURRENT_SOURCE_DIR} ${SOURCES})
else (buildEmbeddedPlugins)
    add_library(${pluginName} SHARED ${SOURCES})

    link_directories(${UG_ROOT_PATH}/bin/plugins)
    target_link_libraries(${pluginName} PUBLIC ug4)
    target_link_libraries(${pluginName} PUBLIC Poroelasticity)
    target_link_libraries(${pluginName} PUBLIC XBraidForUG4)
endif (buildEmbeddedPlugins)





