cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(basic_HEADERS
    assertions.h
    basic_common.h
    basic_types.h
    logger.h
    Config.h
    )

set(basic_SOURCES
    assertions.cpp
    basic_types.cpp
    )

set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
find_package(Glog REQUIRED)
include_directories(${GLOG_INCLUDE_DIRS})

add_library(${PROJECT_NAME} SHARED ${basic_SOURCES} ${basic_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "Facade")


target_compile_definitions(${PROJECT_NAME} PRIVATE BASIC_EXPORTS)


