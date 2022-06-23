cmake_minimum_required(VERSION 3.1...3.15)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project( Facade )

# set (CMAKE_CXX_STANDARD 11)
set_property(GLOBAL PROPERTY USE_FOLDERS ON)

### Configuration
SET(CMAKE_CXX_FLAGS_DEBUG "$ENV{CXXFLAGS} -O0 -Wall -g2 -ggdb") 
SET(CMAKE_CXX_FLAGS_RELEASE "$ENV{CXXFLAGS} -O3 -Wall")
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
set(FACADE_ROOT         ${CMAKE_CURRENT_LIST_DIR})
set(FACADE_INCLUDE_DIR  ${FACADE_ROOT}  ${CMAKE_CURRENT_BINARY_DIR})
set(FACADE_SOURCE_DIR   ${FACADE_ROOT})


set(FACADE_scip_DIR                ${FACADE_ROOT}/3rd_scip)
set(FACADE_soplex_DIR              ${FACADE_ROOT}/3rd_soplex)

# add external target
add_library(3rd_party INTERFACE)
target_include_directories(3rd_party INTERFACE 3rd_party/)

add_subdirectory(3rd_scip)
add_subdirectory(3rd_soplex)

add_subdirectory(basic)
add_subdirectory(math)
add_subdirectory(Facade)

 

