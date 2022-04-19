# UseCompass.cmake
#
# Builds the GCOptimization library
#
# This will define the target
#
#     libGCO

add_library(libGCO
  ${CMAKE_SOURCE_DIR}/3rd_party/gco-v3.0/GCoptimization.cpp
  ${CMAKE_SOURCE_DIR}/3rd_party/gco-v3.0/graph.cpp
  ${CMAKE_SOURCE_DIR}/3rd_party/gco-v3.0/LinkedBlockList.cpp
  ${CMAKE_SOURCE_DIR}/3rd_party/gco-v3.0/maxflow.cpp
)
target_include_directories(libGCO PUBLIC ${CMAKE_SOURCE_DIR}/3rd_party/)
message(STATUS "Build libGCO library")

