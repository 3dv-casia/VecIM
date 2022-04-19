cmake_minimum_required(VERSION 3.1)

get_filename_component(PROJECT_NAME ${CMAKE_CURRENT_SOURCE_DIR} NAME)
project(${PROJECT_NAME})


set(math_HEADERS    
    math_common.h   
    linear_program.h
    linear_program_solver.h
    )

set(math_SOURCES    
    linear_program.cpp
    linear_program_solver.cpp
    linear_program_solver_SCIP.cpp
    )



add_library(${PROJECT_NAME} SHARED ${math_SOURCES} ${math_HEADERS})
set_target_properties(${PROJECT_NAME} PROPERTIES
                      FOLDER "Facade")


target_include_directories(${PROJECT_NAME} PRIVATE
  ${FACADE_scip_DIR}
  ${FACADE_scip_DIR}/scip
  )

target_compile_definitions(${PROJECT_NAME} PRIVATE MATH_EXPORTS)



target_link_libraries( ${PROJECT_NAME} PRIVATE basic 3rd_scip  3rd_soplex boost_system)
