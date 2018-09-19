cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/utils/include_guard.cmake)
include_guard(utils)

# Utility modules.
add_library(utils utils.f90)

# Link to subsidiary libraries.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/precision precision/)
target_link_libraries(utils precision)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io_basic io_basic/)
target_link_libraries(utils io_basic)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/abstract abstract/)
target_link_libraries(utils abstract)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io io/)
target_link_libraries(utils io)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/algebra algebra/)
target_link_libraries(utils algebra)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/random random/)
target_link_libraries(utils random)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/arguments arguments/)
target_link_libraries(utils arguments)

target_include_directories(utils PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
