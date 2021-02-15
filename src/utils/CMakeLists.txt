cmake_minimum_required (VERSION 3.5 FATAL_ERROR)
include(${CMAKE_SOURCE_DIR}/cmake/cmake_utils.cmake)
cmake_include_guard(utils)

# Utility modules.
add_library(utils_modules utils.f90)
add_library(utils ${DUMMY_FILE})
target_link_libraries(utils utils_modules)

# Link to subsidiary libraries.
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/foundations foundations/)
target_link_libraries(utils_modules foundations_modules)
target_link_libraries(utils foundations_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/macros macros/)
target_link_libraries(utils_modules macros_modules)
target_link_libraries(utils macros_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/abstract abstract/)
target_link_libraries(utils_modules abstract_modules)
target_link_libraries(utils abstract_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/io io/)
target_link_libraries(utils_modules io_modules)
target_link_libraries(utils io_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/algebra algebra/)
target_link_libraries(utils_modules algebra_modules)
target_link_libraries(utils algebra_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/random random/)
target_link_libraries(utils_modules random_modules)
target_link_libraries(utils random_submodules)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/arguments arguments/)
target_link_libraries(utils_modules arguments_modules)
target_link_libraries(utils arguments_submodules)

target_include_directories(utils_modules PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
