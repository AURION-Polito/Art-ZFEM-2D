cmake_minimum_required(VERSION 3.12)
project(Art_ZFEM_2D LANGUAGES CXX C VERSION 1.0)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/gedim/GeDiM)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/PolyDiM)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Elliptic_PCC_2D)
add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/Elliptic_PCC_DFN)

include(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/gedim/cmake/cppcheck.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/gedim/cmake/clang-format.cmake)
include(${CMAKE_CURRENT_SOURCE_DIR}/PolyDiM/gedim/cmake/prepend_license.cmake)
