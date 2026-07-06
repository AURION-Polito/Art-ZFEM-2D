list(APPEND Parabolic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND Parabolic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND Parabolic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)
list(APPEND Parabolic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)

list(APPEND Parabolic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND Parabolic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND Parabolic_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(Parabolic_PCC_2D_sources ${Parabolic_PCC_2D_sources} PARENT_SCOPE)
set(Parabolic_PCC_2D_headers ${Parabolic_PCC_2D_headers} PARENT_SCOPE)
set(Parabolic_PCC_2D_includes ${Parabolic_PCC_2D_includes} PARENT_SCOPE)
