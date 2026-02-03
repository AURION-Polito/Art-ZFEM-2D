list(APPEND Elliptic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND Elliptic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND Elliptic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)
list(APPEND Elliptic_PCC_2D_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)

list(APPEND Elliptic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND Elliptic_PCC_2D_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND Elliptic_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(Elliptic_PCC_2D_sources ${Elliptic_PCC_2D_sources} PARENT_SCOPE)
set(Elliptic_PCC_2D_headers ${Elliptic_PCC_2D_headers} PARENT_SCOPE)
set(Elliptic_PCC_2D_includes ${Elliptic_PCC_2D_includes} PARENT_SCOPE)
