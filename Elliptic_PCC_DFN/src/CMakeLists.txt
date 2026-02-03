add_subdirectory(PDETools)

list(APPEND Elliptic_PCC_DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_configuration.hpp)
list(APPEND Elliptic_PCC_DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.hpp)
list(APPEND Elliptic_PCC_DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/assembler.hpp)
list(APPEND Elliptic_PCC_DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/test_definition.hpp)

list(APPEND Elliptic_PCC_DFN_sources ${CMAKE_CURRENT_SOURCE_DIR}/assembler.cpp)
list(APPEND Elliptic_PCC_DFN_sources ${CMAKE_CURRENT_SOURCE_DIR}/program_utilities.cpp)

list(APPEND Elliptic_PCC_DFN_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(Elliptic_PCC_DFN_sources ${Elliptic_PCC_DFN_sources} PARENT_SCOPE)
set(Elliptic_PCC_DFN_headers ${Elliptic_PCC_DFN_headers} PARENT_SCOPE)
set(Elliptic_PCC_DFN_includes ${Elliptic_PCC_DFN_includes} PARENT_SCOPE)
