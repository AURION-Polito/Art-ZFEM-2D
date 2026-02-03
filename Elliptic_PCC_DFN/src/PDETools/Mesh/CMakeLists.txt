list(APPEND Elliptic_PCC_DFN_headers ${CMAKE_CURRENT_SOURCE_DIR}/Extended_PDE_Mesh_Utilities.hpp)

list(APPEND Elliptic_PCC_DFN_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(Elliptic_PCC_DFN_sources ${Elliptic_PCC_DFN_sources} PARENT_SCOPE)
set(Elliptic_PCC_DFN_headers ${Elliptic_PCC_DFN_headers} PARENT_SCOPE)
set(Elliptic_PCC_DFN_includes ${Elliptic_PCC_DFN_includes} PARENT_SCOPE)
