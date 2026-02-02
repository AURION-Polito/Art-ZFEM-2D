add_subdirectory(Mesh)

list(APPEND elliptic_extended_PCC_2D_includes ${CMAKE_CURRENT_SOURCE_DIR})

set(elliptic_extended_PCC_2D_sources ${elliptic_extended_PCC_2D_sources} PARENT_SCOPE)
set(elliptic_extended_PCC_2D_headers ${elliptic_extended_PCC_2D_headers} PARENT_SCOPE)
set(elliptic_extended_PCC_2D_includes ${elliptic_extended_PCC_2D_includes} PARENT_SCOPE)
