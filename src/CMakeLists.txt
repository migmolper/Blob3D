set(BLob3D_src
  ${CMAKE_CURRENT_SOURCE_DIR}/Simulation.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Variables.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Blobs/Blobs.cpp  
  ${CMAKE_CURRENT_SOURCE_DIR}/Blobs/Blob-Function.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Blobs/Ghosts.cpp  
  ${CMAKE_CURRENT_SOURCE_DIR}/Blobs/Neighbors.cpp  
  ${CMAKE_CURRENT_SOURCE_DIR}/Blobs/IO/dump-inputs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Blobs/IO/vtk-outputs.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Boundaries/sphere.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Potentials/Advection-Diff-OpenMP.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Solvers/Mass-Transport-JKO-TAO.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Mesh/Boundary-Conditions.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Mesh/Coordinates.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Mesh/Create-Mesh.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Mesh/InOut-Mesh.cpp
)

add_library(blob3D STATIC ${BLob3D_src})
add_library(BLOB3D::blob3D ALIAS blob3D)

target_include_directories(blob3D
  PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<BUILD_INTERFACE:${EIGEN_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

target_link_libraries(blob3D PUBLIC ${BLOB3D_LINK_LIBS})
target_compile_definitions(blob3D PUBLIC GH3TH_INTEGRAL)
target_compile_definitions(blob3D PUBLIC NUMERICAL_DERIVATIVES=0)
if(BLOB3D_SPLINES_DIR_DEFAULT)
  target_compile_definitions(
    blob3D PUBLIC BLOB3D_SPLINES_DIR_DEFAULT="${BLOB3D_SPLINES_DIR_DEFAULT}")
endif()

include(GNUInstallDirs)

install(TARGETS blob3D
  EXPORT BLOB3DTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
  DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/BLOB3D
  FILES_MATCHING PATTERN "*.hpp"
)