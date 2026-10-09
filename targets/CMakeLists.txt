add_executable(wann3dNL wann3dNL.cpp)
target_link_libraries(wann3dNL PRIVATE WANN_lib)

add_executable(radial radialFlow.cpp)
target_link_libraries(radial PRIVATE WANN_lib)

add_executable(visualize visualizeMesh.cpp)
target_link_libraries(visualize PRIVATE WANN_lib)

add_executable(twoWells twoWells.cpp)
target_link_libraries(twoWells PRIVATE WANN_lib)
