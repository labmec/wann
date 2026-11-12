add_executable(wann3dNL wann3dNL.cpp)
target_link_libraries(wann3dNL PRIVATE WANN_lib)

add_executable(radialFlowTest radialFlow.cpp)
target_link_libraries(radialFlowTest PRIVATE WANN_lib)

add_executable(adaptiveRefinement adaptiveRefinement.cpp)
target_link_libraries(adaptiveRefinement PRIVATE WANN_lib)
