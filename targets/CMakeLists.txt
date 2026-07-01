add_executable(wann3d wann3d.cpp)
target_link_libraries(wann3d PRIVATE WANN_lib)

add_executable(wann3dRef wann3dRef.cpp)
target_link_libraries(wann3dRef PRIVATE WANN_lib)

add_executable(wann3dNL wann3dNL.cpp)
target_link_libraries(wann3dNL PRIVATE WANN_lib)

add_executable(radial radialFlowTest.cpp)
target_link_libraries(radial PRIVATE WANN_lib)

add_executable(wannPS wannParameterStudy.cpp)
target_link_libraries(wannPS PRIVATE WANN_lib)

add_executable(adptRef adaptiveRefinement.cpp)
target_link_libraries(adptRef PRIVATE WANN_lib)

add_executable(goalDarcy goalOrientedDarcy.cpp)
target_link_libraries(goalDarcy PRIVATE WANN_lib)

add_executable(goalWann goalOrientedWann.cpp)
target_link_libraries(goalWann PRIVATE WANN_lib)

add_executable(oldwann3d old-wann3d.cpp)
target_link_libraries(oldwann3d PRIVATE WANN_lib)

add_executable(hanging-nodes-test hanging-nodes-test.cpp)
target_link_libraries(hanging-nodes-test PRIVATE WANN_lib)

add_executable(test-divide test-divide.cpp)
target_link_libraries(test-divide PRIVATE WANN_lib)

add_executable(test-sideOrient-1d test-sideOrient-1d.cpp)
target_link_libraries(test-sideOrient-1d PRIVATE WANN_lib)

add_executable(test-nonlinear-3d test-nonlinear-3d.cpp)
target_link_libraries(test-nonlinear-3d PRIVATE WANN_lib)

add_executable(testDarcyNL testDarcyNL.cpp)
target_link_libraries(testDarcyNL PRIVATE WANN_lib)

add_executable(test-nonlinear-1d test-nonlinear-1d.cpp)
target_link_libraries(test-nonlinear-1d PRIVATE WANN_lib)