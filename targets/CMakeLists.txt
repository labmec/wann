add_executable(wann3d wann3d.cpp)
target_link_libraries(wann3d PRIVATE WANN_lib)

add_executable(wann3dRef wann3dRef.cpp)
target_link_libraries(wann3dRef PRIVATE WANN_lib)

add_executable(oldwann3d old-wann3d.cpp)
target_link_libraries(oldwann3d PRIVATE WANN_lib)

add_executable(hanging-nodes-test hanging-nodes-test.cpp)
target_link_libraries(hanging-nodes-test PRIVATE WANN_lib)

add_executable(test-divide test-divide.cpp)
target_link_libraries(test-divide PRIVATE WANN_lib)

add_executable(test-sideOrient-1d test-sideOrient-1d.cpp)
target_link_libraries(test-sideOrient-1d PRIVATE WANN_lib)

add_executable(test-nonlinear test-nonlinear.cpp)
target_link_libraries(test-nonlinear PRIVATE WANN_lib)