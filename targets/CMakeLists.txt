add_executable(wann3d wann3d.cpp)
target_link_libraries(wann3d PRIVATE WANN_lib)

add_executable(oldwann3d old-wann3d.cpp)
target_link_libraries(oldwann3d PRIVATE WANN_lib)

add_executable(hanging-nodes-test hanging-nodes-test.cpp)
target_link_libraries(hanging-nodes-test PRIVATE WANN_lib)