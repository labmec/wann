# Standalone check executables, run manually (not registered with CTest).
set(WANN_CHECK_TARGETS
    hanging-nodes-test
    test-divide
    test-sideOrient-1d
    test-nonlinear-3d
    testDarcyNL
    test-nonlinear-1d
)

foreach(check_target IN LISTS WANN_CHECK_TARGETS)
    add_executable(${check_target} ${check_target}.cpp)
    target_link_libraries(${check_target} PRIVATE WANN_lib)
endforeach()
