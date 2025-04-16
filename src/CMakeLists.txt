# let us create a config file with the directory containing the input files
set(INPUTSDIR ${PROJECT_SOURCE_DIR}/input)
set(CONFIG_H "${CMAKE_CURRENT_SOURCE_DIR}/input_config.h")
configure_file ("${CMAKE_CURRENT_SOURCE_DIR}/input_config.h.in" ${CONFIG_H})

target_include_directories(WANN_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

file(GLOB WANN_SRC CONFIGURE_DEPENDS "*.cpp" "*.c" "*.hpp")
file(GLOB WANN_HEADERS CONFIGURE_DEPENDS "*.h")

target_sources(WANN_lib PRIVATE ${WANN_SRC} ${WANN_HEADERS})