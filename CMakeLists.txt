cmake_minimum_required(VERSION 3.14)

project (WANN LANGUAGES C CXX)


set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# Making pz "seeable" in Xcode
# ----------------
file(GLOB_RECURSE PZ_SRC
     "../neopz/*.h"
     "../neopz/*.cpp"
     "../NeoPzDevelop/*.h"
     "../NeoPzDevelop/*.cpp"
)
add_custom_target( GlobalFiles SOURCES ${PZ_SRC})
# ----------------

add_subdirectory(3d-well)

# creating the projects
add_executable(wann wann.cpp)
target_link_libraries(wann PUBLIC NeoPZ::pz)
