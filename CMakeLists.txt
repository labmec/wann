cmake_minimum_required(VERSION 3.14)

project (WANN LANGUAGES C CXX)

set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Finding the pz package
find_package(NeoPZ REQUIRED HINTS ../NeoPZ_install ../neopz_install)

# Mac OS X specific code
if (APPLE)
  add_definitions(-DMACOSX)
endif()

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

add_library(WANN_lib SHARED "")
target_link_libraries(WANN_lib PUBLIC NeoPZ::pz)
target_include_directories(WANN_lib PUBLIC ${PZ_INCLUDE_DIRS})
set_target_properties(WANN_lib PROPERTIES LINKER_LANGUAGE CXX)

add_subdirectory(src)
add_subdirectory(targets)
