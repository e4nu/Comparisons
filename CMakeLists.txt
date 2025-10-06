cmake_minimum_required(VERSION 3.17 FATAL_ERROR)
project(e4nuComparisons)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Add your source files
set(SOURCES
    src/Utils.h
)

# Add executable target
add_executable(comparison_1p1pi src/comparison_1p1pi.cxx ${SOURCES})

find_package(ROOT REQUIRED)
include_directories(${ROOT_INCLUDE_DIRS})
target_link_libraries(comparison_1p1pi PUBLIC ROOT::Physics)
target_link_libraries(comparison_1p1pi PUBLIC ROOT::Tree)
target_link_libraries(comparison_1p1pi PUBLIC ROOT::Hist)
target_link_libraries(comparison_1p1pi PUBLIC ROOT::Gpad)

LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_LIST_DIR}/cmake/Modules)

include(CPM)
CPMAddPackage(
  NAME NuHepMC_CPPUtils
  GIT_TAG v1-RC15
  GIT_REPOSITORY "https://github.com/NuHepMC/cpputils.git" 
  OPTIONS "NuHepMC_CPPUtils_BUILTIN_HEPMC3 ON"
  OPTIONS "HEPMC3_ENABLE_ROOTIO OFF"
  OPTIONS "CMAKE_CXX_STANDARD 17"
)
target_link_libraries(comparison_1p1pi PUBLIC NuHepMC::CPPUtils)

