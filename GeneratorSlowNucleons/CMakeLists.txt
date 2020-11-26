cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(GeneratorSlowNucleons)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS EG)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.)

set(HEADERS GeneratorSlowNucleons.h SlowNucleonModel.h SlowNucleonModelExp.h)

ROOT_GENERATE_DICTIONARY(G__GeneratorSlowNucleons ${HEADERS} LINKDEF GeneratorSlowNucleonsLinkDef.h)

#---Create a shared library with geneated dictionary
add_library(GeneratorSlowNucleons SHARED GeneratorSlowNucleons.cxx SlowNucleonModel.cxx SlowNucleonModelExp.cxx G__GeneratorSlowNucleons.cxx)
target_link_libraries(GeneratorSlowNucleons ${ROOT_LIBRARIES})


set_target_properties(GeneratorSlowNucleons
  PROPERTIES
  PUBLIC_HEADER "${HEADERS}" )


install(TARGETS GeneratorSlowNucleons
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libGeneratorSlowNucleons_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libGeneratorSlowNucleons.rootmap
      DESTINATION lib)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(GeneratorSlowNucleons PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)
