cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(GeneratorSpectators)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS EG)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.)

set(HEADERS GeneratorSpectators.h)

ROOT_GENERATE_DICTIONARY(G__GeneratorSpectators ${HEADERS} LINKDEF GeneratorSpectatorsLinkDef.h)

#---Create a shared library with generated dictionary
add_library(GeneratorSpectators SHARED GeneratorSpectators.cxx G__GeneratorSpectators.cxx)
target_link_libraries(GeneratorSpectators ${ROOT_LIBRARIES})


set_target_properties(GeneratorSpectators
  PROPERTIES
  PUBLIC_HEADER "${HEADERS}" )


install(TARGETS GeneratorSpectators
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libGeneratorSpectators_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libGeneratorSpectators.rootmap
      DESTINATION lib)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(GeneratorSpectators PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)
