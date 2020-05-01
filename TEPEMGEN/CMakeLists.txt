cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(TEPEMGEN C CXX Fortran)

# You need to tell CMake where to find the ROOT installation. This can be done in a number of ways:
#   - ROOT built with classic configure/make use the provided $ROOTSYS/etc/cmake/FindROOT.cmake
#   - ROOT built with CMake. Add in CMAKE_PREFIX_PATH the installation prefix for ROOT
list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

#---Locate the ROOT package and defines a number of variables (e.g. ROOT_INCLUDE_DIRS)
find_package(ROOT REQUIRED COMPONENTS EG)

#---Define useful ROOT functions and macros (e.g. ROOT_GENERATE_DICTIONARY)
include(${ROOT_USE_FILE})

link_directories($ENV{PYTHIA6_ROOT}/lib)

include_directories(${CMAKE_CURRENT_SOURCE_DIR}/.)

set(HEADERS 
  TGenQEDBg.h TGenEpEmv1.h TEpEmGen.h TEcommon.h)

ROOT_GENERATE_DICTIONARY(G__TEPEMGEN ${HEADERS} LINKDEF TEPEMGENLinkDef.h)

#---Create a shared library with geneated dictionary
add_library(TEPEMGEN SHARED TGenQEDBg.cxx TGenEpEmv1.cxx TEpEmGen.cxx epemgen.f diffcross.f dtrint.f G__TEPEMGEN.cxx)
target_link_libraries(TEPEMGEN ${ROOT_LIBRARIES} pythia6 MICROCERN)


set_target_properties(TEPEMGEN
  PROPERTIES
  PUBLIC_HEADER "${HEADERS}" )


install(TARGETS TEPEMGEN
        LIBRARY DESTINATION lib
        PUBLIC_HEADER DESTINATION include)

if (${ROOT_VERSION} VERSION_GREATER "6.0")
   install(
      FILES
      ${CMAKE_CURRENT_BINARY_DIR}/libTEPEMGEN_rdict.pcm
      ${CMAKE_CURRENT_BINARY_DIR}/libTEPEMGEN.rootmap
      DESTINATION lib)
endif (${ROOT_VERSION} VERSION_GREATER "6.0")

if(${CMAKE_SYSTEM} MATCHES Darwin)
    set_target_properties(TEPEMGEN PROPERTIES LINK_FLAGS "-undefined dynamic_lookup")
endif(${CMAKE_SYSTEM} MATCHES Darwin)
